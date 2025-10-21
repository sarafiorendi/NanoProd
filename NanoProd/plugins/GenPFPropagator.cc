/*
 * Add propagated coords to GEN PF Cands
 */

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

#include <iostream>
#include <fstream>
#include <cstring>


constexpr float epsilon = 0.001;
/** Hard-wired numbers defining the surfaces on which the crystal front faces lie. */
// https://github.com/cms-sw/cmssw/blob/6cad17bc389ac53d9b2193d99f5a76713a7a99f4/RecoEgamma/EgammaPhotonAlgos/src/ConversionTrackEcalImpactPoint.cc#L73-L82
constexpr float barrelRadius = 129.f;       // p81, p50, ECAL TDR
constexpr float barrelHalfLength = 270.9f;  // p81, p50, ECAL TDR
constexpr float endcapRadius = 171.1f;      // fig 3.26, p81, ECAL TDR
constexpr float endcapZ = 320.5f;           // fig 3.26, p81, ECAL TDR


static BoundCylinder* initBarrel() {
  Surface::RotationType rot;  // unit rotation matrix

  return new Cylinder(
      barrelRadius,
      Surface::PositionType(0, 0, 0),
      rot,
      new SimpleCylinderBounds(barrelRadius - epsilon, barrelRadius + epsilon, -barrelHalfLength, barrelHalfLength));
}

static BoundDisk* initNegativeEcalEndcap() {
  Surface::RotationType rot;  // unit rotation matrix
  return new BoundDisk(Surface::PositionType(0, 0, -endcapZ), 
                       rot, 
                       new SimpleDiskBounds(0, endcapRadius, -epsilon, epsilon));
}      

static BoundDisk* initPositiveEcalEndcap() {
  Surface::RotationType rot;  // unit rotation matrix
  return new BoundDisk(Surface::PositionType(0, 0, endcapZ), 
                       rot, 
                       new SimpleDiskBounds(0, endcapRadius, -epsilon, epsilon));
}      


class GenPFPropagator : public edm::stream::EDProducer<> {
public:
    explicit GenPFPropagator(const edm::ParameterSet&);
    ~GenPFPropagator(){};

private:
    PropagatorWithMaterial* forwardPropagatorECALpion_;//, *forwardPropagatorECALpion_, *forwardPropagatorECALpionRK_;

    static ReferenceCountingPointer<BoundCylinder> theBarrel_;
    static ReferenceCountingPointer<BoundDisk> thePositiveEndcap_;
    static ReferenceCountingPointer<BoundDisk> theNegativeEndcap_;

    static const BoundCylinder& barrel() { return *theBarrel_; }
    static const BoundDisk& diskPlus() { return *thePositiveEndcap_; }
    static const BoundDisk& diskMinus() { return *theNegativeEndcap_; }

    void produce(edm::Event&, const edm::EventSetup&) override;

    const edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartToken_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;
};


GenPFPropagator::GenPFPropagator(const edm::ParameterSet& cfg)
    : 
      genPartToken_{consumes<std::vector<pat::PackedGenParticle>>( cfg.getParameter<edm::InputTag>("src") )},
      idealMagneticFieldRecordToken_(esConsumes())
{  
    theBarrel_ = initBarrel();
    thePositiveEndcap_ = initPositiveEcalEndcap();
    theNegativeEndcap_ = initNegativeEcalEndcap();
 
    produces<edm::ValueMap<float>>("etaAtEcal"); 
    produces<edm::ValueMap<float>>("phiAtEcal"); 
}

void GenPFPropagator::produce(edm::Event& event, const edm::EventSetup& setup) {

    edm::Handle<std::vector<pat::PackedGenParticle>> genParticles;
    event.getByToken(genPartToken_, genParticles);
   
    // Get magnetic field
    edm::ESHandle<MagneticField> bFieldHandle;
    bFieldHandle = setup.getHandle(idealMagneticFieldRecordToken_);
    const MagneticField* bField = bFieldHandle.product();

    forwardPropagatorECALpion_ = new PropagatorWithMaterial(alongMomentum, 0.139 , bField, 6, false, -1, true); // pion mass

    const size_t genPart_size = genParticles->size();
    std::vector <Float_t> v_eta_ecal(genPart_size, -99);
    std::vector <Float_t> v_phi_ecal(genPart_size, -99);

    float gen_eta_ecal = -99;
    float gen_phi_ecal = -99;
    float gen_charge_ = -99;
    
    std::cout << "[GenPFPropagator]" << std::endl;
    for(size_t igen = 0; igen < genParticles->size();igen++){
      
        // reset output variables
        gen_eta_ecal = -99;
        gen_phi_ecal = -99;

        const pat::PackedGenParticle genPF = genParticles->at(igen);
        int pdgId = genPF.pdgId();
        int status = genPF.status();
        if (abs(pdgId) != 211 ) continue;
//         if (abs(pdgId) != 211 || status != 1) continue;
        std::cout << "\t [GenPFPropagator] pion with pT " << genPF.pt() << std::endl;
  
        math::XYZTLorentzVector p4 = genPF.p4();
        math::XYZPoint vertex = genPF.vertex();
        GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z()); 
        gen_charge_ = genPF.charge();
        FreeTrajectoryState initialState(vtxPos,
                                         GlobalVector(p4.px(), p4.py(), p4.pz()),
                                         gen_charge_, bField);
      
        TrajectoryStateOnSurface genStateAtECAL_ = forwardPropagatorECALpion_->propagate(initialState, barrel());
        if (!genStateAtECAL_.isValid() || (genStateAtECAL_.isValid() && fabs(genStateAtECAL_.globalPosition().eta()) > 1.479f)) {
           if (genPF.eta() > 0.) {
             genStateAtECAL_ = forwardPropagatorECALpion_->propagate(initialState, diskPlus());
           } else {
             genStateAtECAL_ = forwardPropagatorECALpion_->propagate(initialState, diskMinus());
           }
        }              
        if (genStateAtECAL_.isValid()) {
          gen_eta_ecal = genStateAtECAL_.globalPosition().eta();
          gen_phi_ecal = genStateAtECAL_.globalPosition().phi();      
        }      
            
        v_eta_ecal.at(igen) = gen_eta_ecal;
        v_phi_ecal.at(igen) = gen_phi_ecal;
    } // end loop on GenPart  

    // now save muon coordinates at ECAL surface
    std::unique_ptr<edm::ValueMap<float>> vm_eta(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_eta(*vm_eta);
    filler_eta.insert(genParticles, v_eta_ecal.begin(), v_eta_ecal.end());
    filler_eta.fill();
    event.put(std::move(vm_eta), "etaAtEcal"); 

    std::unique_ptr<edm::ValueMap<float>> vm_phi(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_phi(*vm_phi);
    filler_phi.insert(genParticles, v_phi_ecal.begin(), v_phi_ecal.end());
    filler_phi.fill();
    event.put(std::move(vm_phi), "phiAtEcal"); 
}


ReferenceCountingPointer<BoundCylinder> GenPFPropagator::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenPFPropagator::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenPFPropagator::theNegativeEndcap_ = nullptr;

DEFINE_FWK_MODULE(GenPFPropagator);
