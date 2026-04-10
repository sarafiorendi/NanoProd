/*
 * Add coordinates of PF candidates once their trajectory is propagated to the ECAL surface
 */

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


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


class PFCandidatePropagator : public edm::stream::EDProducer<> {
public:
    explicit PFCandidatePropagator(const edm::ParameterSet&);
    ~PFCandidatePropagator(){};

private:
    PropagatorWithMaterial *forwardPropagatorECALpion_;

    static ReferenceCountingPointer<BoundCylinder> theBarrel_;
    static ReferenceCountingPointer<BoundDisk> thePositiveEndcap_;
    static ReferenceCountingPointer<BoundDisk> theNegativeEndcap_;

    static const BoundCylinder& barrel() { return *theBarrel_; }
    static const BoundDisk& diskPlus() { return *thePositiveEndcap_; }
    static const BoundDisk& diskMinus() { return *theNegativeEndcap_; }

    void produce(edm::Event&, const edm::EventSetup&) override;

    const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;
    
};


PFCandidatePropagator::PFCandidatePropagator(const edm::ParameterSet& cfg)
    : pfToken_{consumes<pat::PackedCandidateCollection>( cfg.getParameter<edm::InputTag>("src") )},    
      idealMagneticFieldRecordToken_(esConsumes())
{  
    theBarrel_ = initBarrel();
    thePositiveEndcap_ = initPositiveEcalEndcap();
    theNegativeEndcap_ = initNegativeEcalEndcap();
 
    produces<edm::ValueMap<float>>("etaAtEcal"); 
    produces<edm::ValueMap<float>>("phiAtEcal"); 
}

void PFCandidatePropagator::produce(edm::Event& event, const edm::EventSetup& setup) {

    edm::Handle<pat::PackedCandidateCollection> pfCands;
    event.getByToken(pfToken_, pfCands);
    if (!pfCands.isValid()) {
        std::cout << "pfCands not valid" << std::endl;
        return;
    }

    // Get magnetic field
    edm::ESHandle<MagneticField> bFieldHandle;
    bFieldHandle = setup.getHandle(idealMagneticFieldRecordToken_);
    const MagneticField* bField = bFieldHandle.product();

    forwardPropagatorECALpion_ = new PropagatorWithMaterial(alongMomentum, 0.139 , bField, 6, false, -1, true); // pion mass
    
    const size_t pf_size = pfCands->size();
    std::vector <Float_t> v_eta_ecal(pf_size, -99);
    std::vector <Float_t> v_phi_ecal(pf_size, -99);
    
    float pf_eta_at_ecal = -99;
    float pf_phi_at_ecal = -99;
    
    // loop on the muons   
    for(size_t pfIndex = 0; pfIndex < pf_size; ++pfIndex)
    {
        // reset outputs
        pf_eta_at_ecal = -99;
        pf_phi_at_ecal = -99;

        const auto& the_pf = pfCands->at(pfIndex);
        if (abs(the_pf.pdgId()) != 211) continue;
        if (!the_pf.hasTrackDetails()) continue;
      
        reco::Track pfTrack = the_pf.pseudoTrack() ;
        reco::TransientTrack pfTransientTrack(pfTrack, &(*bFieldHandle));
        if (pfTransientTrack.isValid()) {
            FreeTrajectoryState pfTSOS = pfTransientTrack.initialFreeState();
            TrajectoryStateOnSurface stateAtECAL_ = forwardPropagatorECALpion_->propagate(pfTSOS, barrel());
            if (!stateAtECAL_.isValid() || (stateAtECAL_.isValid() && fabs(stateAtECAL_.globalPosition().eta()) > 1.479f)) {
                if (the_pf.eta() > 0.) {
                    stateAtECAL_ = forwardPropagatorECALpion_->propagate(pfTSOS, diskPlus());
                } else {
                    stateAtECAL_ = forwardPropagatorECALpion_->propagate(pfTSOS, diskMinus());
                }
            }
            if (stateAtECAL_.isValid()) {  
                float eta_ecal_ = stateAtECAL_.globalPosition().eta();
                float phi_ecal_ = stateAtECAL_.globalPosition().phi();  
                pf_eta_at_ecal = eta_ecal_;
                pf_phi_at_ecal = phi_ecal_;
            }    
     
        }
        v_eta_ecal.at(pfIndex) = pf_eta_at_ecal;
        v_phi_ecal.at(pfIndex) = pf_phi_at_ecal;
    } // end loop on pf cands 


    // now save pf coordinates at ECAL surface
    std::unique_ptr<edm::ValueMap<float>> vm_eta(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_eta(*vm_eta);
    filler_eta.insert(pfCands, v_eta_ecal.begin(), v_eta_ecal.end());
    filler_eta.fill();
    event.put(std::move(vm_eta), "etaAtEcal"); 

    std::unique_ptr<edm::ValueMap<float>> vm_phi(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_phi(*vm_phi);
    filler_phi.insert(pfCands, v_phi_ecal.begin(), v_phi_ecal.end());
    filler_phi.fill();
    event.put(std::move(vm_phi), "phiAtEcal"); 
}

ReferenceCountingPointer<BoundCylinder> PFCandidatePropagator::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> PFCandidatePropagator::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> PFCandidatePropagator::theNegativeEndcap_ = nullptr;

DEFINE_FWK_MODULE(PFCandidatePropagator);
