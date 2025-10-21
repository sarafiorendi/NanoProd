/*
 * Add propagated coords to GEN level muons, only if coming from a stau -> tau decay
 */

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

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


class GenMuonPropagator : public edm::stream::EDProducer<> {
public:
    explicit GenMuonPropagator(const edm::ParameterSet&);
    ~GenMuonPropagator(){};
    bool isAncestor(const reco::Candidate*, const reco::Candidate *);

private:
    PropagatorWithMaterial* forwardPropagatorECAL_, *forwardPropagatorECALpion_, *forwardPropagatorECALpionRK_;

    static ReferenceCountingPointer<BoundCylinder> theBarrel_;
    static ReferenceCountingPointer<BoundDisk> thePositiveEndcap_;
    static ReferenceCountingPointer<BoundDisk> theNegativeEndcap_;

    static const BoundCylinder& barrel() { return *theBarrel_; }
    static const BoundDisk& diskPlus() { return *thePositiveEndcap_; }
    static const BoundDisk& diskMinus() { return *theNegativeEndcap_; }

    void produce(edm::Event&, const edm::EventSetup&) override;

    const edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;
};


GenMuonPropagator::GenMuonPropagator(const edm::ParameterSet& cfg)
    : 
      genPartToken_{consumes<std::vector<reco::GenParticle>>( cfg.getParameter<edm::InputTag>("src") )},
      idealMagneticFieldRecordToken_(esConsumes())
{  
    theBarrel_ = initBarrel();
    thePositiveEndcap_ = initPositiveEcalEndcap();
    theNegativeEndcap_ = initNegativeEcalEndcap();
 
    produces<edm::ValueMap<float>>("etaAtEcal"); 
    produces<edm::ValueMap<float>>("phiAtEcal"); 
}

void GenMuonPropagator::produce(edm::Event& event, const edm::EventSetup& setup) {

    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    event.getByToken(genPartToken_, genParticles);
   
    // Get magnetic field
    edm::ESHandle<MagneticField> bFieldHandle;
    bFieldHandle = setup.getHandle(idealMagneticFieldRecordToken_);
    const MagneticField* bField = bFieldHandle.product();

    forwardPropagatorECAL_ = new PropagatorWithMaterial(alongMomentum, 0.1057, bField); // muon mass
    forwardPropagatorECALpion_ = new PropagatorWithMaterial(alongMomentum, 0.139 , bField, 6, false, -1, true); // pion mass
    forwardPropagatorECALpionRK_ = new PropagatorWithMaterial(alongMomentum, 0.139 , bField, 6, true, -1, true); // pion mass

    const size_t genPart_size = genParticles->size();
    std::vector <Float_t> v_eta_ecal(genPart_size, -99);
    std::vector <Float_t> v_phi_ecal(genPart_size, -99);

    float gen_eta_ecal = -99;
    float gen_phi_ecal = -99;
    float gen_charge_ = -99;
    
    for(size_t igen = 0; igen < genParticles->size();igen++){
      const reco::Candidate * genStau = &(*genParticles)[igen];
      int stau_pdgId = genStau->pdgId();
      if (abs(stau_pdgId) != 1000015 || genStau->status() != 22 ) continue;
      
      for (size_t jgen = 0; jgen < genParticles->size();jgen++) {
      
        // reset output variables
        gen_eta_ecal = -99;
        gen_phi_ecal = -99;

        const reco::GenParticle genMu = genParticles->at(jgen);
        int pdgId = genMu.pdgId();
        int status = genMu.status();
        if (abs(pdgId) != 13 || status != 1) continue;
  
        const reco::Candidate * motherInPrunedCollection = genMu.mother(0) ;
        if(motherInPrunedCollection != nullptr && isAncestor( genStau , motherInPrunedCollection)) {
  
          math::XYZTLorentzVector p4 = genMu.p4();
          math::XYZPoint vertex = genMu.vertex();
          GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z()); 
          gen_charge_ = genMu.charge();
          FreeTrajectoryState initialState(vtxPos,
                                           GlobalVector(p4.px(), p4.py(), p4.pz()),
                                           gen_charge_, bField);
      
          TrajectoryStateOnSurface genStateAtECAL_ = forwardPropagatorECAL_->propagate(initialState, barrel());
          if (!genStateAtECAL_.isValid() || (genStateAtECAL_.isValid() && fabs(genStateAtECAL_.globalPosition().eta()) > 1.479f)) {
             if (genMu.eta() > 0.) {
               genStateAtECAL_ = forwardPropagatorECAL_->propagate(initialState, diskPlus());
             } else {
               genStateAtECAL_ = forwardPropagatorECAL_->propagate(initialState, diskMinus());
             }
          }              
          if (genStateAtECAL_.isValid()) {
            gen_eta_ecal = genStateAtECAL_.globalPosition().eta();
            gen_phi_ecal = genStateAtECAL_.globalPosition().phi();      
          }      
            
        } // end if muon is from stau
        v_eta_ecal.at(jgen) = gen_eta_ecal;
        v_phi_ecal.at(jgen) = gen_phi_ecal;
      } // end loop on gen part = muons
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

// -----------------------
bool GenMuonPropagator::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
  {
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}


ReferenceCountingPointer<BoundCylinder> GenMuonPropagator::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenMuonPropagator::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenMuonPropagator::theNegativeEndcap_ = nullptr;

DEFINE_FWK_MODULE(GenMuonPropagator);
