/*
 * Add propagated coords to GEN level muons, only if coming from a stau -> tau decay or if running on a cosmics MC sample
 */

#include <memory>

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
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"


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

private:

    static ReferenceCountingPointer<BoundCylinder> theBarrel_;
    static ReferenceCountingPointer<BoundDisk> thePositiveEndcap_;
    static ReferenceCountingPointer<BoundDisk> theNegativeEndcap_;

    static const BoundCylinder& barrel() { return *theBarrel_; }
    static const BoundDisk& diskPlus() { return *thePositiveEndcap_; }
    static const BoundDisk& diskMinus() { return *theNegativeEndcap_; }

    void produce(edm::Event&, const edm::EventSetup&) override;

    const edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_;
    const PropagateToMuonSetup genSt2propSetup_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;
};


GenMuonPropagator::GenMuonPropagator(const edm::ParameterSet& cfg)
    : 
      genPartToken_{consumes<std::vector<reco::GenParticle>>( cfg.getParameter<edm::InputTag>("src") )},
      genSt2propSetup_(cfg.getParameter<edm::ParameterSet>("genMuPropagator2nd"), consumesCollector()),
      idealMagneticFieldRecordToken_(esConsumes())
{  
    if (!theBarrel_){
      theBarrel_ = initBarrel();
      thePositiveEndcap_ = initPositiveEcalEndcap();
      theNegativeEndcap_ = initNegativeEcalEndcap();
    }  
 
    produces<edm::ValueMap<float>>("etaAtEcal"); 
    produces<edm::ValueMap<float>>("phiAtEcal"); 
    produces<edm::ValueMap<float>>("etaAtMB2"); 
    produces<edm::ValueMap<float>>("phiAtMB2"); 
}

void GenMuonPropagator::produce(edm::Event& event, const edm::EventSetup& setup) {

    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    event.getByToken(genPartToken_, genParticles);
   
    // Get magnetic field
    edm::ESHandle<MagneticField> bFieldHandle;
    bFieldHandle = setup.getHandle(idealMagneticFieldRecordToken_);
    const MagneticField* bField = bFieldHandle.product();
    
    PropagatorWithMaterial forwardPropagatorECAL(alongMomentum, 0.1057, bField);

    auto const genSt2prop = genSt2propSetup_.init(setup);

    const size_t genPart_size = genParticles->size();
    std::vector <Float_t> v_eta_ecal(genPart_size, -99);
    std::vector <Float_t> v_phi_ecal(genPart_size, -99);
    std::vector <Float_t> v_eta_mb2(genPart_size, -99);
    std::vector <Float_t> v_phi_mb2(genPart_size, -99);

    for (size_t jgen = 0; jgen < genPart_size; jgen++) {
      
      const reco::GenParticle& genMu = genParticles->at(jgen);
      int pdgId = genMu.pdgId();
      int status = genMu.status();
      if (abs(pdgId) != 13 || status != 1) continue;
  
      math::XYZTLorentzVector p4 = genMu.p4();
      math::XYZPoint vertex = genMu.vertex();
      GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z()); 
      FreeTrajectoryState initialState(vtxPos,
                                       GlobalVector(p4.px(), p4.py(), p4.pz()),
                                       genMu.charge(), bField);
      
      TrajectoryStateOnSurface genStateAtECAL_ = forwardPropagatorECAL.propagate(initialState, barrel());
      if (!genStateAtECAL_.isValid() || fabs(genStateAtECAL_.globalPosition().eta()) > 1.479f) {
        if (genMu.eta() > 0.) {
          genStateAtECAL_ = forwardPropagatorECAL.propagate(initialState, diskPlus());
        } else {
          genStateAtECAL_ = forwardPropagatorECAL.propagate(initialState, diskMinus());
        }
      }              
      if (genStateAtECAL_.isValid()) {
        v_eta_ecal[jgen] = genStateAtECAL_.globalPosition().eta();
        v_phi_ecal[jgen] = genStateAtECAL_.globalPosition().phi();      
      }

      TrajectoryStateOnSurface genStateAtMB2 = genSt2prop.extrapolate(initialState);
      if (genStateAtMB2.isValid()){
        v_eta_mb2[jgen] = genStateAtMB2.globalPosition().eta();
        v_phi_mb2[jgen] = genStateAtMB2.globalPosition().phi(); 
      }
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

    // now save muon coordinates at MB2 surface
    std::unique_ptr<edm::ValueMap<float>> vm_etamb2(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_etamb2(*vm_etamb2);
    filler_etamb2.insert(genParticles, v_eta_mb2.begin(), v_eta_mb2.end());
    filler_etamb2.fill();
    event.put(std::move(vm_etamb2), "etaAtMB2"); 

    std::unique_ptr<edm::ValueMap<float>> vm_phimb2(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_phimb2(*vm_phimb2);
    filler_phimb2.insert(genParticles, v_phi_mb2.begin(), v_phi_mb2.end());
    filler_phimb2.fill();
    event.put(std::move(vm_phimb2), "phiAtMB2"); 
}

// -----------------------

ReferenceCountingPointer<BoundCylinder> GenMuonPropagator::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenMuonPropagator::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenMuonPropagator::theNegativeEndcap_ = nullptr;

DEFINE_FWK_MODULE(GenMuonPropagator);
