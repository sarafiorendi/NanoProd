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

typedef std::pair<TrajectoryStateOnSurface, double> TsosPath;

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


template<typename T>
void putValueMap(edm::Event& event,
                 const edm::Handle<std::vector<reco::GenParticle>>& src,
                 const std::vector<T>& values,
                 const std::string& label)
{
    auto vm = std::make_unique<edm::ValueMap<T>>();
    typename edm::ValueMap<T>::Filler filler(*vm);
    filler.insert(src, values.begin(), values.end());
    filler.fill();
    event.put(std::move(vm), label);
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
    edm::ESGetToken<Propagator, TrackingComponentsRecord> propAlongToken_;
    edm::ESGetToken<Propagator, TrackingComponentsRecord> propOppositeToken_;
    const edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
    
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;
};


GenMuonPropagator::GenMuonPropagator(const edm::ParameterSet& cfg)
    : 
      genPartToken_{consumes<std::vector<reco::GenParticle>>( cfg.getParameter<edm::InputTag>("src") )},
      genSt2propSetup_(cfg.getParameter<edm::ParameterSet>("genMuPropagator2nd"), consumesCollector()),
      propAlongToken_{esConsumes<Propagator, TrackingComponentsRecord>(cfg.getParameter<edm::ESInputTag>("propagatorAlong"))},
      propOppositeToken_{esConsumes<Propagator, TrackingComponentsRecord>(cfg.getParameter<edm::ESInputTag>("propagatorOpposite"))},
      muonsToken_{consumes<std::vector<pat::Muon>>(  cfg.getParameter<edm::InputTag>("reco"))},
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
    produces<edm::ValueMap<float>>("xAtMB2"); 
    produces<edm::ValueMap<float>>("yAtMB2"); 
    produces<edm::ValueMap<float>>("zAtMB2"); 
    produces<edm::ValueMap<float>>("pxAtMB2"); 
    produces<edm::ValueMap<float>>("pyAtMB2"); 
    produces<edm::ValueMap<float>>("pzAtMB2"); 
    produces<edm::ValueMap<float>>("propEtaAtMB2"); 
    produces<edm::ValueMap<float>>("propPhiAtMB2"); 
    produces<edm::ValueMap<float>>("initr"); 
    produces<edm::ValueMap<float>>("initz"); 
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
    auto const& propagatorAlong = setup.getData(propAlongToken_);
    auto const& propagatorOpposite = setup.getData(propOppositeToken_);

    const size_t genPart_size = genParticles->size();
    std::vector <Float_t> v_eta_ecal(genPart_size, -99);
    std::vector <Float_t> v_phi_ecal(genPart_size, -99);
    std::vector <Float_t> v_eta_mb2(genPart_size, -99);
    std::vector <Float_t> v_phi_mb2(genPart_size, -99);
    std::vector <Float_t> v_x_mb2(genPart_size, -9999);
    std::vector <Float_t> v_y_mb2(genPart_size, -9999);
    std::vector <Float_t> v_z_mb2(genPart_size, -9999);
    std::vector <Float_t> v_px_mb2(genPart_size, -99);
    std::vector <Float_t> v_py_mb2(genPart_size, -99);
    std::vector <Float_t> v_pz_mb2(genPart_size, -99);
    std::vector <Float_t> v_prop_eta_mb2(genPart_size, -99);
    std::vector <Float_t> v_prop_phi_mb2(genPart_size, -99);
    std::vector <Float_t> v_init_r(genPart_size, -99);
    std::vector <Float_t> v_init_z(genPart_size, -99);

    for (size_t jgen = 0; jgen < genPart_size; jgen++) {
      
      const reco::GenParticle& genMu = genParticles->at(jgen);
      int pdgId = genMu.pdgId();
      int status = genMu.status();
      if (abs(pdgId) != 13 || status != 1) continue;
  
//       std::cout << "genParticle \t\t\t\t\t\t\t\t " << genMu.eta() ;// << std::endl;
//       std::cout << "\t " << genMu.phi() ;// << std::endl;
//       std::cout << "\t " << genMu.pt() << std::endl;
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
        v_prop_eta_mb2[jgen] = genStateAtMB2.globalPosition().eta();
        v_prop_phi_mb2[jgen] = genStateAtMB2.globalPosition().phi(); 
      }
      
      // alternative way
      const Surface::RotationType dummyRot;
      float radius_cyl = 523.854; //820.;
      float z_cyl = 659; // 1500;
      Cylinder::CylinderPointer theTargetCylinder =
        Cylinder::build(Surface::PositionType(0., 0., 0.), dummyRot, radius_cyl);
//       std::cout << "initial state: 	" << initialState.position().perp() ; // << std::endl;
//       std::cout <<"\t " << initialState.position().z() ;// << std::endl;
//       std::cout <<"\t " << initialState.position().eta() ;// << std::endl;
//       std::cout <<"\t " << initialState.position().phi() ; // << std::endl;
//       std::cout <<"\t " << initialState.momentum().eta() ;// << std::endl;
//       std::cout <<"\t " << initialState.momentum().phi() << std::endl;
      v_init_r[jgen] = initialState.position().perp();
      v_init_z[jgen] = initialState.position().z(); 
  
      bool isInsideInitial = 
        initialState.position().perp() < radius_cyl && initialState.position().z() >= -z_cyl && initialState.position().z() <= z_cyl;
      bool isExternalToSurface = 
        initialState.position().perp() > radius_cyl || initialState.position().z() >= -z_cyl || initialState.position().z() <= z_cyl;
      bool isPointingTwrOrigin = 
        initialState.momentum().dot(GlobalVector(initialState.position().x(), initialState.position().y(), initialState.position().z())) < 0  ;
      
//       const Propagator* selectedPropagator = isInsideInitial ? &propagatorAlong : &propagatorAlong;
      const Propagator* selectedPropagator = (isExternalToSurface == isPointingTwrOrigin) ? &propagatorAlong : &propagatorOpposite ;
      TsosPath tsosPath = selectedPropagator->propagateWithPath(initialState, *theTargetCylinder);
      if (!tsosPath.first.isValid()) {
//         std::cout << "not valid alternative prop" << std::endl;
//           return false;
      }
      else{
//         std::cout << "tsos \t\t" <<
//           tsosPath.first.globalPosition().perp() << "\t " <<
//           tsosPath.first.globalPosition().z() << "\t " <<
//           tsosPath.first.globalPosition().eta() << "\t " <<
//           tsosPath.first.globalPosition().phi() << "\t " <<
//           tsosPath.first.globalMomentum().eta() << "\t " <<
//           tsosPath.first.globalMomentum().phi() << "\t " <<
//           std::endl;
          v_eta_mb2[jgen] = tsosPath.first.globalPosition().eta();
          v_phi_mb2[jgen] = tsosPath.first.globalPosition().phi(); 
          v_x_mb2[jgen] = tsosPath.first.globalPosition().x();
          v_y_mb2[jgen] = tsosPath.first.globalPosition().y(); 
          v_z_mb2[jgen] = tsosPath.first.globalPosition().z(); 
          v_px_mb2[jgen] = tsosPath.first.globalMomentum().x();
          v_py_mb2[jgen] = tsosPath.first.globalMomentum().y(); 
          v_pz_mb2[jgen] = tsosPath.first.globalMomentum().z(); 
      }
//       if (genStateAtMB2.isValid()){
//         std::cout << "original MB2:   " << 
//           genStateAtMB2.globalPosition().perp() << "\t " <<
//           genStateAtMB2.globalPosition().z() << "\t " <<
//           genStateAtMB2.globalPosition().eta() << "\t " <<
//           genStateAtMB2.globalPosition().phi() << "\t " <<
//           genStateAtMB2.globalMomentum().eta() << "\t " <<
//           genStateAtMB2.globalMomentum().phi() << "\t " <<
//           std::endl;
//       }
//       if (checkFinalZ) {
//           bool withinZRange = tsosPath.first.globalPosition().z() >= minZ &&
//                               tsosPath.first.globalPosition().z() <= maxZ;
//           return withinZRange;
//       }
    } // end loop on GenPart
    
    // tmp loop on reco
//     edm::Handle<std::vector<pat::Muon>> recoMuons;
//     event.getByToken(muonsToken_, recoMuons);
//     const size_t muons_size = recoMuons->size();
//     for(size_t muIndex = 0; muIndex < muons_size; ++muIndex)
//     {
//       const auto& mu = recoMuons->at(muIndex);
//         std::cout << "reco muon: \t\t\t\t\t\t\t\t   " << 
//           mu.eta() << "\t " <<
//           mu.phi() << "\t " <<
//           mu.pt() << "\t " <<
// //           mu.phi() << "\t " <<
//           std::endl;
//     }
//         std::cout << "\t\t" <<
//           std::endl;

    // now save muon coordinates at ECAL surface
    putValueMap(event, genParticles, v_eta_ecal, "etaAtEcal");
    putValueMap(event, genParticles, v_phi_ecal, "phiAtEcal");
    // now save muon coordinates at MB2 surface
    putValueMap(event, genParticles, v_eta_mb2, "etaAtMB2");
    putValueMap(event, genParticles, v_phi_mb2, "phiAtMB2");
    putValueMap(event, genParticles, v_x_mb2, "xAtMB2");
    putValueMap(event, genParticles, v_y_mb2, "yAtMB2");
    putValueMap(event, genParticles, v_z_mb2, "zAtMB2");
    putValueMap(event, genParticles, v_px_mb2, "pxAtMB2");
    putValueMap(event, genParticles, v_py_mb2, "pyAtMB2");
    putValueMap(event, genParticles, v_pz_mb2, "pzAtMB2");

    putValueMap(event, genParticles, v_prop_eta_mb2, "propEtaAtMB2");
    putValueMap(event, genParticles, v_prop_phi_mb2, "propPhiAtMB2");
    putValueMap(event, genParticles, v_init_r, "initr");
    putValueMap(event, genParticles, v_init_z, "initz");
}

// -----------------------

ReferenceCountingPointer<BoundCylinder> GenMuonPropagator::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenMuonPropagator::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenMuonPropagator::theNegativeEndcap_ = nullptr;

DEFINE_FWK_MODULE(GenMuonPropagator);
