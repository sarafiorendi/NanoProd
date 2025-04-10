// Take two leptons and fit them to a vertex
// Produces a collection of di-lepton vertices
// spin-off of BParkingNano framework and V0Fitter
// #include "FWCore/Framework/interface/stream/EDProducer.h"
// #include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "TrackingTools/DetLayers/interface/DetLayer.h"
// #include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
// #include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

// to create a root out file for testing
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TROOT.h"


#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

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

// const ReferenceCountingPointer<BoundCylinder> MuonTrajectoryPropagator::theBarrel_ = initBarrel();


class MuonTrajectoryPropagator : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  explicit MuonTrajectoryPropagator(const edm::ParameterSet &cfg):
    src_{consumes<std::vector<reco::GenParticle>>( cfg.getParameter<edm::InputTag>("src") )},
    recoSrc_{consumes<std::vector<pat::Muon>>( cfg.getParameter<edm::InputTag>("recoSrc") )},
    beamSpot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
    idealMagneticFieldRecordToken_(esConsumes()),
    st1propSetup_(cfg.getParameter<edm::ParameterSet>("muPropagator1st"), consumesCollector()),
    st2propSetup_(cfg.getParameter<edm::ParameterSet>("muPropagator2nd"), consumesCollector()),
    genSt2propSetup_(cfg.getParameter<edm::ParameterSet>("genMuPropagator2nd"), consumesCollector())
  {
    usesResource("TFileService");
    edm::Service<TFileService> fs_;
    outTree_ = fs_->make<TTree>("muonTree", "muonTree");
    outTree_->Branch("eventNumber", &eventNumber_, "eventNumber/L");
    outTree_->Branch("eta_mb2", &eta_mb2_, "eta_mb2/F");
    outTree_->Branch("phi_mb2", &phi_mb2_, "phi_mb2/F");
    outTree_->Branch("eta_mb1", &eta_mb1_, "eta_mb1/F");
    outTree_->Branch("phi_mb1", &phi_mb1_, "phi_mb1/F");
    outTree_->Branch("eta_ecal", &eta_ecal_, "eta_ecal/F");
    outTree_->Branch("phi_ecal", &phi_ecal_, "phi_ecal/F");
    outTree_->Branch("eta", &eta_, "eta/F");
    outTree_->Branch("phi", &phi_, "phi/F");
    outTree_->Branch("pt", &pt_, "pt/F");
    outTree_->Branch("d0", &d0_, "phi/F");
    outTree_->Branch("IP3D", &ip3d_, "phi/F");
    outTree_->Branch("isSTA", &isSTA_, "isSTA/I");
    outTree_->Branch("x_ecal", &x_ecal_, "x_ecal/F");
    outTree_->Branch("y_ecal", &y_ecal_, "y_ecal/F");
    outTree_->Branch("x_mb1",  &x_mb1_,  "x_mb1/F");
    outTree_->Branch("y_mb1",  &y_mb1_,  "y_mb1/F");
    outTree_->Branch("x_mb2",  &x_mb2_,  "x_mb2/F");
    outTree_->Branch("y_mb2",  &y_mb2_,  "y_mb2/F");


    outGENTree_ = fs_->make<TTree>("genMuonTree", "genMuonTree");
    outGENTree_->Branch("eventNumber", &eventNumber_, "eventNumber/L");
    outGENTree_->Branch("gen_eta_mb2", &gen_eta_mb2_, "gen_eta_mb2/F");
    outGENTree_->Branch("gen_phi_mb2", &gen_phi_mb2_, "gen_phi_mb2/F");
    outGENTree_->Branch("gen_eta_ecal", &gen_eta_ecal_, "gen_eta_ecal/F");
    outGENTree_->Branch("gen_phi_ecal", &gen_phi_ecal_, "gen_phi_ecal/F");
    outGENTree_->Branch("gen_eta", &gen_eta_, "eta/F");
    outGENTree_->Branch("gen_phi", &gen_phi_, "phi/F");
    outGENTree_->Branch("gen_pt", &gen_pt_, "gen_pt/F");
    outGENTree_->Branch("gen_d0", &gen_d0_, "gen_d0/F");
    outGENTree_->Branch("gen_lxy", &gen_lxy_, "gen_lxy/F");
    
    theBarrel_ = initBarrel();


//     outGENTree_->Branch("IP3D", &ip3d_, "phi/F");
//     outGENTree_->Branch("isSTA", &isSTA_, "isSTA/I");
    
//       produces<pat::CompositeCandidateCollection>();
  }

  ~MuonTrajectoryPropagator() override {}

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

//   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);


private:

  const edm::EDGetTokenT<std::vector<reco::GenParticle>> src_;
  const edm::EDGetTokenT<std::vector<pat::Muon>> recoSrc_;
  
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;

//   edm::ESHandle<MuonDetLayerGeometry> muonGeometry_;
  // to use PropagateToMuon class
  // https://github.com/cms-sw/cmssw/blob/8f2b35471f5cc4641db661c87eb6535195554380/MuonAnalysis/MuonAssociators/src/PropagateToMuon.cc
  const PropagateToMuonSetup st1propSetup_;
  const PropagateToMuonSetup st2propSetup_;
  const PropagateToMuonSetup genSt2propSetup_;
  
  PropagatorWithMaterial* forwardPropagatorECAL_;

  static ReferenceCountingPointer<BoundCylinder> theBarrel_;

  static const BoundCylinder& barrel() { return *theBarrel_; }


//   edm::ESGetToken<MuonDetLayerGeometry, MuonRecoGeometryRecord> muonGeomToken_;

  TTree* outTree_;

  long unsigned int eventNumber_;

  float eta_mb2_;
  float phi_mb2_;
  float eta_mb1_;
  float phi_mb1_;
  float eta_ecal_;
  float phi_ecal_;
  float eta_;
  float phi_;
  float pt_;
  float d0_;
  float ip3d_;
  int isSTA_;
  float x_ecal_; 
  float y_ecal_; 
  float x_mb1_; 
  float y_mb1_; 
  float x_mb2_; 
  float y_mb2_; 

  TTree* outGENTree_;
  float gen_eta_mb2_;
  float gen_phi_mb2_;
  float gen_eta_ecal_;
  float gen_phi_ecal_;
  float gen_eta_;
  float gen_phi_;
  float gen_pt_;
  float gen_d0_;
  float gen_ip3d_;
  float gen_lxy_;

};

void MuonTrajectoryPropagator::analyze(const edm::Event &evt, const edm::EventSetup & iSetup) {

  eventNumber_ = evt.id().event();

  // Get magnetic field
  edm::ESHandle<MagneticField> bFieldHandle;
  bFieldHandle = iSetup.getHandle(idealMagneticFieldRecordToken_);
  const MagneticField* bField = bFieldHandle.product();

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  evt.getByToken(beamSpot_, theBeamSpotHandle);
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
  math::XYZPoint bsPosition(theBeamSpot->position());

  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  evt.getByToken(src_, genParticles);

  forwardPropagatorECAL_ = new PropagatorWithMaterial(alongMomentum, 0.1057, bField); // change mass to muon mass

  auto const st1prop = st1propSetup_.init(iSetup);
  auto const st2prop = st2propSetup_.init(iSetup);
  auto const genSt2prop = genSt2propSetup_.init(iSetup);

  for (const auto& genParticle : *genParticles) {
    int pdgId = genParticle.pdgId();
    int status = genParticle.status();
    if (abs(pdgId) != 13 ||  status != 1 ) continue;
    
    math::XYZTLorentzVector p4 = genParticle.p4();
    math::XYZPoint vertex = genParticle.vertex();
    GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z()); 
    int charge = genParticle.charge();

//     std::cout << "Particle PDG ID: " << pdgId << "\n";
//     std::cout << "Status: " << status << "\n";
//     std::cout << "Momentum (px, py, pz, E): " 
//                   << p4.px() << ", " 
//                   << p4.py() << ", " 
//                   << p4.pz() << ", " 
//                   << p4.energy() << "\n";
    
    
    double bField_at_vtx = bField->inTesla(vtxPos).z(); // mag field strength at the vertex location
    double rinv = p4.pt() / (0.3 * charge * bField_at_vtx); // inv radius of curvature (GeV, Tesla, cm)
//     std::cout << "bField_at_vtx: " << bField_at_vtx 
//               << ", rinv: " << rinv << std::endl;
    

    FreeTrajectoryState initialState(vtxPos,
                                     GlobalVector(p4.px(), p4.py(), p4.pz()),
                                     charge, bField);
    TrajectoryStateOnSurface genStateAtMB2 = genSt2prop.extrapolate(initialState);
    if (!genStateAtMB2.isValid())
      continue;

    double etaGen = genStateAtMB2.globalPosition().eta();
    double phiGen = genStateAtMB2.globalPosition().phi();      
    float gen_d0 = (genParticle.vy()*genParticle.px() - genParticle.vx()*genParticle.py())/genParticle.pt();
    float gen_lxy = sqrt(pow(genParticle.vx() - bsPosition.x(),2) + pow(genParticle.vy() - bsPosition.y(),2));

    gen_eta_mb2_ = etaGen;
    gen_phi_mb2_ = phiGen;
    gen_eta_ = genParticle.eta();
    gen_phi_ = genParticle.phi();
    gen_pt_ = genParticle.pt();
    gen_d0_ = gen_d0;
    gen_lxy_ = gen_lxy;

    double x_gen = genStateAtMB2.globalPosition().x();      
    double y_gen = genStateAtMB2.globalPosition().y();
    double arc_tan = atan2(y_gen,x_gen);
    std::cout << arc_tan << " " << gen_phi_mb2_ << std::endl;

    TrajectoryStateOnSurface genStateAtECAL_ = forwardPropagatorECAL_->propagate(initialState, barrel());
    if (!genStateAtECAL_.isValid()) continue;
 
    gen_eta_ecal_ = genStateAtECAL_.globalPosition().eta();
    gen_phi_ecal_ = genStateAtECAL_.globalPosition().phi();      
      
    outGENTree_->Fill();  
                                      
//     const auto &fts = trajectoryStateTransform::outerFreeState((initialState), bField);
//     TrajectoryStateOnSurface tsosAtVtx = extrapolator.extrapolate(fts, pos);
//     TrajectoryStateOnSurface propagatedState = propagator->propagate(initialState, targetSurface);

  }
  
  edm::Handle<std::vector<pat::Muon>> recoMuons;
  evt.getByToken(recoSrc_, recoMuons);

//   const MuonDetLayerGeometry& muonGeom = iSetup.getData(muonGeomToken_);

  reco::TrackRef track;
  for (const auto& mu : *recoMuons) {
  
    TrajectoryStateOnSurface stateAtMB2 = st2prop.extrapolate(mu);
    if (!stateAtMB2.isValid())
      continue;

    TrajectoryStateOnSurface stateAtMB1 = st1prop.extrapolate(mu);
    if (!stateAtMB1.isValid())
      continue;

    double etaTk2 = stateAtMB2.globalPosition().eta();
    double phiTk2 = stateAtMB2.globalPosition().phi();      

    double etaTk1 = stateAtMB1.globalPosition().eta();
    double phiTk1 = stateAtMB1.globalPosition().phi();      

    eta_mb2_ = etaTk2;
    phi_mb2_ = phiTk2;
    eta_mb1_ = etaTk1;
    phi_mb1_ = phiTk1;
    eta_ =  mu.eta();
    phi_ =  mu.phi();
    pt_ =  mu.pt();
    isSTA_ =  mu.isStandAloneMuon();
    // IP from https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/interface/Muon.h#L237
    d0_ = mu.dB(pat::Muon::BS2D);
    ip3d_= mu.dB(pat::Muon::BS3D);

    double x_reco = stateAtMB2.globalPosition().x();      
    double y_reco = stateAtMB2.globalPosition().y();
    double arc_tan = atan2(y_reco,x_reco);
//     std::cout << "reco: " << arc_tan << " " << phi_mb2_ << std::endl;

    double px_reco = stateAtMB2.globalMomentum().x();      
    double py_reco = stateAtMB2.globalMomentum().y();
    double px_arc_tan = atan2(py_reco,px_reco);
//     std::cout << "pxpy: " << px_arc_tan << std::endl;

    x_mb1_ = stateAtMB1.globalPosition().x();      
    y_mb1_ = stateAtMB1.globalPosition().y();      
    x_mb2_ = stateAtMB2.globalPosition().x();      
    y_mb2_ = stateAtMB2.globalPosition().y();      


    // at ECAL surface
    bool noTrack = false;
    TString type = "";
    if (mu.isGlobalMuon()) {track = mu.globalTrack(); type = "global";}
    else if (mu.isStandAloneMuon()) {track = mu.standAloneMuon(); type = "standalone";}
    else if (mu.isTrackerMuon()) {track = mu.innerTrack(); type = "tracker";}
    else noTrack = true;
    
    if (!noTrack){
      
      const reco::TransientTrack muTransientTrack(track, &(*bFieldHandle));
      if (!muTransientTrack.isValid()) continue;
      FreeTrajectoryState innerMuTSOS = muTransientTrack.initialFreeState();

      TrajectoryStateOnSurface stateAtECAL_ = forwardPropagatorECAL_->propagate(innerMuTSOS, barrel());
      if (!stateAtECAL_.isValid()) continue;

      eta_ecal_ = stateAtECAL_.globalPosition().eta();
      phi_ecal_ = stateAtECAL_.globalPosition().phi();      

      x_ecal_ = stateAtECAL_.globalPosition().x();      
      y_ecal_ = stateAtECAL_.globalPosition().y();      
    }

    outTree_->Fill();  
    

    //Extrapolate to layer 4
//     const DetLayer *dtLay = muonGeom.allDTLayers()[3];
//     barrelCylinder_ = dynamic_cast<const BoundCylinder *>(&dtLay->surface());
//     barrelHalfLength_ = barrelCylinder_->bounds().length() / 2;

//     if (!noTrack && mu.pt() > 10)
//     {
//         FreeTrajectoryState start;
//         start = trajectoryStateTransform::outerFreeState(*track, bField);
//         TrajectoryStateOnSurface tsos = prop.propagate(start, *barrelCylinder_);
//     }  
  }
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void MuonTrajectoryPropagator::beginJob() {
}  

// ------------ method called once each job just after ending the event loop  ------------
void MuonTrajectoryPropagator::endJob() {
  outTree_->GetDirectory()->cd();
  outTree_->Write();
  outGENTree_->GetDirectory()->cd();
  outGENTree_->Write();
}

ReferenceCountingPointer<BoundCylinder> MuonTrajectoryPropagator::theBarrel_ = nullptr;


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonTrajectoryPropagator);
