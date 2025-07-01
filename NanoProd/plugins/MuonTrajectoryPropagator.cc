#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/TrackerBounds.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

// to create a root out file for testing
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TROOT.h"



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

static BoundSurface* initSurface() {
  Surface::RotationType rot;  // unit rotation matrix

  return new Disk(Surface::PositionType(0, 0, 0),
                  rot,
                  new SimpleDiskBounds(0, 2*TrackerBounds::radius(), -epsilon, epsilon));
}


class MuonTrajectoryPropagator : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  explicit MuonTrajectoryPropagator(const edm::ParameterSet &cfg):
    ttkToken_{esConsumes(edm::ESInputTag{"", "TransientTrackBuilder"})},
    src_{consumes<std::vector<reco::GenParticle>>( cfg.getParameter<edm::InputTag>("src") )},
    recoSrc_{consumes<std::vector<pat::Muon>>( cfg.getParameter<edm::InputTag>("recoSrc") )},
    packedGenToken_{consumes<edm::View<pat::PackedGenParticle>>( cfg.getParameter<edm::InputTag>("packedGen") )},
    beamSpot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
    tracksToken_{consumes<reco::TrackCollection>( cfg.getParameter<edm::InputTag>("tracksForIso") )},
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
    outTree_->Branch("isGlobal", &isGlobal_, "isGlobal/I");
    outTree_->Branch("isTracker", &isTracker_, "isTracker/I");
    outTree_->Branch("charge", &charge_, "charge/I");
    outTree_->Branch("x_ecal", &x_ecal_, "x_ecal/F");
    outTree_->Branch("y_ecal", &y_ecal_, "y_ecal/F");
    outTree_->Branch("x_mb1",  &x_mb1_,  "x_mb1/F");
    outTree_->Branch("y_mb1",  &y_mb1_,  "y_mb1/F");
    outTree_->Branch("x_mb2",  &x_mb2_,  "x_mb2/F");
    outTree_->Branch("y_mb2",  &y_mb2_,  "y_mb2/F");
    outTree_->Branch("x_origin",  &x_origin_,  "x_origin/F");
    outTree_->Branch("y_origin",  &y_origin_,  "y_origin/F");
    outTree_->Branch("momentum_eta_mb2", &momentum_eta_mb2_, "momentum_eta_mb2/F");
    outTree_->Branch("momentum_phi_mb2", &momentum_phi_mb2_, "momentum_phi_mb2/F");
    outTree_->Branch("momentum_eta_ecal", &momentum_eta_ecal_, "momentum_eta_ecal/F");
    outTree_->Branch("momentum_phi_ecal", &momentum_phi_ecal_, "momentum_phi_ecal/F");
    outTree_->Branch("tkIso", &tkIso_, "tkIso/F");
    outTree_->Branch("myTkIso", &myTkIso_, "myTkIso/F");
    outTree_->Branch("tkIsoNewTk", &tkIsoNewTk_, "tkIsoNewTk/F");
    outTree_->Branch("tkIsoNewDR", &tkIsoNewDR_, "tkIsoNewDR/F");


    outGENTree_ = fs_->make<TTree>("genMuonTree", "genMuonTree");
    outGENTree_->Branch("eventNumber", &eventNumber_, "eventNumber/L");
    outGENTree_->Branch("gen_eta_mb2", &gen_eta_mb2_, "gen_eta_mb2/F");
    outGENTree_->Branch("gen_phi_mb2", &gen_phi_mb2_, "gen_phi_mb2/F");
    outGENTree_->Branch("gen_eta_ecal", &gen_eta_ecal_, "gen_eta_ecal/F");
    outGENTree_->Branch("gen_phi_ecal", &gen_phi_ecal_, "gen_phi_ecal/F");
    outGENTree_->Branch("gen_eta_origin", &gen_eta_origin_, "gen_eta_origin/F");
    outGENTree_->Branch("gen_phi_origin", &gen_phi_origin_, "gen_phi_origin/F");
    outGENTree_->Branch("gen_eta", &gen_eta_, "eta/F");
    outGENTree_->Branch("gen_phi", &gen_phi_, "phi/F");
    outGENTree_->Branch("gen_pt", &gen_pt_, "gen_pt/F");
    outGENTree_->Branch("gen_d0", &gen_d0_, "gen_d0/F");
    outGENTree_->Branch("gen_lxy", &gen_lxy_, "gen_lxy/F");
    outGENTree_->Branch("gen_charge", &gen_charge_, "gen_charge/I");
    outGENTree_->Branch("x_gen_origin", &x_gen_origin_, "x_gen_origin/F");
    outGENTree_->Branch("y_gen_origin", &y_gen_origin_, "y_gen_origin/F");

    outGENTree_->Branch("gen_momentum_eta_mb2", &gen_momentum_eta_mb2_, "gen_momentum_eta_mb2/F");
    outGENTree_->Branch("gen_momentum_phi_mb2", &gen_momentum_phi_mb2_, "gen_momentum_phi_mb2/F");
    outGENTree_->Branch("gen_momentum_eta_ecal", &gen_momentum_eta_ecal_, "gen_momentum_eta_ecal/F");
    outGENTree_->Branch("gen_momentum_phi_ecal", &gen_momentum_phi_ecal_, "gen_momentum_phi_ecal/F");
    
    trkTree_ = fs_->make<TTree>("trkTree", "trkTree");
    trkTree_->Branch("eventNumber", &eventNumber_, "eventNumber/L");
    trkTree_->Branch("pt", &trk_pt_, "trk_pt/F");
    trkTree_->Branch("eta", &trk_eta_, "trk_eta/F");
    trkTree_->Branch("phi", &trk_phi_, "trk_phi/F");
    trkTree_->Branch("trk_dxy", &trk_dxy_, "trk_dxy/F");
    trkTree_->Branch("eta_ecal", &trk_eta_ecal_, "trk_eta_ecal/F");
    trkTree_->Branch("phi_ecal", &trk_phi_ecal_, "trk_phi_ecal/F");    
    trkTree_->Branch("success", &success_, "success/F");
    trkTree_->Branch("useRK", &useRK_, "useRK/F");

    theBarrel_ = initBarrel();
    thePositiveEndcap_ = initPositiveEcalEndcap();
    theNegativeEndcap_ = initNegativeEcalEndcap();
    theSurfaceAtZero_ = initSurface();

//     outGENTree_->Branch("IP3D", &ip3d_, "phi/F");
    
  }

  ~MuonTrajectoryPropagator() override {}

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;
  bool isAncestor(const reco::Candidate*, const reco::Candidate *);


private:

  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> src_;
  const edm::EDGetTokenT<std::vector<pat::Muon>> recoSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
  
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;

//   edm::ESHandle<MuonDetLayerGeometry> muonGeometry_;
  // to use PropagateToMuon class
  // https://github.com/cms-sw/cmssw/blob/8f2b35471f5cc4641db661c87eb6535195554380/MuonAnalysis/MuonAssociators/src/PropagateToMuon.cc
  const PropagateToMuonSetup st1propSetup_;
  const PropagateToMuonSetup st2propSetup_;
  const PropagateToMuonSetup genSt2propSetup_;
  
  PropagatorWithMaterial* forwardPropagatorECAL_, *forwardPropagatorECALpion_, *forwardPropagatorECALpionRK_;

  static ReferenceCountingPointer<BoundCylinder> theBarrel_;
  static ReferenceCountingPointer<BoundDisk> thePositiveEndcap_;
  static ReferenceCountingPointer<BoundDisk> theNegativeEndcap_;
  static ReferenceCountingPointer<BoundSurface> theSurfaceAtZero_;

  static const BoundCylinder& barrel() { return *theBarrel_; }
  static const BoundDisk& diskPlus() { return *thePositiveEndcap_; }
  static const BoundDisk& diskMinus() { return *theNegativeEndcap_; }
  static const BoundSurface& surfaceAtZero() { return *theSurfaceAtZero_; }

  void clearVars();
  void clearGenVars();
  void clearTrkVars();

  // output ttrees and their variables
  long unsigned int eventNumber_;

  TTree* outTree_;
  float eta_mb2_;
  float phi_mb2_;
  float eta_mb1_;
  float phi_mb1_;
  float eta_ecal_;
  float phi_ecal_;
  float momentum_eta_mb2_;
  float momentum_phi_mb2_;
  float momentum_eta_ecal_;
  float momentum_phi_ecal_;
  float eta_;
  float phi_;
  float pt_;
  float d0_;
  float ip3d_;
  int isSTA_;
  int isGlobal_;
  int isTracker_;
  float x_ecal_; 
  float y_ecal_; 
  float x_mb1_; 
  float y_mb1_; 
  float x_mb2_; 
  float y_mb2_; 
  float x_origin_;
  float y_origin_;
  float tkIso_;
  float myTkIso_;
  float tkIsoNewTk_;
  float tkIsoNewDR_;
  int charge_;

  TTree* outGENTree_;
  float gen_eta_mb2_;
  float gen_phi_mb2_;
  float gen_eta_ecal_;
  float gen_phi_ecal_;
  float gen_eta_origin_;
  float gen_phi_origin_;
  float gen_eta_;
  float gen_phi_;
  float gen_pt_;
  float gen_d0_;
  float gen_ip3d_;
  float gen_lxy_;
  float gen_momentum_eta_mb2_;
  float gen_momentum_phi_mb2_;
  float gen_momentum_eta_ecal_;
  float gen_momentum_phi_ecal_;
  float x_gen_origin_;
  float y_gen_origin_;
  int gen_charge_;
  
  TTree* trkTree_;
  float trk_pt_;
  float trk_eta_;
  float trk_phi_;
  float trk_eta_ecal_;
  float trk_phi_ecal_;
  float trk_dxy_;
  float success_;
  float useRK_;
  
};

void MuonTrajectoryPropagator::analyze(const edm::Event &evt, const edm::EventSetup & iSetup) {

  eventNumber_ = evt.id().event();
  
  clearVars();
  clearGenVars();

  // Get magnetic field
  edm::ESHandle<MagneticField> bFieldHandle;
  bFieldHandle = iSetup.getHandle(idealMagneticFieldRecordToken_);
  const MagneticField* bField = bFieldHandle.product();

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  evt.getByToken(beamSpot_, theBeamSpotHandle);
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
  math::XYZPoint bsPosition(theBeamSpot->position());
  
  const TransientTrackBuilder* theTTBuilder = &iSetup.getData(ttkToken_);

  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  evt.getByToken(src_, genParticles);

  edm::Handle<edm::View<pat::PackedGenParticle> > packed;
  evt.getByToken(packedGenToken_,packed);
  
  forwardPropagatorECAL_ = new PropagatorWithMaterial(alongMomentum, 0.1057, bField); // muon mass
  forwardPropagatorECALpion_ = new PropagatorWithMaterial(alongMomentum, 0.139 , bField, 6, false, -1, true); // pion mass
  forwardPropagatorECALpionRK_ = new PropagatorWithMaterial(alongMomentum, 0.139 , bField, 6, true, -1, true); // pion mass

  auto const st1prop = st1propSetup_.init(iSetup);
  auto const st2prop = st2propSetup_.init(iSetup);
  auto const genSt2prop = genSt2propSetup_.init(iSetup);

  for(size_t igen = 0; igen < genParticles->size();igen++){
    const reco::Candidate * genStau = &(*genParticles)[igen];
    int stau_pdgId = genStau->pdgId();
    if (abs(stau_pdgId) != 1000015 || genStau->status() != 22 ) continue;
    
    for (const reco::GenParticle &genMu : *genParticles) {
      int pdgId = genMu.pdgId();
      int status = genMu.status();
      if (abs(pdgId) != 13 || status != 1) continue;

      clearGenVars();
      const reco::Candidate * motherInPrunedCollection = genMu.mother(0) ;
      if(motherInPrunedCollection != nullptr && isAncestor( genStau , motherInPrunedCollection)) {

        math::XYZTLorentzVector p4 = genMu.p4();
        math::XYZPoint vertex = genMu.vertex();
        GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z()); 
        gen_charge_ = genMu.charge();

        float gen_d0 = -(genMu.vx()) * std::sin(float(genMu.p4().Phi())) +
                        (genMu.vy()) * std::cos(float(genMu.p4().Phi()));
    //     float gen_d0 = (genMu.vy()*genMu.px() - genMu.vx()*genMu.py())/genMu.pt();
        float gen_lxy = sqrt(pow(genMu.vx() - bsPosition.x(),2) + pow(genMu.vy() - bsPosition.y(),2));
    
        gen_eta_ = genMu.eta();
        gen_phi_ = genMu.phi();
        gen_pt_ = genMu.pt();
        gen_d0_ = gen_d0;
        gen_lxy_ = gen_lxy;
    
    //     double bField_at_vtx = bField->inTesla(vtxPos).z(); // mag field strength at the vertex location
    //     double rinv = p4.pt() / (0.3 * gen_charge_ * bField_at_vtx); // inv radius of curvature (GeV, Tesla, cm)
    
        FreeTrajectoryState initialState(vtxPos,
                                         GlobalVector(p4.px(), p4.py(), p4.pz()),
                                         gen_charge_, bField);
        TrajectoryStateOnSurface genStateAtMB2 = genSt2prop.extrapolate(initialState);

        if (genStateAtMB2.isValid()){
          gen_eta_mb2_ = genStateAtMB2.globalPosition().eta();
          gen_phi_mb2_ = genStateAtMB2.globalPosition().phi(); 
          gen_momentum_phi_mb2_ = genStateAtMB2.globalMomentum().phi();
          gen_momentum_eta_mb2_ = genStateAtMB2.globalMomentum().eta();
//           double x_gen = genStateAtMB2.globalPosition().x();      
//           double y_gen = genStateAtMB2.globalPosition().y();
//         double arc_tan = atan2(y_gen,x_gen);
//         std::cout << arc_tan << " " << gen_phi_mb2_ << std::endl;
        }
    
        TrajectoryStateOnSurface genStateAtECAL_ = forwardPropagatorECAL_->propagate(initialState, barrel());
        if (genStateAtECAL_.isValid()) {
     
          gen_eta_ecal_ = genStateAtECAL_.globalPosition().eta();
          gen_phi_ecal_ = genStateAtECAL_.globalPosition().phi();      
          gen_momentum_eta_ecal_ = genStateAtECAL_.globalMomentum().eta();
          gen_momentum_phi_ecal_ = genStateAtECAL_.globalMomentum().phi();      
      
          TrajectoryStateOnSurface genStateAtOriginSurface_ = forwardPropagatorECAL_->propagate(initialState, surfaceAtZero());
          if (genStateAtOriginSurface_.isValid()) {
            gen_eta_origin_ = genStateAtOriginSurface_.globalPosition().eta();
            gen_phi_origin_ = genStateAtOriginSurface_.globalPosition().phi();      
            x_gen_origin_ = genStateAtOriginSurface_.globalPosition().x();      
            y_gen_origin_ = genStateAtOriginSurface_.globalPosition().y();  
          }
        }      
          
        outGENTree_->Fill();  
      } // end if muon is from stau
    } // end loop on gen part = muons
  }// end loop on staus  


  edm::Handle<std::vector<pat::Muon>> recoMuons;
  evt.getByToken(recoSrc_, recoMuons);

  edm::Handle<reco::TrackCollection> isoTracks;
  evt.getByToken(tracksToken_, isoTracks);
  if (!isoTracks.isValid())
    std::cout << "no track collection" << std::endl;
    
  reco::TrackRef track;
  for (const auto& mu : *recoMuons) {
 
    clearVars();
 
    eta_ =  mu.eta();
    phi_ =  mu.phi();
    pt_ =  mu.pt();
    isSTA_ =  mu.isStandAloneMuon();
    isGlobal_ =  mu.isGlobalMuon();
    isTracker_ =  mu.isTrackerMuon();
    charge_ = mu.charge();
    // IP from https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/interface/Muon.h#L237
    d0_ = mu.dB(pat::Muon::BS2D);
    ip3d_= mu.dB(pat::Muon::BS3D);

    TrajectoryStateOnSurface stateAtMB2 = st2prop.extrapolate(mu);
    if (stateAtMB2.isValid()){
      eta_mb2_ = stateAtMB2.globalPosition().eta();
      phi_mb2_ = stateAtMB2.globalPosition().phi();
      momentum_phi_mb2_ = stateAtMB2.globalMomentum().phi();
      momentum_eta_mb2_ = stateAtMB2.globalMomentum().eta();
      x_mb2_ = stateAtMB2.globalPosition().x();      
      y_mb2_ = stateAtMB2.globalPosition().y();      
    }

    TrajectoryStateOnSurface stateAtMB1 = st1prop.extrapolate(mu);
    if (stateAtMB1.isValid()){
      eta_mb1_ = stateAtMB1.globalPosition().eta();
      phi_mb1_ = stateAtMB1.globalPosition().phi();      
      x_mb1_ = stateAtMB1.globalPosition().x();      
      y_mb1_ = stateAtMB1.globalPosition().y();      
    }

    // at ECAL surface
    bool noTrack = false;
    TString type = "";
    if (mu.isGlobalMuon()) {track = mu.globalTrack(); type = "global";}
    else if (mu.isStandAloneMuon()) {track = mu.standAloneMuon(); type = "standalone";}
    else if (mu.isTrackerMuon()) {track = mu.innerTrack(); type = "tracker";}
    else noTrack = true;
    
    float my_iso = -99.;
    float my_iso_newTk = -99.;
    float my_iso_newDR = -99;
    double theDiff_z = 0.5;
    double theDiff_r = 0.2;
    double theDR_Max = 0.3;  
    double vtx_z = mu.vz();

    if (mu.pt() > 20){
      my_iso = 0;
      // recalc trk based isolation
      for (const auto& itrack : *isoTracks) {
  //       if (deltaR(mu,itrack) > 0.32 )
  //         continue;
  // //       std::cout << "This track has: pt= " << itrack.pt() << ", eta= " << itrack.eta() 
  // //                 << ", phi= " << itrack.phi() 
  // //                 << ", deltaZ= " << fabs(vtx_z - itrack.vz())
  // //                 << ", D0= " << fabs(itrack.dxy(bsPosition)) 
  // //                 << " dr: " << deltaR(mu,itrack) << std::endl;
        if (deltaR(mu,itrack) > theDR_Max || deltaR(mu,itrack) < 0.01)
          continue;
        if (fabs(vtx_z - itrack.vz()) > theDiff_z )  
          continue;
        if (fabs(itrack.dxy(bsPosition)) > theDiff_r )  
          continue;
  // //       std::cout << "\t will use pt= " << itrack.pt() << ", eta= " << itrack.eta() << ", phi= " << itrack.phi() << std::endl;
        my_iso += itrack.pt();
      }
    }      
    myTkIso_ = my_iso;
    tkIso_ = mu.isolationR03().sumPt;

    
    // recalc trk based isolation by projecting both muon and trks to the ecal surface
    if (!noTrack ){
      my_iso_newDR = 0;
      my_iso_newTk = 0;
      
      reco::TransientTrack muTransientTrack = theTTBuilder->build(track);      
      if (muTransientTrack.isValid()) {
        FreeTrajectoryState innerMuTSOS = muTransientTrack.initialFreeState();
        TrajectoryStateOnSurface stateAtECAL_ = forwardPropagatorECAL_->propagate(innerMuTSOS, barrel());
        if (!stateAtECAL_.isValid() || (stateAtECAL_.isValid() && fabs(stateAtECAL_.globalPosition().eta()) > 1.479f)) {
           if (mu.eta() > 0.) {
             stateAtECAL_ = forwardPropagatorECAL_->propagate(innerMuTSOS, diskPlus());
           } else {
             stateAtECAL_ = forwardPropagatorECAL_->propagate(innerMuTSOS, diskMinus());
           }
        }              

        if (stateAtECAL_.isValid()) {  
          eta_ecal_ = stateAtECAL_.globalPosition().eta();
          phi_ecal_ = stateAtECAL_.globalPosition().phi();      
          x_ecal_ = stateAtECAL_.globalPosition().x();      
          y_ecal_ = stateAtECAL_.globalPosition().y();      
          momentum_phi_ecal_ = stateAtECAL_.globalMomentum().phi();
          momentum_eta_ecal_ = stateAtECAL_.globalMomentum().eta();
    
          // test propagating at a plane centered in the origin
          TrajectoryStateOnSurface stateAtOriginSurface_ = forwardPropagatorECAL_->propagate(innerMuTSOS, surfaceAtZero());
          if (stateAtOriginSurface_.isValid()) {
            x_origin_ = stateAtOriginSurface_.globalPosition().x();      
            y_origin_ = stateAtOriginSurface_.globalPosition().y();      
          }
          // recompute isolation at ECAL
          // does it make sense to propagate the tracks to the disk if the muon is in the endcap? 
          // though shouldn't happend given the dr cone
          
          if (mu.pt() > 20){
            for (const auto& itrack : *isoTracks) {
              clearTrkVars();

              if (fabs(vtx_z - itrack.vz()) > theDiff_z || fabs(itrack.dxy(bsPosition)) > theDiff_r)
                continue;
      
              // first project track to ECAL
              const reco::TransientTrack trkTransientTrack(itrack, &(*bFieldHandle));
              if (!trkTransientTrack.isValid()) continue;
              // ** up to here, same as myTrk  ** 
              FreeTrajectoryState trackTSOS = trkTransientTrack.initialFreeState();
              TrajectoryStateOnSurface trkStateAtECAL_ = forwardPropagatorECALpion_->propagate(trackTSOS, barrel());
              trk_pt_ = itrack.pt();
              trk_eta_ = itrack.eta();
              trk_phi_ = itrack.phi();
              trk_dxy_ = itrack.dxy(bsPosition);
              
              if (!trkStateAtECAL_.isValid() || (trkStateAtECAL_.isValid() && fabs(trkStateAtECAL_.globalPosition().eta()) > 1.479f)) {
                //endcap propagator
//           TrajectoryStateOnSurface innermostState = muTransientTrack.innermostMeasurementState();

                if (itrack.eta() > 0.) {
//                 if (trkTransientTrack.innermostMeasurementState().globalPosition().eta() > 0.) {
                  trkStateAtECAL_ = forwardPropagatorECALpion_->propagate(trackTSOS, diskPlus());
                } else {
                  trkStateAtECAL_ = forwardPropagatorECALpion_->propagate(trackTSOS, diskMinus());
                }
              }  
              useRK_ = 0;
              if (!trkStateAtECAL_.isValid()) {
                trkStateAtECAL_ = forwardPropagatorECALpionRK_->propagate(trackTSOS, barrel());
                useRK_ = 1;
                if (!trkStateAtECAL_.isValid()){
                  success_ = 0;
                  trkTree_->Fill();  
                  continue;
                }
              }  
              success_ = 1;
              trk_eta_ecal_ = trkStateAtECAL_.globalPosition().eta();
              trk_phi_ecal_ = trkStateAtECAL_.globalPosition().phi();      
              trkTree_->Fill();  
              
              // first compute isolation using only tracks that are propagated to ECAL 
              // still using standard coords
              float dr_tmp = deltaR(mu.eta(), mu.phi(), itrack.eta(), itrack.phi());
              if ( dr_tmp <= theDR_Max && dr_tmp > 0.01)
                my_iso_newTk += itrack.pt();
              
              // then compute isolation using propagated info 
              float dr_at_ecal = deltaR(eta_ecal_, phi_ecal_, trk_eta_ecal_, trk_phi_ecal_);
              if ( dr_at_ecal <= theDR_Max && dr_at_ecal > 0.01)
                my_iso_newDR += itrack.pt();
            }
          } // end mupT > 20  
        }
      }  
    } // end muon track
    tkIsoNewTk_ = my_iso_newTk;
    tkIsoNewDR_ = my_iso_newDR;
    outTree_->Fill();  
    
  } // end loop on reco muons
}

// -----------------------
bool MuonTrajectoryPropagator::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
  {
//     std::cout << "in isAncestor " << particle->mother(i)->pdgId() << std::endl;
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
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
  trkTree_->GetDirectory()->cd();
  trkTree_->Write();
}

void MuonTrajectoryPropagator::clearVars() {

  eta_mb2_ = -999;
  phi_mb2_= -999;
  eta_mb1_= -999;
  phi_mb1_= -999;
  eta_ecal_= -999;
  phi_ecal_= -999;
  momentum_eta_mb2_= -999;
  momentum_phi_mb2_= -999;
  momentum_eta_ecal_= -999;
  momentum_phi_ecal_= -999;
  eta_= -999;
  phi_= -999;
  pt_= -999;
  d0_= -999;
  ip3d_= -999;
  isSTA_= -999;
  isGlobal_= -999;
  isTracker_= -999;
  x_ecal_= -999; 
  y_ecal_= -999; 
  x_mb1_= -999; 
  y_mb1_= -999; 
  x_mb2_= -999; 
  y_mb2_= -999; 
  x_origin_= -999;
  y_origin_= -999;
  tkIso_= -999;
  myTkIso_= -999;
  tkIsoNewTk_= -999;
  tkIsoNewDR_= -999;  
  charge_ = -999;
}

void MuonTrajectoryPropagator::clearGenVars() {
  
  gen_eta_mb2_= -999;
  gen_phi_mb2_= -999;
  gen_eta_ecal_= -999;
  gen_phi_ecal_= -999;
  gen_eta_origin_= -999;
  gen_phi_origin_= -999;
  gen_eta_= -999;
  gen_phi_= -999;
  gen_pt_= -999;
  gen_d0_= -999;
  gen_ip3d_= -999;
  gen_lxy_= -999;
  gen_momentum_eta_mb2_= -999;
  gen_momentum_phi_mb2_= -999;
  gen_momentum_eta_ecal_= -999;
  gen_momentum_phi_ecal_= -999;
  x_gen_origin_= -999;
  y_gen_origin_= -999;
  gen_charge_= -999;
}

void MuonTrajectoryPropagator::clearTrkVars() {
  trk_pt_ = -999;
  trk_eta_ = -999;
  trk_phi_ = -999;
  trk_eta_ecal_ = -999;
  trk_phi_ecal_ = -999;
  trk_dxy_ = -999;
  useRK_ = -999;
  success_ = -999;
}


ReferenceCountingPointer<BoundCylinder> MuonTrajectoryPropagator::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> MuonTrajectoryPropagator::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> MuonTrajectoryPropagator::theNegativeEndcap_ = nullptr;
ReferenceCountingPointer<BoundSurface> MuonTrajectoryPropagator::theSurfaceAtZero_ = nullptr;


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonTrajectoryPropagator);
