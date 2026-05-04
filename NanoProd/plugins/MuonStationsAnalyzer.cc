#include "MuonStationsAnalyzer.h"

MuonStationsAnalyzer::MuonStationsAnalyzer(const edm::ParameterSet& cfg)
    : trackingGeomToken_(esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>()), 
      st1propSetup_(cfg.getParameter<edm::ParameterSet>("muPropagator1st"), consumesCollector()),
      st2propSetup_(cfg.getParameter<edm::ParameterSet>("muPropagator2nd"), consumesCollector())

{
  inputMuonCollection_ = consumes<reco::MuonCollection>(cfg.getParameter<edm::InputTag>("inputMuonCollection"));
  inputDTRecSegment4DCollection_ =
      consumes<DTRecSegment4DCollection>(cfg.getParameter<edm::InputTag>("inputDTRecSegment4DCollection"));
  inputCSCSegmentCollection_ =
      consumes<CSCSegmentCollection>(cfg.getParameter<edm::InputTag>("inputCSCSegmentCollection"));
  useTrackerMuons_ = cfg.getUntrackedParameter<bool>("useTrackerMuons");
  useGlobalMuons_ = cfg.getUntrackedParameter<bool>("useGlobalMuons");
  useTrackerMuonsNotGlobalMuons_ = cfg.getUntrackedParameter<bool>("useTrackerMuonsNotGlobalMuons");
  useGlobalMuonsNotTrackerMuons_ = cfg.getUntrackedParameter<bool>("useGlobalMuonsNotTrackerMuons");
}

MuonStationsAnalyzer::~MuonStationsAnalyzer() {}

void MuonStationsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& setup) {
  using namespace edm;
  using namespace reco;

  iEvent.getByToken(inputMuonCollection_, muonCollectionH_);
  iEvent.getByToken(inputDTRecSegment4DCollection_, dtSegmentCollectionH_);
  iEvent.getByToken(inputCSCSegmentCollection_, cscSegmentCollectionH_);

  evt_n_ = iEvent.id().event();

  geometry_ = setup.getHandle(trackingGeomToken_);
  auto const st1prop = st1propSetup_.init(setup);
  auto const st2prop = st2propSetup_.init(setup);

  for (MuonCollection::const_iterator muon = muonCollectionH_->begin(); muon != muonCollectionH_->end(); ++muon) {
//       std::cout << "muon pT = " << muon->pt() << std::endl;
    if (!muon->isMatchesValid())
      continue;
    muon_pt_  =  muon->pt();
    muon_eta_ =  muon->eta();
    muon_phi_ =  muon->phi();
    muon_is_tracker_ = muon->isTrackerMuon();
    muon_is_sta_ = muon->isStandAloneMuon();
    muon_is_glb_ = muon->isGlobalMuon();
    muon_n_chambers_ =  muon->numberOfChambersCSCorDT();
    muon_n_matches_ =  muon->numberOfMatches();  // number of chambers with matched segments
    muon_n_segments_ = 0;
    muon_n_csc_segments_ = 0;
    muon_n_dt_segments_ = 0;

    if (muon->outerTrack().isNonnull()) {
      muon_n_hits_out_ = muon->outerTrack()->numberOfValidHits();
    }  

    // this is only to count the number of segments
    for (std::vector<MuonChamberMatch>::const_iterator chamberMatch = muon->matches().begin();
         chamberMatch != muon->matches().end();
         ++chamberMatch) {
        muon_n_segments_ += chamberMatch->segmentMatches.size();
    }

    const reco::Muon &mu = (*muon);
    TrajectoryStateOnSurface stateAtMB1 = st1prop.extrapolate(mu);
    if (stateAtMB1.isValid()){
        muon_eta_at_mb1_ = stateAtMB1.globalPosition().eta();
        muon_phi_at_mb1_ = stateAtMB1.globalPosition().phi();
    }

    TrajectoryStateOnSurface stateAtMB2 = st2prop.extrapolate(mu);
    if (stateAtMB2.isValid()){
        muon_eta_at_mb2_ = stateAtMB2.globalPosition().eta();
        muon_phi_at_mb2_ = stateAtMB2.globalPosition().phi();
    }
    muonTree_->Fill();
  }

  for (DTRecSegment4DCollection::const_iterator segment = dtSegmentCollectionH_->begin();
       segment != dtSegmentCollectionH_->end();
       ++segment) {
    LocalPoint segmentLocalPosition = segment->localPosition();
    LocalVector segmentLocalDirection = segment->localDirection();
    LocalError segmentLocalPositionError = segment->localPositionError();
    LocalError segmentLocalDirectionError = segment->localDirectionError();
    bool segmentFound = false;
    

    DetId id = segment->geographicalId();
    DTChamberId chambId(id);
    seg_station_ = chambId.station();
    
    const GeomDet* det = geometry_->idToDet(id);
    if (!det) continue;
    GlobalPoint gpos = det->surface().toGlobal(segmentLocalPosition);
    seg_x_ = gpos.x();
    seg_y_ = gpos.y();
    seg_z_ = gpos.z();

    seg_eta_ = gpos.eta();
    seg_phi_ = gpos.phi();

    this_muon_pt_ = -1;
    for (MuonCollection::const_iterator muon = muonCollectionH_->begin(); muon != muonCollectionH_->end(); ++muon) {
      if (!muon->isMatchesValid())
        continue;
//       std::cout << "muon pT = " << muon->pt() << std::endl;
//       muon_pt_  =  muon->pt();
//       muon_eta_ =  muon->eta();
//       muon_phi_ =  muon->phi();
//       muon_n_chambers_ =  muon->numberOfChambersCSCorDT();
//       muon_n_matches_ =  muon->numberOfMatches();  // number of chambers with matched segments
//       muon_n_segments_ = 0;
//       muon_n_csc_segments_ = 0;
//       muon_n_dt_segments_ = 0;


      for (std::vector<MuonChamberMatch>::const_iterator chamberMatch = muon->matches().begin();
           chamberMatch != muon->matches().end();
           ++chamberMatch) {
        for (std::vector<MuonSegmentMatch>::const_iterator segmentMatch = chamberMatch->segmentMatches.begin();
             segmentMatch != chamberMatch->segmentMatches.end();
             ++segmentMatch) {
             
          if (fabs(segmentMatch->x - segmentLocalPosition.x()) < 1E-6 &&
              fabs(segmentMatch->y - segmentLocalPosition.y()) < 1E-6 &&
              fabs(segmentMatch->dXdZ - segmentLocalDirection.x() / segmentLocalDirection.z()) < 1E-6 &&
              fabs(segmentMatch->dYdZ - segmentLocalDirection.y() / segmentLocalDirection.z()) < 1E-6 &&
              fabs(segmentMatch->xErr - sqrt(segmentLocalPositionError.xx())) < 1E-6 &&
              fabs(segmentMatch->yErr - sqrt(segmentLocalPositionError.yy())) < 1E-6 &&
              fabs(segmentMatch->dXdZErr - sqrt(segmentLocalDirectionError.xx())) < 1E-6 &&
              fabs(segmentMatch->dYdZErr - sqrt(segmentLocalDirectionError.yy())) < 1E-6) {

//               std::cout << "\t matched segment at local x = " << segmentLocalPosition.x() << std::endl;
//               std::cout << "\t matched segment at global x = " << gpos.x() << std::endl;
            segmentFound = true;
            this_muon_pt_ = muon->pt();
            break;
          }
        }  // segmentMatch
        if (segmentFound)
          break;
      }  // chamberMatch
      if (segmentFound)
        break;
    }  // muon

    if (segmentFound) isUsed_ = 1;
    else isUsed_ = 0;
    isCSC_ = 0;  
    isDT_ = 1;  
    outTree_->Fill();  
  }  // dt segment
  

  for (CSCSegmentCollection::const_iterator segment = cscSegmentCollectionH_->begin();
       segment != cscSegmentCollectionH_->end();
       ++segment) {
    LocalPoint segmentLocalPosition = segment->localPosition();
    LocalVector segmentLocalDirection = segment->localDirection();
    LocalError segmentLocalPositionError = segment->localPositionError();
    LocalError segmentLocalDirectionError = segment->localDirectionError();
    bool segmentFound = false;

    DetId id = segment->geographicalId();
    const GeomDet* det = geometry_->idToDet(id);
    if (!det) continue;
    GlobalPoint gpos = det->surface().toGlobal(segmentLocalPosition);
    seg_x_ = gpos.x();
    seg_y_ = gpos.y();
    seg_z_ = gpos.z();

    seg_eta_ = gpos.eta();
    seg_phi_ = gpos.phi();

    this_muon_pt_ = -1;
    for (MuonCollection::const_iterator muon = muonCollectionH_->begin(); muon != muonCollectionH_->end(); ++muon) {
      if (!muon->isMatchesValid())
        continue;
//       muon_pt_  =  muon->pt();
//       muon_eta_ =  muon->eta();
//       muon_phi_ =  muon->phi();
//       muon_n_chambers_ =  muon->numberOfChambersCSCorDT();
//       muon_n_matches_  =  muon->numberOfMatches();  // number of chambers with matched segments
//       muon_n_segments_ = 0;
//       muon_n_csc_segments_ = 0;
//       muon_n_dt_segments_ = 0;

      // this is only to count the number of segments
//       for (std::vector<MuonChamberMatch>::const_iterator chamberMatch = muon->matches().begin();
//            chamberMatch != muon->matches().end();
//            ++chamberMatch) {
//           muon_n_segments_ += chamberMatch->segmentMatches.size();
//       }
// 
      for (std::vector<MuonChamberMatch>::const_iterator chamberMatch = muon->matches().begin();
           chamberMatch != muon->matches().end();
           ++chamberMatch) {
        for (std::vector<MuonSegmentMatch>::const_iterator segmentMatch = chamberMatch->segmentMatches.begin();
             segmentMatch != chamberMatch->segmentMatches.end();
             ++segmentMatch) {
          if (fabs(segmentMatch->x - segmentLocalPosition.x()) < 1E-6 &&
              fabs(segmentMatch->y - segmentLocalPosition.y()) < 1E-6 &&
              fabs(segmentMatch->dXdZ - segmentLocalDirection.x() / segmentLocalDirection.z()) < 1E-6 &&
              fabs(segmentMatch->dYdZ - segmentLocalDirection.y() / segmentLocalDirection.z()) < 1E-6 &&
              fabs(segmentMatch->xErr - sqrt(segmentLocalPositionError.xx())) < 1E-6 &&
              fabs(segmentMatch->yErr - sqrt(segmentLocalPositionError.yy())) < 1E-6 &&
              fabs(segmentMatch->dXdZErr - sqrt(segmentLocalDirectionError.xx())) < 1E-6 &&
              fabs(segmentMatch->dYdZErr - sqrt(segmentLocalDirectionError.yy())) < 1E-6) {
            segmentFound = true;
            this_muon_pt_ = muon->pt();
            break;
          }
        }  // segmentMatch
        if (segmentFound)
          break;
      }  // chamberMatch
      if (segmentFound)
        break;
    }  // muon

    if (segmentFound)
      isUsed_ = 1;
    else  
      isUsed_ = 0;
    isCSC_ = 1;  
    isDT_ = 0;  
    outTree_->Fill();  
  }  // csc segment
}

void MuonStationsAnalyzer::beginJob() {
  outTree_ = fs_->make<TTree>("segTree", "segTree");
  outTree_->Branch("evt_n", &evt_n_, "evt_n/I");
  outTree_->Branch("seg_x", &seg_x_, "seg_x/F");
  outTree_->Branch("seg_y", &seg_y_, "seg_x/F");
  outTree_->Branch("seg_z", &seg_z_, "seg_x/F");
  outTree_->Branch("seg_eta", &seg_eta_, "seg_eta/F");
  outTree_->Branch("seg_phi", &seg_phi_, "seg_phi/F");
  outTree_->Branch("seg_station", &seg_station_, "seg_station/I");
  outTree_->Branch("matched_muon_pt", &this_muon_pt_, "matched_muon_pt/F");
//   outTree_->Branch("muon_pt", &muon_pt_, "muon_pt/F");
//   outTree_->Branch("muon_eta", &muon_eta_, "muon_eta/F");
//   outTree_->Branch("muon_phi", &muon_phi_, "muon_phi/F");
//   outTree_->Branch("muon_n_chambers", &muon_n_chambers_, "muon_n_chambers/I");
//   outTree_->Branch("muon_n_matches", &muon_n_matches_, "muon_n_matches/I");
//   outTree_->Branch("muon_n_segments", &muon_n_segments_, "muon_n_segments/I");
//   outTree_->Branch("muon_n_csc_segments", &muon_n_csc_segments_, "muon_n_csc_segments/I");
//   outTree_->Branch("muon_n_dt_segments",  &muon_n_dt_segments_, "muon_n_dt_segments/I");
  outTree_->Branch("isUsed", &isUsed_, "isUsed/I");
  outTree_->Branch("isDT",  &isDT_,  "isDT/I");
  outTree_->Branch("isCSC", &isCSC_, "isCSC/I");


  muonTree_ = fs_->make<TTree>("muonTree", "muonTree");
  muonTree_->Branch("evt_n", &evt_n_, "evt_n/I");
  muonTree_->Branch("muon_pt", &muon_pt_, "muon_pt/F");
  muonTree_->Branch("muon_eta", &muon_eta_, "muon_eta/F");
  muonTree_->Branch("muon_phi", &muon_phi_, "muon_phi/F");
  muonTree_->Branch("muon_n_hits_out", &muon_n_hits_out_, "muon_n_hits_out/F");
  
  muonTree_->Branch("muon_n_chambers", &muon_n_chambers_, "muon_n_chambers/I");
  muonTree_->Branch("muon_n_matches", &muon_n_matches_, "muon_n_matches/I");
  muonTree_->Branch("muon_n_segments", &muon_n_segments_, "muon_n_segments/I");
  muonTree_->Branch("muon_n_csc_segments", &muon_n_csc_segments_, "muon_n_csc_segments/I");
  muonTree_->Branch("muon_n_dt_segments",  &muon_n_dt_segments_, "muon_n_dt_segments/I");
  muonTree_->Branch("muon_is_tracker",  &muon_is_tracker_, "muon_is_tracker/I");
  muonTree_->Branch("muon_is_glb",  &muon_is_glb_, "muon_is_glb/I");
  muonTree_->Branch("muon_is_sta",  &muon_is_sta_, "muon_is_sta/I");
  muonTree_->Branch("muon_eta_at_mb1",  &muon_eta_at_mb1_, "muon_eta_at_mb1/F");
  muonTree_->Branch("muon_phi_at_mb1",  &muon_phi_at_mb1_, "muon_phi_at_mb1/F");
  muonTree_->Branch("muon_eta_at_mb2",  &muon_eta_at_mb2_, "muon_eta_at_mb2/F");
  muonTree_->Branch("muon_phi_at_mb2",  &muon_phi_at_mb2_, "muon_phi_at_mb2/F");
  
}
void MuonStationsAnalyzer::endJob() {
  outTree_->GetDirectory()->cd();
  outTree_->Write();

  muonTree_->GetDirectory()->cd();
  muonTree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonStationsAnalyzer);
