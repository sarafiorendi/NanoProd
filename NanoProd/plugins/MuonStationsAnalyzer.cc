#include "MuonStationsAnalyzer.h"

MuonStationsAnalyzer::MuonStationsAnalyzer(const edm::ParameterSet& cfg)
    : trackingGeomToken_(esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>()), 
      dtGeomToken_(esConsumes<DTGeometry, MuonGeometryRecord>()),
      st1propSetup_(cfg.getParameter<edm::ParameterSet>("muPropagator1st"), consumesCollector()),
      st2propSetup_(cfg.getParameter<edm::ParameterSet>("muPropagator2nd"), consumesCollector()),
      beamSpotToken_(consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") ))
{
  inputMuonCollection_ = consumes<reco::MuonCollection>(cfg.getParameter<edm::InputTag>("inputMuonCollection"));
  inputDTRecSegment4DCollection_ =
      consumes<DTRecSegment4DCollection>(cfg.getParameter<edm::InputTag>("inputDTRecSegment4DCollection"));
  inputCSCSegmentCollection_ =
      consumes<CSCSegmentCollection>(cfg.getParameter<edm::InputTag>("inputCSCSegmentCollection"));
  
  hits_detids_ = new std::vector<int>;
  hits_detids_ -> clear();

  invalid_hits_detids_ = new std::vector<int>;
  invalid_hits_detids_ -> clear();
}

MuonStationsAnalyzer::~MuonStationsAnalyzer() {
  delete hits_detids_;
  delete invalid_hits_detids_;
}

void MuonStationsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& setup) {
  using namespace edm;
  using namespace reco;

  iEvent.getByToken(inputMuonCollection_, muonCollectionH_);
  iEvent.getByToken(inputDTRecSegment4DCollection_, dtSegmentCollectionH_);
  iEvent.getByToken(inputCSCSegmentCollection_, cscSegmentCollectionH_);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamSpotToken_, beamSpot);

  evt_n_ = iEvent.id().event();

  geometry_ = setup.getHandle(trackingGeomToken_);
  dtGeom_ = setup.getHandle(dtGeomToken_);
  
  auto const st1prop = st1propSetup_.init(setup);
  auto const st2prop = st2propSetup_.init(setup);

  for (const auto& muon : *muonCollectionH_) {
//       std::cout << "muon pT = " << muon->pt() << std::endl;
    if (!muon.isMatchesValid())
      continue;
    clearMuonVars();  

    muon_pt_  =  muon.pt();
    muon_eta_ =  muon.eta();
    muon_phi_ =  muon.phi();
    if (beamSpot.isValid() && muon.muonBestTrack().isNonnull()) {
      muon_dz_ = muon.muonBestTrack()->dz(beamSpot->position());
    }

    muon_is_tracker_ = muon.isTrackerMuon();
    muon_is_sta_ = muon.isStandAloneMuon();
    muon_is_glb_ = muon.isGlobalMuon();
    muon_n_chambers_ =  muon.numberOfChambersCSCorDT();
    muon_n_matches_ =  muon.numberOfMatches();  // number of chambers with matched segments
    muon_n_matched_stations_ =  muon.numberOfMatchedStations();  // number of chambers with matched segments

//     std::cout << muon_pt_ << "  \t n_matches = " 
//               << muon_n_matches_ << "  \t n_matched stations = " 
//               << muon_n_matched_stations_ << "  \t" 
//               << std::endl;

    if (muon.outerTrack().isNonnull()) {
      muon_n_hits_out_ = muon.outerTrack()->numberOfValidHits();
      
      int ihit = 0;
      const reco::Track& sta = *muon.outerTrack();
      for (auto hitIt = sta.recHitsBegin(); hitIt != sta.recHitsEnd(); ++hitIt, ++ihit) {
        const TrackingRecHit& hit = **hitIt;

        const DetId detId = hit.geographicalId();
        if (detId.det() != DetId::Muon || detId.subdetId() != MuonSubdetId::DT) continue;
//         DTLayerId layerId(detId.rawId());
//         DTChamberId chamberId = layerId.chamberId();
//         std::cout
//             << "    outerTrack hit " << ihit
//             << "  valid=" << hit.isValid()
//             << "  det=" << detId.det()
//             << "  subdet=" << detId.subdetId()
//             << "  rawId=" << detId.rawId() 
//             << "  chamberId=" << chamberId.rawId() 
//             << "\n" ;

        if (!hit.isValid()) {
          invalid_hits_detids_->push_back(detId.rawId());
          continue;
        }  
        hits_detids_->push_back(detId.rawId());

        const GeomDet* geomDet = dtGeom_->idToDet(detId);
        if (!geomDet) {
          edm::LogPrint("MuonOuterTrackDTHitDump") << "      No GeomDet found for DT hit rawId=" << detId.rawId();
          continue;
        }
//         const auto localPos = hit.localPosition();
//         const auto globalPos = geomDet->surface().toGlobal(localPos);
//         std::cout 
//           << std::fixed << std::setprecision(4)
//           << "      DT hit local(x,z)= (" << localPos.x() << ", " << localPos.z() << ")"
//           << " global(x,y,z)= (" << globalPos.x() << ", " << globalPos.y() << ", " << globalPos.z() << ")\n";
      } // end loop on recHits
//       std::cout << "\t  \t number of hits for this muon = " << muon_n_hits_out_ << std::endl;
    }
    
    // this is only to count the number of segments
    for (std::vector<MuonChamberMatch>::const_iterator chamberMatch = muon.matches().begin();
         chamberMatch != muon.matches().end();
         ++chamberMatch) {
        muon_n_segments_ += chamberMatch->segmentMatches.size();
    }

    TrajectoryStateOnSurface stateAtMB1 = st1prop.extrapolate(muon);
    if (stateAtMB1.isValid()){
        muon_eta_at_mb1_ = stateAtMB1.globalPosition().eta();
        muon_phi_at_mb1_ = stateAtMB1.globalPosition().phi();
    }

    TrajectoryStateOnSurface stateAtMB2 = st2prop.extrapolate(muon);
    if (stateAtMB2.isValid()){
        muon_eta_at_mb2_ = stateAtMB2.globalPosition().eta();
        muon_phi_at_mb2_ = stateAtMB2.globalPosition().phi();
        x_at_mb2_ = stateAtMB2.globalPosition().x();
        y_at_mb2_ = stateAtMB2.globalPosition().y();
        z_at_mb2_ = stateAtMB2.globalPosition().z();   
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
    
    int valid_seg = segment->isValid();

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
    seg_detid_ = id.rawId();
    
    seg_eta_dir_ = det->surface().toGlobal(segmentLocalDirection).eta();
    seg_phi_dir_ = det->surface().toGlobal(segmentLocalDirection).phi();

    std::cout <<  " DT segment on detId " << id.rawId() << "   station = " << seg_station_ << std::endl;   
    std::cout <<  " \t eta/phi " << seg_eta_<< " / " << seg_phi_ << std::endl;   
    std::cout <<  " \t valid? " << valid_seg  << std::endl;   
    auto segment_hits = segment->recHits();
    for (auto ihit : segment_hits){
      std::cout <<  " \t\t hit: " << ihit->geographicalId().rawId() << " - valid? " << ihit->isValid()  << std::endl;   

    }
    if (segment->hasPhi())
      std::cout <<  " \t\t phi segment detID: " << segment->phiSegment()->geographicalId().rawId()  << std::endl;   
    if (segment->hasZed())
    std::cout <<  " \t\t r-z segment detID: " << segment->zSegment()->geographicalId().rawId()  << std::endl;   
    

    this_muon_pt_ = -1;
    for (MuonCollection::const_iterator muon = muonCollectionH_->begin(); muon != muonCollectionH_->end(); ++muon) {
      if (!muon->isMatchesValid())
        continue;

      for (std::vector<MuonChamberMatch>::const_iterator chamberMatch = muon->matches().begin();
           chamberMatch != muon->matches().end();
           ++chamberMatch) {
        for (std::vector<MuonSegmentMatch>::const_iterator segmentMatch = chamberMatch->segmentMatches.begin();
             segmentMatch != chamberMatch->segmentMatches.end();
             ++segmentMatch) {
             
//           std::cout <<  " \t checking segments of muon with pT " << muon->pt() << std::endl;   
//           if (segmentMatch->isMask(MuonSegmentMatch::BestInChamberByDR))  std::cout <<  " \t\tBestInChamberByDR" << std::endl;   
//           if (segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByDR)) std::cout <<  " \t\tBelongsToTrackByDR" << std::endl;   
//           if (segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByCleaning)) std::cout <<  " \t\tBelongsToTrackByCleaning\n" << std::endl;   

          // previously was -6
          if (fabs(segmentMatch->x - segmentLocalPosition.x()) < 1E-4 &&
              fabs(segmentMatch->y - segmentLocalPosition.y()) < 1E-4 &&
              fabs(segmentMatch->dXdZ - segmentLocalDirection.x() / segmentLocalDirection.z()) < 1E-4 &&
              fabs(segmentMatch->dYdZ - segmentLocalDirection.y() / segmentLocalDirection.z()) < 1E-4 &&
              fabs(segmentMatch->xErr - sqrt(segmentLocalPositionError.xx())) < 1E-4 &&
              fabs(segmentMatch->yErr - sqrt(segmentLocalPositionError.yy())) < 1E-4 &&
              fabs(segmentMatch->dXdZErr - sqrt(segmentLocalDirectionError.xx())) < 1E-4 &&
              fabs(segmentMatch->dYdZErr - sqrt(segmentLocalDirectionError.yy())) < 1E-4) {

            segmentFound = true;
            this_muon_pt_ = muon->pt();
//             if (segmentFound) std::cout <<  "\t\t this muon is matched to our segment (pt = " << this_muon_pt_ <<  ")" << std::endl;   
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
//     CSCChamber chambId(id);
    seg_station_ = -1;
//     seg_station_ = chambId.station();

   if (!det) continue;
    GlobalPoint gpos = det->surface().toGlobal(segmentLocalPosition);
    seg_x_ = gpos.x();
    seg_y_ = gpos.y();
    seg_z_ = gpos.z();

    seg_eta_ = gpos.eta();
    seg_phi_ = gpos.phi();

//     std::cout <<  " CSC segment" << std::endl;   
////     std::cout <<  " CSC segment on detId " << chambId << std::endl;   
//     std::cout <<  " \t eta/phi " << seg_eta_<< " / " << seg_phi_ << std::endl;   

    this_muon_pt_ = -1;
    for (MuonCollection::const_iterator muon = muonCollectionH_->begin(); muon != muonCollectionH_->end(); ++muon) {
      if (!muon->isMatchesValid())
        continue;
      for (std::vector<MuonChamberMatch>::const_iterator chamberMatch = muon->matches().begin();
           chamberMatch != muon->matches().end();
           ++chamberMatch) {
        for (std::vector<MuonSegmentMatch>::const_iterator segmentMatch = chamberMatch->segmentMatches.begin();
             segmentMatch != chamberMatch->segmentMatches.end();
             ++segmentMatch) {
//           std::cout <<  " checking segments of muon " << muon->pt() << std::endl;   
//           if (segmentMatch->isMask(MuonSegmentMatch::BestInChamberByDR))  std::cout <<  " \tBestInChamberByDR" << std::endl;   
//           if (segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByDR)) std::cout <<  " \tBelongsToTrackByDR" << std::endl;   
//           if (segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByCleaning)) std::cout <<  " \tBelongsToTrackByCleaning\n" << std::endl;   
          if (fabs(segmentMatch->x - segmentLocalPosition.x()) < 1E-4 &&
              fabs(segmentMatch->y - segmentLocalPosition.y()) < 1E-4 &&
              fabs(segmentMatch->dXdZ - segmentLocalDirection.x() / segmentLocalDirection.z()) < 1E-4 &&
              fabs(segmentMatch->dYdZ - segmentLocalDirection.y() / segmentLocalDirection.z()) < 1E-4 &&
              fabs(segmentMatch->xErr - sqrt(segmentLocalPositionError.xx())) < 1E-4 &&
              fabs(segmentMatch->yErr - sqrt(segmentLocalPositionError.yy())) < 1E-4 &&
              fabs(segmentMatch->dXdZErr - sqrt(segmentLocalDirectionError.xx())) < 1E-4 &&
              fabs(segmentMatch->dYdZErr - sqrt(segmentLocalDirectionError.yy())) < 1E-4) {
            segmentFound = true;
            this_muon_pt_ = muon->pt();
//             if (segmentFound) std::cout <<  "\t is matched to our segment " << this_muon_pt_ << std::endl;   
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

void MuonStationsAnalyzer::clearMuonVars(){

  muon_pt_ = -999;
  muon_eta_ = -999;
  muon_phi_ = -999;
  muon_dz_ = -999;
  muon_n_hits_out_ = -999;
  muon_n_chambers_ = -999;
  muon_n_matches_ = -999;
  muon_n_csc_segments_ = -999;
  muon_n_dt_segments_ = -999;
  muon_is_tracker_ = -999;
  muon_is_glb_ = -999;
  muon_is_sta_ = -999;
  muon_eta_at_mb1_ = -999;
  muon_phi_at_mb1_ = -999;
  muon_eta_at_mb2_ = -999;
  muon_phi_at_mb2_ = -999;

  // this must be 0 because we count it in a loop
  muon_n_segments_ = 0;

  x_at_mb2_ = -999;
  y_at_mb2_ = -999;
  z_at_mb2_ = -999;
  
  hits_detids_->clear();
  invalid_hits_detids_->clear();
}

void MuonStationsAnalyzer::beginJob() {
  outTree_ = fs_->make<TTree>("segTree", "segTree");
  outTree_->Branch("evt_n", &evt_n_, "evt_n/I");
  outTree_->Branch("seg_x", &seg_x_, "seg_x/F");
  outTree_->Branch("seg_y", &seg_y_, "seg_y/F");
  outTree_->Branch("seg_z", &seg_z_, "seg_z/F");
  outTree_->Branch("seg_eta", &seg_eta_, "seg_eta/F");
  outTree_->Branch("seg_phi", &seg_phi_, "seg_phi/F");
  outTree_->Branch("seg_eta_dir", &seg_eta_dir_, "seg_eta_dir/F");
  outTree_->Branch("seg_phi_dir", &seg_phi_dir_, "seg_phi_dir/F");
  outTree_->Branch("seg_station", &seg_station_, "seg_station/I");
  outTree_->Branch("seg_detid", &seg_detid_, "seg_detid/I");
  outTree_->Branch("matched_muon_pt", &this_muon_pt_, "matched_muon_pt/F");
  outTree_->Branch("isUsed", &isUsed_, "isUsed/I");
  outTree_->Branch("isDT",  &isDT_,  "isDT/I");
  outTree_->Branch("isCSC", &isCSC_, "isCSC/I");


  muonTree_ = fs_->make<TTree>("muonTree", "muonTree");
  muonTree_->Branch("evt_n", &evt_n_, "evt_n/I");
  muonTree_->Branch("muon_pt", &muon_pt_, "muon_pt/F");
  muonTree_->Branch("muon_eta", &muon_eta_, "muon_eta/F");
  muonTree_->Branch("muon_phi", &muon_phi_, "muon_phi/F");
  muonTree_->Branch("muon_dz", &muon_dz_, "muon_dz/F");
  muonTree_->Branch("muon_n_hits_out", &muon_n_hits_out_, "muon_n_hits_out/I");
  
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

  muonTree_->Branch("x_at_mb2",  &x_at_mb2_, "x_at_mb2/F");
  muonTree_->Branch("y_at_mb2",  &y_at_mb2_, "y_at_mb2/F");
  muonTree_->Branch("z_at_mb2",  &z_at_mb2_, "z_at_mb2/F");
  
  muonTree_->Branch("hits_detids",  &hits_detids_);
  muonTree_->Branch("invalid_hits_detids",  &invalid_hits_detids_);
  
}
void MuonStationsAnalyzer::endJob() {
//   outTree_->GetDirectory()->cd();
//   outTree_->Write();

//   muonTree_->GetDirectory()->cd();
//   muonTree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonStationsAnalyzer);
