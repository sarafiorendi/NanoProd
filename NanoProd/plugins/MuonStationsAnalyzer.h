// system include files
#include <string>

// user include files
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

#include "TTree.h"
#include "TROOT.h"

class MuonStationsAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MuonStationsAnalyzer(const edm::ParameterSet&);
  ~MuonStationsAnalyzer() override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  void beginJob() override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::MuonCollection> inputMuonCollection_;
  edm::EDGetTokenT<DTRecSegment4DCollection> inputDTRecSegment4DCollection_;
  edm::EDGetTokenT<CSCSegmentCollection> inputCSCSegmentCollection_;
  bool useTrackerMuons_;
  bool useGlobalMuons_;
  bool useTrackerMuonsNotGlobalMuons_;
  bool useGlobalMuonsNotTrackerMuons_;

  edm::Handle<reco::MuonCollection> muonCollectionH_;
  edm::Handle<DTRecSegment4DCollection> dtSegmentCollectionH_;
  edm::Handle<CSCSegmentCollection> cscSegmentCollectionH_;
  edm::ESHandle<GlobalTrackingGeometry> geometry_;

  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> trackingGeomToken_;

  const PropagateToMuonSetup st1propSetup_, st2propSetup_;

  edm::Service<TFileService> fs_;
  TTree* outTree_;
  TTree* muonTree_;

  unsigned long evt_n_;
  float seg_x_, seg_y_, seg_z_, seg_eta_, seg_phi_, this_muon_pt_;
  int isUsed_, isCSC_, isDT_, seg_station_;

  float muon_pt_, muon_eta_, muon_phi_;
  int muon_n_matches_, muon_n_chambers_, muon_n_segments_, muon_n_csc_segments_, muon_n_dt_segments_;
  int muon_is_tracker_, muon_is_sta_, muon_is_glb_, muon_n_hits_out_;
  float muon_eta_at_mb1_, muon_phi_at_mb1_;
  float muon_eta_at_mb2_, muon_phi_at_mb2_;
};

