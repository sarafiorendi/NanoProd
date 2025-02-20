// -*- C++ -*-
//
// Package:    todel/DileptonAnalyzer
// Class:      DileptonAnalyzer
//
/**\class DileptonAnalyzer DileptonAnalyzer.cc todel/DileptonAnalyzer/plugins/DileptonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sara Fiorendi
//         Created:  Fri, 14 Feb 2025 14:02:50 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TTree.h"
#include "TROOT.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class DileptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DileptonAnalyzer(const edm::ParameterSet&);
  ~DileptonAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  typedef ROOT::Math::SVector<double, 3> SVector3;

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::Service<TFileService> fs_;
  TTree* outTree_;

  float vtx_x_;
  float vtx_y_;
  float vtx_z_;
  float vtx_r_;
  float vtx_m_;
  float vtx_normchi2_;
  float bs_x_;
  float bs_y_;
  float vtx_r_bs_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::CompositeCandidateCollection> vtxToken_;  //used to select what vtxs to read from configuration file
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DileptonAnalyzer::DileptonAnalyzer(const edm::ParameterSet& iConfig)
    : vtxToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vtxs"))),
      beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot")))
{
  //now do what ever initialization is needed
}

DileptonAnalyzer::~DileptonAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void DileptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByToken(beamSpot_, theBeamSpotHandle);
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
  math::XYZPoint bsPosition(theBeamSpot->position());

  bs_x_ =  bsPosition.x() ;
  bs_y_ =  bsPosition.y() ;

  for (const auto& vtx : iEvent.get(vtxToken_)) {
    vtx_x_ =  vtx.userFloat("vtx_x") ;
    vtx_y_ =  vtx.userFloat("vtx_y") ;
    vtx_z_ =  vtx.userFloat("vtx_z") ;
    vtx_r_ =  vtx.userFloat("vtx_lxy") ;
    vtx_normchi2_ =  vtx.userFloat("vtx_normChi2") ;
    vtx_m_ =  vtx.mass() ;
    SVector3 distVecXY(vtx_x_ - bs_x_, vtx_y_ - bs_y_, 0.);
    vtx_r_bs_ = ROOT::Math::Mag(distVecXY);
    outTree_->Fill();  
  }
}

// ------------ method called once each job just before starting event loop  ------------
void DileptonAnalyzer::beginJob() {
  outTree_ = fs_->make<TTree>("vtxTree", "vtxTree");
  outTree_->Branch("vtx_x", &vtx_x_, "vtx_x/F");
  outTree_->Branch("vtx_y", &vtx_y_, "vtx_y/F");
  outTree_->Branch("vtx_z", &vtx_z_, "vtx_z/F");
  outTree_->Branch("vtx_r", &vtx_r_, "vtx_r/F");
  outTree_->Branch("vtx_m", &vtx_m_, "vtx_m/F");
  outTree_->Branch("vtx_normchi2", &vtx_normchi2_, "vtx_normchi2/F");
  outTree_->Branch("bs_x", &bs_x_, "bs_x/F");
  outTree_->Branch("bs_y", &bs_y_, "bs_y/F");
  outTree_->Branch("vtx_r_bs", &vtx_r_bs_, "vtx_r_bs/F");
}  

// ------------ method called once each job just after ending the event loop  ------------
void DileptonAnalyzer::endJob() {
  outTree_->GetDirectory()->cd();
  outTree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DileptonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'vtxs' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DileptonAnalyzer);