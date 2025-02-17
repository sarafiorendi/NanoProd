// Take two leptons and fit them to a vertex
// Produces a collection of di-lepton vertices
// spin-off of BParkingNano framework and V0Fitter
#include "FWCore/Framework/interface/stream/EDProducer.h"
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
#include "DataFormats/Math/interface/deltaR.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


template<typename Lepton>
class DileptonVertex : public edm::stream::EDProducer<> {

public:
  typedef std::vector<Lepton> LeptonCollection;
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  explicit DileptonVertex(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<LeptonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    beamSpot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
    idealMagneticFieldRecordToken_(esConsumes()) {
      produces<pat::CompositeCandidateCollection>();
  }

  ~DileptonVertex() override {}

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  const StringCutObjectSelector<Lepton> l1_selection_; // cut on leading lepton
  const StringCutObjectSelector<Lepton> l2_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<LeptonCollection> src_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;

};

template<typename Lepton>
void DileptonVertex<Lepton>::produce(edm::Event &evt, const edm::EventSetup & iSetup) {

  // Get magnetic field
  edm::ESHandle<MagneticField> bFieldHandle;
  bFieldHandle = iSetup.getHandle(idealMagneticFieldRecordToken_);

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  evt.getByToken(beamSpot_, theBeamSpotHandle);
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
  math::XYZPoint bsPosition(theBeamSpot->position());

  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(src_, leptons);

  std::unique_ptr<pat::CompositeCandidateCollection> output(new pat::CompositeCandidateCollection());

  for (size_t l1_idx = 0; l1_idx < leptons->size(); ++l1_idx) {
    edm::Ptr<Lepton> l1_ptr(leptons, l1_idx);
    if (!l1_selection_(*l1_ptr)) continue;

    for (size_t l2_idx = l1_idx + 1; l2_idx < leptons->size(); ++l2_idx) {
      edm::Ptr<Lepton> l2_ptr(leptons, l2_idx);
      if (!l2_selection_(*l2_ptr)) continue;

      pat::CompositeCandidate lepton_pair;
      lepton_pair.setP4(l1_ptr->p4() + l2_ptr->p4());
      lepton_pair.setCharge(l1_ptr->charge() + l2_ptr->charge());
      lepton_pair.addUserFloat("lep_deltaR", reco::deltaR(*l1_ptr, *l2_ptr));

      // associate to the dilepton the pointer to single leptons
      lepton_pair.addUserCand("l1", l1_ptr );
      lepton_pair.addUserCand("l2", l2_ptr );
      // preselection before fitting 
      if ( !pre_vtx_selection_(lepton_pair) ) continue; 
      
      // build transient tracks needed for the fit 
      const reco::TransientTrack l1_tt(*l1_ptr->bestTrack(), &(*bFieldHandle));
      const reco::TransientTrack l2_tt(*l2_ptr->bestTrack(), &(*bFieldHandle));

      auto const& l1_impact = l1_tt.impactPointTSCP();
      auto const& l2_impact = l2_tt.impactPointTSCP();
      if (!l1_impact.isValid() || !l2_impact.isValid())
        continue;
      FreeTrajectoryState const& l1State = l1_impact.theState();
      FreeTrajectoryState const& l2State = l2_impact.theState();
      ClosestApproachInRPhi cApp;
      cApp.calculate(l1State, l2State);
      if (!cApp.status())
        continue;      

      // float dca = std::abs(cApp.distance()); // not used

      // the point of closest approach should at least be in the sensitive volume
      GlobalPoint cxPt = cApp.crossingPoint();
      const double cxPtR2 = cxPt.x() * cxPt.x() + cxPt.y() * cxPt.y();
      if (cxPtR2 > 120. * 120. || std::abs(cxPt.z()) > 300.)
        continue;

      // fill vector of TransientTracks to be fit, and fit
      std::vector<reco::TransientTrack> transTracks;
      transTracks.reserve(2);
      transTracks.push_back(l1_tt);
      transTracks.push_back(l2_tt);
      const GlobalError dummyError(1.0e-3, 0.0, 1.0e-3, 0.0, 0.0, 1.0e-3);
      TransientVertex theRecoVertex(cxPt, dummyError, transTracks, 1.0e-3);
      KalmanVertexFitter theKalmanFitter( false ); // useRefTracks_
      theRecoVertex = theKalmanFitter.vertex(transTracks);
      if (!theRecoVertex.isValid())
        continue;

      // from https://github.com/cms-sw/cmssw/blob/69407d5cf6621ac745631dd97ec7d4ef10ea019b/RecoVertex/V0Producer/python/generalV0Candidates_cfi.py#L4
      reco::Vertex theVtx = theRecoVertex;
      GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());
      lepton_pair.setVertex(reco::Candidate::Point((vtxPos)));

     // 2D decay significance
      SVector3 distVecXY(vtxPos.x() - bsPosition.x(), vtxPos.y() - bsPosition.y(), 0.);
      double distMagXY = ROOT::Math::Mag(distVecXY);
      lepton_pair.addUserFloat("vtx_normChi2", theVtx.normalizedChi2());
//       lepton_pair.addUserFloat("sv_prob", fitter.prob());
      lepton_pair.addUserFloat("vtx_x", lepton_pair.vx());
      lepton_pair.addUserFloat("vtx_y", lepton_pair.vy());
      lepton_pair.addUserFloat("vtx_z", lepton_pair.vz());
      lepton_pair.addUserFloat("vtx_lxy", distMagXY);

      // cut on the vertex info
      if ( !post_vtx_selection_(lepton_pair) ) continue;
      output->push_back(lepton_pair);
    }
  }

  evt.put(std::move(output));
}

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
typedef DileptonVertex<pat::Muon> DimuonVertex;
typedef DileptonVertex<pat::Electron> DielectronVertex;

template<>
void DileptonVertex<pat::Muon>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("slimmedMuons"));
  desc.add<std::string>("lep1Selection", "pt > 2");
  desc.add<std::string>("lep2Selection", "pt > 2");
  desc.add<std::string>("preVtxSelection", "charge() == 0");
  desc.add<std::string>("postVtxSelection", "");
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  descriptions.add("dimuonVertices", desc);
}

template<>
void DileptonVertex<pat::Electron>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("slimmedElectrons"));
  desc.add<std::string>("lep1Selection", "pt > 2");
  desc.add<std::string>("lep2Selection", "pt > 2");
  desc.add<std::string>("preVtxSelection", "charge() == 0");
  desc.add<std::string>("postVtxSelection", "");
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  descriptions.add("dieleVertices", desc);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DimuonVertex);
DEFINE_FWK_MODULE(DielectronVertex);
