/*
 * Add recomputed isolation to displaced muons
 */

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
//#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

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


class DisplacedMuonIsolation : public edm::stream::EDProducer<> {
public:
    explicit DisplacedMuonIsolation(const edm::ParameterSet&);
    ~DisplacedMuonIsolation(){};

private:
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
    PropagatorWithMaterial* forwardPropagatorECAL_, *forwardPropagatorECALpion_, *forwardPropagatorECALpionRK_;

    static ReferenceCountingPointer<BoundCylinder> theBarrel_;
    static ReferenceCountingPointer<BoundDisk> thePositiveEndcap_;
    static ReferenceCountingPointer<BoundDisk> theNegativeEndcap_;

    static const BoundCylinder& barrel() { return *theBarrel_; }
    static const BoundDisk& diskPlus() { return *thePositiveEndcap_; }
    static const BoundDisk& diskMinus() { return *theNegativeEndcap_; }

    void produce(edm::Event&, const edm::EventSetup&) override;

    const edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
    const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;
    
    double theDiff_z_;
    double theDiff_r_;
    double theDR_Max_;  
    double theDR_Min_;  
    double thePt_Min_;


};


DisplacedMuonIsolation::DisplacedMuonIsolation(const edm::ParameterSet& cfg)
    : ttkToken_{esConsumes(edm::ESInputTag{"", "TransientTrackBuilder"})},
      muonsToken_{consumes<std::vector<pat::Muon>>( cfg.getParameter<edm::InputTag>("muons") )},
      tracksToken_{consumes<reco::TrackCollection>( cfg.getParameter<edm::InputTag>("tracksForIso") )},    
      beamSpotToken_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
      idealMagneticFieldRecordToken_(esConsumes()),
      theDiff_z_( cfg.getUntrackedParameter<double>("deltaZ")),
      theDiff_r_( cfg.getUntrackedParameter<double>("deltaDxy")),
      theDR_Max_( cfg.getUntrackedParameter<double>("maxDeltaR")),
      theDR_Min_( cfg.getUntrackedParameter<double>("minDeltaR")),
      thePt_Min_( cfg.getUntrackedParameter<double>("minTrkPt"))
{  
    theBarrel_ = initBarrel();
    thePositiveEndcap_ = initPositiveEcalEndcap();
    theNegativeEndcap_ = initNegativeEcalEndcap();
 
    produces<edm::ValueMap<float>>("isoNewTk");
    produces<edm::ValueMap<float>>("isoNewDR"); 
    produces<edm::ValueMap<float>>("isoNewDRDz0p2");
    produces<edm::ValueMap<float>>("isoNewDRDz0p2Dxy0p1"); 
}

void DisplacedMuonIsolation::produce(edm::Event& event, const edm::EventSetup& setup) {

    edm::Handle<std::vector<pat::Muon>> recoMuons;
    event.getByToken(muonsToken_, recoMuons);

    edm::Handle<reco::TrackCollection> isoTracks;
    event.getByToken(tracksToken_, isoTracks);
    if (!isoTracks.isValid()) return;

    edm::Handle<reco::BeamSpot> theBeamSpotHandle;
    event.getByToken(beamSpotToken_, theBeamSpotHandle);
    const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
    math::XYZPoint bsPosition(theBeamSpot->position());

    // Get magnetic field
    edm::ESHandle<MagneticField> bFieldHandle;
    bFieldHandle = setup.getHandle(idealMagneticFieldRecordToken_);
    const MagneticField* bField = bFieldHandle.product();

    forwardPropagatorECAL_ = new PropagatorWithMaterial(alongMomentum, 0.1057, bField); // muon mass
    forwardPropagatorECALpion_ = new PropagatorWithMaterial(alongMomentum, 0.139 , bField, 6, false, -1, true); // pion mass
    forwardPropagatorECALpionRK_ = new PropagatorWithMaterial(alongMomentum, 0.139 , bField, 6, true, -1, true); // pion mass

    const TransientTrackBuilder* theTTBuilder = &setup.getData(ttkToken_);
    
    const size_t muons_size = recoMuons->size();
    std::vector <Float_t> v_iso0(muons_size, -9);
    std::vector <Float_t> v_iso1(muons_size, -9);
    std::vector <Float_t> v_iso2(muons_size, -9);
    std::vector <Float_t> v_iso3(muons_size, -9);
    
    float my_iso_newTk = -99.;
    float my_iso_newDR = -99;
    float my_iso_newDR_dz0p2 = -99;
    float my_iso_newDR_dz0p2_dxy0p1 = -99;
    
    // loop on the muons   
    reco::TrackRef muonTrack;
    for(size_t muIndex = 0; muIndex < muons_size; ++muIndex)
    {
      const auto& mu = recoMuons->at(muIndex);
      
      // project muon traj at ECAL surface
      bool noTrack = false;
      TString type = "";
      if (mu.isGlobalMuon()) {muonTrack = mu.globalTrack(); type = "global";}
      else if (mu.isStandAloneMuon()) {muonTrack = mu.standAloneMuon(); type = "standalone";}
      else if (mu.isTrackerMuon()) {muonTrack = mu.innerTrack(); type = "tracker";}
      else noTrack = true;
      
      if (!noTrack ){

        // reset iso variable
        my_iso_newTk = 0;
        my_iso_newDR = 0;
        my_iso_newDR_dz0p2 = 0;
        my_iso_newDR_dz0p2_dxy0p1 = 0;

        reco::TransientTrack muTransientTrack = theTTBuilder->build(muonTrack);      
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
            float eta_ecal_ = stateAtECAL_.globalPosition().eta();
            float phi_ecal_ = stateAtECAL_.globalPosition().phi();      
          
            // loop on tracks to build the isolation
            for (const auto& itrack : *isoTracks) {
            
              if (itrack.pt() < thePt_Min_) 
                continue;

              if (fabs(mu.vz() - itrack.vz()) > theDiff_z_ || fabs(itrack.dxy(bsPosition)) > theDiff_r_)
                continue;
      
              // project track to ECAL
              const reco::TransientTrack trkTransientTrack(itrack, &(*bFieldHandle));
              if (!trkTransientTrack.isValid()) continue;
              FreeTrajectoryState trackTSOS = trkTransientTrack.initialFreeState();
              TrajectoryStateOnSurface trkStateAtECAL_ = forwardPropagatorECALpion_->propagate(trackTSOS, barrel());
              if (!trkStateAtECAL_.isValid() || (trkStateAtECAL_.isValid() && fabs(trkStateAtECAL_.globalPosition().eta()) > 1.479f)) {
                if (itrack.eta() > 0.) {
                  trkStateAtECAL_ = forwardPropagatorECALpion_->propagate(trackTSOS, diskPlus());
                } else {
                  trkStateAtECAL_ = forwardPropagatorECALpion_->propagate(trackTSOS, diskMinus());
                }
              }  
              if (!trkStateAtECAL_.isValid()) {
                trkStateAtECAL_ = forwardPropagatorECALpionRK_->propagate(trackTSOS, barrel());
                if (!trkStateAtECAL_.isValid()){
                  continue;
                }
              }  
              float trk_eta_ecal_ = trkStateAtECAL_.globalPosition().eta();
              float trk_phi_ecal_ = trkStateAtECAL_.globalPosition().phi();      
              
              // first compute isolation using only tracks that are propagated to ECAL 
              // still using standard coords
              float dr_tmp = deltaR(mu.eta(), mu.phi(), itrack.eta(), itrack.phi());
              if ( dr_tmp <= theDR_Max_ && dr_tmp > theDR_Min_)
                my_iso_newTk += itrack.pt();
              
              // then compute isolation using propagated info 
              float dr_at_ecal = deltaR(eta_ecal_, phi_ecal_, trk_eta_ecal_, trk_phi_ecal_);
              if ( dr_at_ecal <= theDR_Max_ && dr_at_ecal > theDR_Min_){
                my_iso_newDR += itrack.pt();
                
                if ( fabs(mu.vz() - itrack.vz()) < 0.2){
                  my_iso_newDR_dz0p2 += itrack.pt();
                  if ( fabs(itrack.dxy(bsPosition)) < 0.1)
                    my_iso_newDR_dz0p2_dxy0p1 += itrack.pt();
                }    
              }      

            } // end loop on tracks
          } // end if muon prop is valid
        } // end if mu tt is valid 
      } // end if mu has a track      
      
      v_iso0.at(muIndex) = my_iso_newTk;
      v_iso1.at(muIndex) = my_iso_newDR;
      v_iso2.at(muIndex) = my_iso_newDR_dz0p2;
      v_iso3.at(muIndex) = my_iso_newDR_dz0p2_dxy0p1;
    } // end loop on muons 

    std::unique_ptr<edm::ValueMap<float>> vm_iso0(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_iso0(*vm_iso0);
    filler_iso0.insert(recoMuons, v_iso0.begin(), v_iso0.end());
    filler_iso0.fill();
    event.put(std::move(vm_iso0), "isoNewTk"); 
    
    std::unique_ptr<edm::ValueMap<float>> vm_iso1(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_iso1(*vm_iso1);
    filler_iso1.insert(recoMuons, v_iso1.begin(), v_iso1.end());
    filler_iso1.fill();
    event.put(std::move(vm_iso1), "isoNewDR"); 

    std::unique_ptr<edm::ValueMap<float>> vm_iso2(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_iso2(*vm_iso2);
    filler_iso2.insert(recoMuons, v_iso2.begin(), v_iso2.end());
    filler_iso2.fill();
    event.put(std::move(vm_iso2), "isoNewDRDz0p2"); 

    std::unique_ptr<edm::ValueMap<float>> vm_iso3(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_iso3(*vm_iso3);
    filler_iso3.insert(recoMuons, v_iso3.begin(), v_iso3.end());
    filler_iso3.fill();
    event.put(std::move(vm_iso3), "isoNewDRDz0p2Dxy0p1"); 
}

ReferenceCountingPointer<BoundCylinder> DisplacedMuonIsolation::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> DisplacedMuonIsolation::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> DisplacedMuonIsolation::theNegativeEndcap_ = nullptr;

DEFINE_FWK_MODULE(DisplacedMuonIsolation);
