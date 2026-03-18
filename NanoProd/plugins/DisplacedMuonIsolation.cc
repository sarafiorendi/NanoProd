/*
 * Add recomputed isolation to displaced muons
 */

#include <memory>
#include <cmath>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
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

template<typename T>
void putValueMap(edm::Event& event,
                 const edm::Handle<std::vector<pat::Muon>>& src,
                 const std::vector<T>& values,
                 const std::string& label)
{
    auto vm = std::make_unique<edm::ValueMap<T>>();
    typename edm::ValueMap<T>::Filler filler(*vm);
    filler.insert(src, values.begin(), values.end());
    filler.fill();
    event.put(std::move(vm), label);
}

class DisplacedMuonIsolation : public edm::stream::EDProducer<> {
public:
    explicit DisplacedMuonIsolation(const edm::ParameterSet&);
    ~DisplacedMuonIsolation(){};

private:
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;

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
    const PropagateToMuonSetup st2propSetup_;
    
    float theDiff_z_;
    float theDiff_r_;
    float theDR_Max_;  
    float theDR_Min_;  
    float thePt_Min_;

    struct MuonPropagationState {
      float eta_at_ecal = -99.f;
      float phi_at_ecal = -99.f;
      float eta_at_mb2 = -99.f;
      float phi_at_mb2 = -99.f;
      float x_at_mb2 = -9999.f;
      float y_at_mb2 = -9999.f;
      float z_at_mb2 = -9999.f;
      float px_at_mb2 = -9999.f;
      float py_at_mb2 = -9999.f;
      float pz_at_mb2 = -9999.f;
    };

};


DisplacedMuonIsolation::DisplacedMuonIsolation(const edm::ParameterSet& cfg)
    : ttkToken_{esConsumes(edm::ESInputTag{"", "TransientTrackBuilder"})},
      muonsToken_{consumes<std::vector<pat::Muon>>( cfg.getParameter<edm::InputTag>("muons") )},
      tracksToken_{consumes<reco::TrackCollection>( cfg.getParameter<edm::InputTag>("tracksForIso") )},    
      beamSpotToken_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
      idealMagneticFieldRecordToken_(esConsumes()),
      st2propSetup_(cfg.getParameter<edm::ParameterSet>("muPropagator2nd"), consumesCollector()),
      theDiff_z_( cfg.getUntrackedParameter<double>("deltaZ")),
      theDiff_r_( cfg.getUntrackedParameter<double>("deltaDxy")),
      theDR_Max_( cfg.getUntrackedParameter<double>("maxDeltaR")),
      theDR_Min_( cfg.getUntrackedParameter<double>("minDeltaR")),
      thePt_Min_( cfg.getUntrackedParameter<double>("minTrkPt"))
{  
    if (!theBarrel_){
      theBarrel_ = initBarrel();
      thePositiveEndcap_ = initPositiveEcalEndcap();
      theNegativeEndcap_ = initNegativeEcalEndcap();
    }  
 
    produces<edm::ValueMap<float>>("isoNewTk");
    produces<edm::ValueMap<float>>("isoNewDR"); 
    produces<edm::ValueMap<float>>("isoNewDRDz0p2");
    produces<edm::ValueMap<float>>("isoNewDRDz0p2Dxy0p1"); 
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

    auto const st2prop = st2propSetup_.init(setup);

    PropagatorWithMaterial forwardPropagatorECAL(alongMomentum, 0.1057, bField); // muon mass
    PropagatorWithMaterial forwardPropagatorECALpion(alongMomentum, 0.139 , bField, 6, false, -1, true); // pion mass
    PropagatorWithMaterial forwardPropagatorECALpionRK(alongMomentum, 0.139 , bField, 6, true, -1, true); // pion mass

    const TransientTrackBuilder* theTTBuilder = &setup.getData(ttkToken_);
    
    const size_t muons_size = recoMuons->size();
    std::vector <Float_t> v_iso0(muons_size, -99);
    std::vector <Float_t> v_iso1(muons_size, -99);
    std::vector <Float_t> v_iso2(muons_size, -99);
    std::vector <Float_t> v_iso3(muons_size, -99);
    std::vector <Float_t> v_eta_ecal(muons_size, -99);
    std::vector <Float_t> v_phi_ecal(muons_size, -99);
    std::vector <Float_t> v_eta_mb2(muons_size, -99);
    std::vector <Float_t> v_phi_mb2(muons_size, -99);
    std::vector <Float_t> v_x_mb2(muons_size, -9999);
    std::vector <Float_t> v_y_mb2(muons_size, -9999);
    std::vector <Float_t> v_z_mb2(muons_size, -9999);
    std::vector <Float_t> v_px_mb2(muons_size, -9999);
    std::vector <Float_t> v_py_mb2(muons_size, -9999);
    std::vector <Float_t> v_pz_mb2(muons_size, -9999);
    
    float my_iso_newTk = -99.;
    float my_iso_newDR = -99;
    float my_iso_newDR_dz0p2 = -99;
    float my_iso_newDR_dz0p2_dxy0p1 = -99;
    
    // loop on the muons   
    reco::TrackRef muonTrack;
    for(size_t muIndex = 0; muIndex < muons_size; ++muIndex)
    {
      const pat::Muon& mu = (*recoMuons)[muIndex];      
      MuonPropagationState muState;

      // reset iso variable
      my_iso_newTk = 0;
      my_iso_newDR = 0;
      my_iso_newDR_dz0p2 = 0;
      my_iso_newDR_dz0p2_dxy0p1 = 0;

      // project muon traj at MB2 surface
      TrajectoryStateOnSurface stateAtMB2 = st2prop.extrapolate(mu);
      if (stateAtMB2.isValid()){
        muState.eta_at_mb2 = stateAtMB2.globalPosition().eta();
        muState.phi_at_mb2 = stateAtMB2.globalPosition().phi();
        muState.x_at_mb2  = stateAtMB2.globalPosition().x();
        muState.y_at_mb2  = stateAtMB2.globalPosition().y();
        muState.z_at_mb2  = stateAtMB2.globalPosition().z();
        muState.px_at_mb2 = stateAtMB2.globalMomentum().x();
        muState.py_at_mb2 = stateAtMB2.globalMomentum().y();
        muState.pz_at_mb2 = stateAtMB2.globalMomentum().z();
      }
      
      // project muon traj at ECAL surface
      bool noTrack = false;
      if (mu.isGlobalMuon()) muonTrack = mu.globalTrack(); 
      else if (mu.isStandAloneMuon()) muonTrack = mu.standAloneMuon(); 
      else if (mu.isTrackerMuon()) muonTrack = mu.innerTrack(); 
      else noTrack = true;
      
      if (!noTrack ){

        reco::TransientTrack muTransientTrack = theTTBuilder->build(muonTrack);      
        if (muTransientTrack.isValid()) {
          FreeTrajectoryState innerMuTSOS = muTransientTrack.initialFreeState();
          TrajectoryStateOnSurface stateAtECAL_ = forwardPropagatorECAL.propagate(innerMuTSOS, barrel());
          if (!stateAtECAL_.isValid() || (std::abs(stateAtECAL_.globalPosition().eta()) > 1.479f)) {
             if (mu.eta() > 0.) {
               stateAtECAL_ = forwardPropagatorECAL.propagate(innerMuTSOS, diskPlus());
             } else {
               stateAtECAL_ = forwardPropagatorECAL.propagate(innerMuTSOS, diskMinus());
             }
          }
          if (stateAtECAL_.isValid()) {  
            muState.eta_at_ecal = stateAtECAL_.globalPosition().eta();
            muState.phi_at_ecal = stateAtECAL_.globalPosition().phi();  
          
            // loop on tracks to build the isolation
            for (const auto& itrack : *isoTracks) {
            
              if (itrack.pt() < thePt_Min_) 
                continue;

              float trk_dxy = std::abs(itrack.dxy(bsPosition));

              if (std::abs(mu.vz() - itrack.vz()) > theDiff_z_ || trk_dxy > theDiff_r_)
                continue;
      
              // project track to ECAL
              const reco::TransientTrack trkTransientTrack(itrack, bField);
              if (!trkTransientTrack.isValid()) continue;
              FreeTrajectoryState trackTSOS = trkTransientTrack.initialFreeState();
              TrajectoryStateOnSurface trkStateAtECAL_ = forwardPropagatorECALpion.propagate(trackTSOS, barrel());
              if (!trkStateAtECAL_.isValid() || (std::abs(trkStateAtECAL_.globalPosition().eta()) > 1.479f)) {
                if (itrack.eta() > 0.) {
                  trkStateAtECAL_ = forwardPropagatorECALpion.propagate(trackTSOS, diskPlus());
                } else {
                  trkStateAtECAL_ = forwardPropagatorECALpion.propagate(trackTSOS, diskMinus());
                }
              }  
              if (!trkStateAtECAL_.isValid()) {
                trkStateAtECAL_ = forwardPropagatorECALpionRK.propagate(trackTSOS, barrel());
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
              float dr_at_ecal = deltaR(muState.eta_at_ecal, muState.phi_at_ecal, trk_eta_ecal_, trk_phi_ecal_);
              if ( dr_at_ecal <= theDR_Max_ && dr_at_ecal > theDR_Min_){
                my_iso_newDR += itrack.pt();
                
                if ( std::abs(mu.vz() - itrack.vz()) < 0.2){
                  my_iso_newDR_dz0p2 += itrack.pt();
                  if ( trk_dxy < 0.1)
                    my_iso_newDR_dz0p2_dxy0p1 += itrack.pt();
                }    
              }      

            } // end loop on tracks
          } // end if muon prop is valid
        } // end if mu tt is valid 
      } // end if mu has a track      
      
      v_iso0[muIndex] = my_iso_newTk;
      v_iso1[muIndex] = my_iso_newDR;
      v_iso2[muIndex] = my_iso_newDR_dz0p2;
      v_iso3[muIndex] = my_iso_newDR_dz0p2_dxy0p1;
      v_eta_ecal[muIndex] = muState.eta_at_ecal;
      v_phi_ecal[muIndex] = muState.phi_at_ecal;
      v_eta_mb2[muIndex]  = muState.eta_at_mb2;
      v_phi_mb2[muIndex]  = muState.phi_at_mb2;
      v_x_mb2[muIndex]    = muState.x_at_mb2;
      v_y_mb2[muIndex]    = muState.y_at_mb2;
      v_z_mb2[muIndex]    = muState.z_at_mb2;
      v_px_mb2[muIndex]   = muState.px_at_mb2;
      v_py_mb2[muIndex]   = muState.py_at_mb2;
      v_pz_mb2[muIndex]   = muState.pz_at_mb2;
    } // end loop on muons 

    putValueMap(event, recoMuons, v_iso0, "isoNewTk");
    putValueMap(event, recoMuons, v_iso1, "isoNewDR");
    putValueMap(event, recoMuons, v_iso2, "isoNewDRDz0p2");
    putValueMap(event, recoMuons, v_iso3, "isoNewDRDz0p2Dxy0p1");
    
    putValueMap(event, recoMuons, v_eta_ecal, "etaAtEcal");
    putValueMap(event, recoMuons, v_phi_ecal, "phiAtEcal");
    
    putValueMap(event, recoMuons, v_eta_mb2, "etaAtMB2");
    putValueMap(event, recoMuons, v_phi_mb2, "phiAtMB2");
    putValueMap(event, recoMuons, v_x_mb2, "xAtMB2");
    putValueMap(event, recoMuons, v_y_mb2, "yAtMB2");
    putValueMap(event, recoMuons, v_z_mb2, "zAtMB2");
    putValueMap(event, recoMuons, v_px_mb2, "pxAtMB2");
    putValueMap(event, recoMuons, v_py_mb2, "pyAtMB2");
    putValueMap(event, recoMuons, v_pz_mb2, "pzAtMB2");
}

ReferenceCountingPointer<BoundCylinder> DisplacedMuonIsolation::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> DisplacedMuonIsolation::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> DisplacedMuonIsolation::theNegativeEndcap_ = nullptr;

DEFINE_FWK_MODULE(DisplacedMuonIsolation);
