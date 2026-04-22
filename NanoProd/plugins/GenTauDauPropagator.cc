/*
 * Add propagated coords to GEN level muons, only if coming from a stau -> tau decay
 */

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

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

#include "TMath.h"

#include <iostream>
#include <fstream>
#include <cstring>


constexpr float epsilon = 0.001;
/** Hard-wired numbers defining the surfaces on which the crystal front faces lie. */
// https://github.com/cms-sw/cmssw/blob/6cad17bc389ac53d9b2193d99f5a76713a7a99f4/RecoEgamma/EgammaPhotonAlgos/src/ConversionTrackEcalImpactPoint.cc#L73-L82
constexpr float barrelRadius = 129.f;       // p81, p50, ECAL TDR
// constexpr float barrelHalfLength = 320.5f;  // p81, p50, ECAL TDR
constexpr float barrelHalfLength = 270.9f;  // p81, p50, ECAL TDR
constexpr float endcapRadius = 171.1f;      // fig 3.26, p81, ECAL TDR
constexpr float endcapZ = 320.5f;           // fig 3.26, p81, ECAL TDR

// tmp from https://github.com/cms-sw/cmssw/blob/d229fece154f10d622a2d913cda332d6875bea27/SimG4CMS/Tracker/src/TkAccumulatingSensitiveDetector.cc#L48-L49
constexpr float trackerRadius = 120.f;     
constexpr float trackerZ = 300.f;          


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


class GenTauDauPropagator : public edm::stream::EDProducer<> {
public:
    explicit GenTauDauPropagator(const edm::ParameterSet&);
    ~GenTauDauPropagator(){};
    bool isAncestor(const reco::Candidate*, const reco::Candidate *);
    float recalculate_phi_at_DV(
        float refx, float refy, float refz, // initial position in cm
        float px, float py, float pz, // initial momentum in GeV
        int charge, // -1 or 1, matches Muon_charge
        float dvx, float dvy // DV coordinates to propagate to
    );

private:
    PropagatorWithMaterial* forwardPropagatorECALpion_;

    static ReferenceCountingPointer<BoundCylinder> theBarrel_;
    static ReferenceCountingPointer<BoundDisk> thePositiveEndcap_;
    static ReferenceCountingPointer<BoundDisk> theNegativeEndcap_;

    static const BoundCylinder& barrel() { return *theBarrel_; }
    static const BoundDisk& diskPlus() { return *thePositiveEndcap_; }
    static const BoundDisk& diskMinus() { return *theNegativeEndcap_; }

    void produce(edm::Event&, const edm::EventSetup&) override;

    const edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> const idealMagneticFieldRecordToken_;
};


GenTauDauPropagator::GenTauDauPropagator(const edm::ParameterSet& cfg)
    : 
      genPartToken_{consumes<std::vector<reco::GenParticle>>( cfg.getParameter<edm::InputTag>("src") )},
      idealMagneticFieldRecordToken_(esConsumes())
{  
    theBarrel_ = initBarrel();
    thePositiveEndcap_ = initPositiveEcalEndcap();
    theNegativeEndcap_ = initNegativeEcalEndcap();
 
    produces<edm::ValueMap<float>>("etaAtEcal"); 
    produces<edm::ValueMap<float>>("phiAtEcal"); 
}

void GenTauDauPropagator::produce(edm::Event& event, const edm::EventSetup& setup) {

    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    event.getByToken(genPartToken_, genParticles);
   
    // Get magnetic field
    edm::ESHandle<MagneticField> bFieldHandle;
    bFieldHandle = setup.getHandle(idealMagneticFieldRecordToken_);
    const MagneticField* bField = bFieldHandle.product();

    forwardPropagatorECALpion_ = new PropagatorWithMaterial(anyDirection, 0.139 , bField, 6, false, -1, true); // pion mass

    const size_t genPart_size = genParticles->size();
    std::vector <Float_t> v_eta_ecal(genPart_size, -99);
    std::vector <Float_t> v_phi_ecal(genPart_size, -99);

    float gen_eta_ecal = -99;
    float gen_phi_ecal = -99;
    float gen_charge_ = -99;
    
    for(size_t igen = 0; igen < genParticles->size();igen++){
      const reco::Candidate * genStau = &(*genParticles)[igen];
      int stau_pdgId = genStau->pdgId();
      if (abs(stau_pdgId) != 1000015 || genStau->status() != 22 ) continue;
      
      for (size_t jgen = 0; jgen < genParticles->size();jgen++) {
      
        // reset output variables
        gen_eta_ecal = -99;
        gen_phi_ecal = -99;

        const reco::GenParticle genPion = genParticles->at(jgen);
        int pdgId = genPion.pdgId();
        int status = genPion.status();
        if (abs(pdgId) != 211 || genPion.pt() < 1) {
          gen_phi_ecal = -59;
          gen_eta_ecal = -59;
          v_eta_ecal.at(jgen) = gen_eta_ecal;
          v_phi_ecal.at(jgen) = gen_phi_ecal;
          continue;
        }  
  
        const reco::Candidate * motherInPrunedCollection = genPion.mother(0) ;
        if(motherInPrunedCollection != nullptr && isAncestor( genStau , motherInPrunedCollection)) {
  
          math::XYZTLorentzVector p4 = genPion.p4();
          math::XYZPoint vertex = genPion.vertex();
          double vertex_radial_position = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
          if (vertex_radial_position > trackerRadius || abs(vertex.z()) > trackerZ) {
            gen_phi_ecal = -59;
            gen_eta_ecal = -59;      
            v_eta_ecal.at(jgen) = gen_eta_ecal;
            v_phi_ecal.at(jgen) = gen_phi_ecal;
            continue;
          }  

          GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z()); 
          gen_charge_ = genPion.charge();
          FreeTrajectoryState initialState(vtxPos,
                                           GlobalVector(p4.px(), p4.py(), p4.pz()),
                                           gen_charge_, bField);
      
//           std::cout << "\n[Main] start prop for pt / eta = " << genPion.pt() << " / " << genPion.eta() << std::endl;
//           std::cout << "                    vx vy vz = " << vertex.x() << " " << vertex.y() << " " << vertex.z() << 
//                        "  radius = " << vertex_radial_position  << std::endl;
          TrajectoryStateOnSurface genStateAtECAL_ = forwardPropagatorECALpion_->propagate(initialState, barrel());
          if (!genStateAtECAL_.isValid() || (genStateAtECAL_.isValid() && fabs(genStateAtECAL_.globalPosition().eta()) > 1.479f)) {
//              std::cout << "\ttrying prop to endcap "  << std::endl;
             if (genStateAtECAL_.isValid() && fabs(genStateAtECAL_.globalPosition().eta()) > 1.479f)
//                std::cout << "\t   in fact propagated eta is = " << genStateAtECAL_.globalPosition().eta() << std::endl;
             if (genPion.eta() > 0.) {
               genStateAtECAL_ = forwardPropagatorECALpion_->propagate(initialState, diskPlus());
//                std::cout << "\t end prop positive endcap: "  ;
//                std::cout << "status = " << genStateAtECAL_.isValid() << std::endl;
             } else {
               genStateAtECAL_ = forwardPropagatorECALpion_->propagate(initialState, diskMinus());
//                std::cout << "\t end prop negative endcap: "  ;
//                std::cout << "status = " << genStateAtECAL_.isValid() << std::endl;
             }
          }              
          if (genStateAtECAL_.isValid()) {
            gen_eta_ecal = genStateAtECAL_.globalPosition().eta();
            gen_phi_ecal = genStateAtECAL_.globalPosition().phi();      
          }   
          else{
//             std::cout << "\t !!! this one has failed: pt /eta " << genPion.pt() << " / " << genPion.eta() << std::endl;
          }   
            
        } // end if muon is from stau
        if (gen_eta_ecal > -99){
          v_eta_ecal.at(jgen) = gen_eta_ecal;
          v_phi_ecal.at(jgen) = gen_phi_ecal;
        }
      } // end loop on gen part = muons
    } // end loop on GenPart  

    // now save muon coordinates at ECAL surface
    std::unique_ptr<edm::ValueMap<float>> vm_eta(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_eta(*vm_eta);
    filler_eta.insert(genParticles, v_eta_ecal.begin(), v_eta_ecal.end());
    filler_eta.fill();
    event.put(std::move(vm_eta), "etaAtEcal"); 

    std::unique_ptr<edm::ValueMap<float>> vm_phi(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_phi(*vm_phi);
    filler_phi.insert(genParticles, v_phi_ecal.begin(), v_phi_ecal.end());
    filler_phi.fill();
    event.put(std::move(vm_phi), "phiAtEcal"); 
}

// -----------------------
float GenTauDauPropagator::recalculate_phi_at_DV(
  float refx, float refy, float refz, // initial position in cm
  float px, float py, float pz, // initial momentum in GeV
  int charge, // -1 or 1, matches Muon_charge
  float dvx, float dvy // DV coordinates to propagate to
) {

  // Performs a helix propagation from 
  // track reference point (initial position)
  // to cylinder containing DV x,y, and then reports
  // the updated phi coordinate the the muon at that point
  float B = 3.8; // b field in Tesla
  float mass = 0.10566; // muon mass in GeV
  float c = 0.29979; // speed of light in m/ns

  float P = pow(px*px + py*py + pz*pz, 0.5);
  float Pxy = pow(px*px + py*py, 0.5);
  float E = pow(P*P + mass*mass, 0.5);

  // relativistic velocities
  float vx = px/E * c;
  float vy = py/E * c;
  float vz = pz/E * c;
  float vxy = pow(vx*vx + vy*vy, 0.5);

  // larmor radius/angular frequency
  float R = 1e3*Pxy/(3.*charge*B);
  float w = vxy/R;

  float t = 0.;
  float curr_x = refx;
  float curr_y = refy;
  float current_rho = -1;
  float target_rho = pow(dvx*dvx + dvy*dvy, 0.5);
  float dt = 0.1;
  int nsteps = 0;
  while(current_rho < target_rho) {
      curr_x = refx + (vy/w)*( 1-TMath::Cos(w*t)) + (vx/w)*TMath::Sin(w*t);
      curr_y = refy + (vx/w)*(-1+TMath::Cos(w*t)) + (vy/w)*TMath::Sin(w*t);
      current_rho = pow(curr_x*curr_x + curr_y*curr_y, 0.5);
      t += dt;
      nsteps++;
      if (nsteps > 10000) {
          //std::cout << "Warning, >10000 steps in propagate_to_cylinder" << std::endl;
          break;
      }
  }
  float curr_vx = vy*TMath::Sin(w*t) + vx*TMath::Cos(w*t);
  float curr_vy = vy*TMath::Cos(w*t) - vx*TMath::Sin(w*t);
  float newphi = atan2(curr_vy, curr_vx);
  return newphi;

}

// -----------------------
bool GenTauDauPropagator::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
  {
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}


ReferenceCountingPointer<BoundCylinder> GenTauDauPropagator::theBarrel_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenTauDauPropagator::thePositiveEndcap_ = nullptr;
ReferenceCountingPointer<BoundDisk> GenTauDauPropagator::theNegativeEndcap_ = nullptr;

DEFINE_FWK_MODULE(GenTauDauPropagator);
