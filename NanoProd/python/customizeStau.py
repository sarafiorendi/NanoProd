import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.NanoAOD.nano_eras_cff import *

from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer
## for 14_0_X
from PhysicsTools.NanoAOD.simplePATJetFlatTableProducer_cfi import simplePATJetFlatTableProducer
from PhysicsTools.NanoAOD.simplePATMuonFlatTableProducer_cfi import simplePATMuonFlatTableProducer

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import *
from PhysicsTools.NanoAOD.jetsAK4_Puppi_cff import *

from PhysicsTools.PatAlgos.tools.puppiJetMETReclusteringTools import puppiAK4METReclusterFromMiniAOD

import os


def customise_run3_jets(process):
    process.linkedObjects.jets="finalJets" 
    
    process.nanoTableTaskCommon.add(process.jetTask)
    process.nanoTableTaskCommon.add(process.jetForMETTask)
    ## remove puppi table otherwise it tries to save the jet table twice
    process.jetPuppiForMETTask.remove(process.corrT1METJetPuppiTable)
    
    process.jetTablesTask.remove(process.bjetNN)
    process.jetTablesTask.remove(process.cjetNN)
    process.nanoTableTaskCommon.replace(process.jetPuppiTablesTask, process.jetTablesTask)
    
    process.jetTable.externalVariables = cms.PSet()
    
    process.ptRatioRelForEle.srcJet="updatedJets"
    process.ptRatioRelForMu.srcJet="updatedJets"

    return process
 

def customize_process_and_associate(process, isMC, useCHSJets = True) :
    # Lost tracks
#     process.lostTrackTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
#         src = cms.InputTag("lostTracks"),
#         #cut = cms.string(""),
#         cut = cms.string("pt > 1"),
#         name= cms.string("LostTrack"),
#         doc = cms.string("Lost tracks"),
#         singleton = cms.bool(False), # the number of entries is variable
#         extension = cms.bool(False), # this is the main table
#         variables = cms.PSet(
#             CandVars,
#         )
#     )

#     from PhysicsTools.NanoAOD.custom_muon_cff import AddVariablesForMuon
#     displacedMuonWithVariables = "displacedMuonWithVariables"
#     setattr(proc, displacedMuonWithVariables, cms.EDProducer("MuonSpecialVariables",
#                     muonSrc=cms.InputTag("slimmedDisplacedMuons"),
#                     vertexSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
#                     trkSrc=cms.InputTag("pfTracks"),
#                     )
#     )

    ## unpack PF candidates and lost track, for isolation
    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
    process.disMuonIsolation = cms.EDProducer(
          "DisplacedMuonIsolation",
          muons = cms.InputTag("slimmedDisplacedMuons", "", "PAT"),
          tracksForIso = cms.InputTag('unpackedTracksAndVertices'),
          beamSpot = cms.InputTag('offlineBeamSpot'),
    )
    d_displacedMuonIsoVars = {
          "tkIso_v0":     ExtVar("disMuonIsolation:score0"       , float, doc = "v0"),
          "tkIso_v1":     ExtVar("disMuonIsolation:score1"       , float, doc = "v1"),
    }
    
    process.disMuonTable = simplePATMuonFlatTableProducer.clone(
        src = cms.InputTag("slimmedDisplacedMuons"),
        name = cms.string("DisMuon"),
        doc = cms.string("Displaced Muon Collection"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(CandVars,
            ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
            tunepRelPt = Var("tunePMuonBestTrack().pt/pt",float,doc="TuneP relative pt, tunePpt/pt",precision=6),
            dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
            dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
            dxybs = Var("dB('BS2D')",float,doc="dxy (with sign) wrt the beam spot, in cm",precision=10),
            dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
            dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
            trkChi2 = Var("? globalTrack().isNonnull() ? globalTrack().normalizedChi2() : ? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from either globalTrack or innerTrack "),
            muonHits = Var("? globalTrack().isNonnull() ? globalTrack().hitPattern().numberOfValidMuonHits() : ?  innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidMuonHits() :-99",float,doc="Number of valid Muon Hits from either globalTrack or innerTrack"),
            pixelHits = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidPixelHits() : -99", float, doc="Numbr of valid pixel hits"),
            validFraction = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().validFraction() : -99", float, doc="Inner Track Valid Fraction"),
            positionChi2 = Var("combinedQuality().chi2LocalPosition", float, doc="chi2 Local Position"),
            trkKink = Var("combinedQuality().trkKink", float, doc="Track Kink"),
            ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
            sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
            segmentComp   = Var("segmentCompatibility()", float, doc = "muon segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
            nStations = Var("numberOfMatchedStations", "uint8", doc = "number of matched stations with default arbitration (segment & track)"),
            nTrackerLayers = Var("?track.isNonnull?innerTrack().hitPattern().trackerLayersWithMeasurement():0", "uint8", doc = "number of layers in the tracker"),
            highPurity = Var("?track.isNonnull?innerTrack().quality('highPurity'):0", bool, doc = "inner track is high purity"),
            jetIdx = Var("?hasUserCand('jet')?userCand('jet').key():-1", "int16", doc="index of the associated jet (-1 if none)"),
            svIdx = Var("?hasUserCand('vertex')?userCand('vertex').key():-1", "int16", doc="index of matching secondary vertex"),
            tkRelIso = Var("isolationR03().sumPt/tunePMuonBestTrack().pt",float,doc="Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt",precision=6),
            #miniPFRelIso_chg = Var("userFloat('miniIsoChg')/pt",float,doc="mini PF relative isolation, charged component"),
            #miniPFRelIso_all = Var("userFloat('miniIsoAll')/pt",float,doc="mini PF relative isolation, total (with scaled rho*EA PU corrections)"),
            pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
            pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
            pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
            #jetRelIso = Var("?userCand('jetForLepJetVar').isNonnull()?(1./userFloat('ptRatio'))-1.:(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)",precision=8),
            #jetPtRelv2 = Var("?userCand('jetForLepJetVar').isNonnull()?userFloat('ptRel'):0",float,doc="Relative momentum of the lepton with respect to the closest jet after subtracting the lepton",precision=8),
            tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0", "uint8", doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
            looseId  = Var("isGlobalMuon||isTrackerMuon",bool, doc="muon is loose muon"),
            isPFcand = Var("isPFMuon",bool,doc="muon is PF candidate"),
            isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon"),
            isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon"),
            isStandalone = Var("isStandAloneMuon",bool,doc="muon is a standalone muon"),
            mediumId = Var("? (globalTrack().isNonnull()) && (isGlobalMuon) && (globalTrack().normalizedChi2() < 3) && (combinedQuality().chi2LocalPosition < 12) && (combinedQuality().trkKink < 20) && (innerTrack().validFraction() > 0.8) && (segmentCompatibility() > 0.303) ? 1 : ? (globalTrack().isNonnull()) && (innerTrack().validFraction() > 0.8) && (segmentCompatibility() > 0.451) ? 1 : 0",bool,doc="cut-based ID, medium WP"),
            mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
            tightId = Var("? (globalTrack().isNonnull()) && (globalTrack().normalizedChi2() < 10) && (globalTrack().hitPattern().numberOfValidMuonHits() > 0) && (numberOfMatchedStations > 1) && (innerTrack().hitPattern().numberOfValidPixelHits() > 0) && (innerTrack().hitPattern().trackerLayersWithMeasurement() > 5) ? 1 : 0","bool",doc="cut-based ID, tight WP"),
            softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"),
            softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
            softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6),
            highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
            pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
            tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),
            miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
            mvaMuID = Var("mvaIDValue()",float,doc="MVA-based ID score ",precision=6),
            #mvaMuID_WP = Var("userFloat('mvaIDMuon_wpMedium') + userFloat('mvaIDMuon_wpTight')","uint8",doc="MVA-based ID selector WPs (1=MVAIDwpMedium,2=MVAIDwpTight)"),
            multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),
            puppiIsoId = Var("passed('PuppiIsoLoose')+passed('PuppiIsoMedium')+passed('PuppiIsoTight')", "uint8", doc="PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight)"),
            triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID"),
            inTimeMuon = Var("passed('InTimeMuon')",bool,doc="inTimeMuon ID"),

            #Sim Variables
#             simType = Var("? simType() ? simType() : -99",int,doc="simType"),
#             simExtType = Var("? simExtType() ? simExtType() : -99",int,doc="simExtType"),
#             simFlavour = Var("? simFlavour() ? simFlavour() : -99",int,doc="simFlavour"),
#             simHeaviestMotherFlavour = Var(" ? simHeaviestMotherFlavour() ? simHeaviestMotherFlavour() : -99",int,doc="simHeaviestMotherFlavour"),
#             simPdgId = Var("? simPdgId() ? simPdgId() : -99",int,doc="simPdgId"),
#             simMotherPdgId = Var("? simMotherPdgId() ? simMotherPdgId() : -99",int,doc="simMotherPdgId"),
# #             simBX = Var("? simBX() ? simBX() : -99",int,doc="simBX")
#             simProdRho = Var("? simProdRho() ? simProdRho(): -99",float,doc="simProdRho"),
#             simProdZ = Var("? simProdZ() ? simProdZ(): -99",float,doc="simProdZ"),
#             simPt = Var("? simPt() ? simPt(): -99",float,doc="simPt"),
#             simEta = Var("? simEta() ? simEta(): -99",float,doc="simEta"),
#             simPhi = Var("? simPhi() ? simPhi(): -99",float,doc='simPhi'),
        ),
    )

    process.disMuonTable.externalVariables = process.disMuonTable.externalVariables.clone(**d_displacedMuonIsoVars)
    process.disMuonTablesTask = cms.Task(process.disMuonTable)


    if isMC:     
        process.disMuonsMCMatchForTable = cms.EDProducer("MCMatcher",    # cut on deltaR, deltaPt/Pt; pick best by deltaR
            src         = process.disMuonTable.src,                      # final reco collection
            matched     = cms.InputTag("finalGenParticles"),             # final mc-truth particle collection
            mcPdgId     = cms.vint32(13),               # one or more PDG ID (13 = mu); absolute values (see below)
            checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
            mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
            maxDeltaR   = cms.double(0.3),              # Minimum deltaR for the match
            maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
            resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
            resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
        )
    
        process.disMuonMCTable = cms.EDProducer("CandMCMatchTableProducer",
            src     = process.disMuonTable.src,
            mcMap   = cms.InputTag("disMuonsMCMatchForTable"),
            objName = process.disMuonTable.name,
            objType = cms.string("Muon"), #cms.string("Muon"),
            branchName = cms.string("genPart"),
            docString = cms.string("MC matching to status==1 muons"),
        )
        
        process.disMuonMCTask = cms.Task(process.disMuonsMCMatchForTable, process.disMuonMCTable)  


    muonTableForID = muonTable.clone()
    muonTableForID.variables.muonHits = Var("? globalTrack().isNonnull() ? globalTrack().hitPattern().numberOfValidMuonHits() : ?  innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidMuonHits() :-99",float,doc="Number of valid Muon Hits from either globalTrack or innerTrack")
    muonTableForID.variables.pixelHits = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidPixelHits() : -99", float, doc="Numbr of valid pixel hits")
    muonTableForID.variables.trkChi2 = Var("? globalTrack().isNonnull() ? globalTrack().normalizedChi2() : ? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from either globalTrack or innerTrack ")
    muonTableForID.variables.positionChi2 = Var("combinedQuality().chi2LocalPosition", float, doc="chi2 Local Position")
    muonTableForID.variables.trkKink = Var("combinedQuality().trkKink", float, doc="Track Kink")
#     muonTableForID.variables.simType = Var("? simType() ? simType() : -99",int,doc="simType")
#     muonTableForID.variables.simExtType = Var("? simExtType() ? simExtType() : -99",int,doc="simExtType")
#     muonTableForID.variables.simFlavour = Var("? simFlavour() ? simFlavour() : -99",int,doc="simFlavour")
#     muonTableForID.variables.simHeaviestMotherFlavour = Var(" ? simHeaviestMotherFlavour() ? simHeaviestMotherFlavour() : -99",int,doc="simHeaviestMotherFlavour")
#     muonTableForID.variables.simPdgId = Var("? simPdgId() ? simPdgId() : -99",int,doc="simPdgId")
#     muonTableForID.variables.simMotherPdgId = Var("? simMotherPdgId() ? simMotherPdgId() : -99",int,doc="simMotherPdgId")
#     muonTableForID.variables.simBX = Var("? simBX() ? simBX() : -99",int,doc="simBX")
#     muonTableForID.variables.simProdRho = Var("? simProdRho() ? simProdRho(): -99",float,doc="simProdRho")
#     muonTableForID.variables.simProdZ = Var("? simProdZ() ? simProdZ(): -99",float,doc="simProdZ")
#     muonTableForID.variables.simPt = Var("? simPt() ? simPt(): -99",float,doc="simPt")
#     muonTableForID.variables.simEta = Var("? simEta() ? simEta(): -99",float,doc="simEta")
#     muonTableForID.variables.simPhi = Var("? simPhi() ? simPhi(): -99",float,doc='simPhi')

    process.globalReplace("muonTable", muonTableForID)

    if isMC:     
        # GenParticles
#         genParticleTable.variables.vertexX        = Var("vertex.X"      , float)
#         genParticleTable.variables.vertexY        = Var("vertex.Y"      , float)
#         genParticleTable.variables.vertexZ        = Var("vertex.Z"      , float)
        genParticleTable.variables.vertexRho      = Var("vertex.Rho"    , float)
        genParticleTable.variables.vertexR        = Var("vertex.R"      , float)
        
    file = "NanoProd/data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"
    if os.path.exists(file):
      file_string = file
    elif os.path.exists( os.path.basename(file) ):
      file_string = os.path.basename(file)
    else:
      file_string = "data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"

    process.disTauTag = cms.EDProducer(
          "DisTauTag",
          ## following line for crab
#           graphPath = cms.string(file_string),
          ## following line for local
          graphPath = cms.string("/afs/cern.ch/work/f/fiorendi/private/displacedTaus/desy/LLStaus_Run2/Production/data/models/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"),
###           graphPath = cms.string(os.getenv('CMSSW_BASE')+'/src/data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb'),
          jets = process.jetTable.src,
          pfCandidates = cms.InputTag('packedPFCandidates'),
          save_inputs  = cms.bool(False)
    )
        
    d_disTauTagVars = {
          "disTauTag_score0":     ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
          "disTauTag_score1":     ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
    }
    
    # Create the task
    if useCHSJets:
      process.jetTable.externalVariables = process.jetTable.externalVariables.clone(**d_disTauTagVars)
    ## for puppi jets, use this!
    else:
#       print ('adding disTau edproducer for PUPPI')
      process.jetPuppiTable.externalVariables = process.jetPuppiTable.externalVariables.clone(**d_disTauTagVars)
   
    process.custom_nanoaod_task = cms.Task(
       process.unpackedTracksAndVertices,
       process.disMuonIsolation,
       process.disMuonTablesTask,
       process.disTauTag,
    )
    
    # Associate the task 
    process.schedule.associate(process.custom_nanoaod_task)
    if isMC:
        process.custom_nanoaod_MC_task = cms.Task(
           process.disMuonMCTask,
        )
        process.schedule.associate(process.custom_nanoaod_MC_task)

#     process.disTauTask = cms.Task(process.disTauTag)
    return process


def BTVCustomNanoAODStaus(process, isMC):
    from PhysicsTools.NanoAOD.custom_btv_addpfcands_cff import addPFCands
#     addPFCands(process,True,False,False)  ## all PF Cands
    addPFCands(process,False,True,False) ## only AK4 cands
    
    ### for MC
    if isMC:
        process.load("PhysicsTools.NanoAOD.btvMC_cff")
        from PhysicsTools.NanoAOD.btvMC_cff import ak4onlyPFCandsMCSequence
        
        process.nanoSequenceMC+=ak4onlyPFCandsMCSequence
    
    return process


## add to the nanoAOD dedicated collections with the possible vertices between two leptons
def addDileptonVertices(process, isMC):

    process.dimuonVtxProducer = cms.EDProducer('DimuonVertex',
        src = cms.InputTag('slimmedMuons'),
        lep1Selection = cms.string('pt > 2.0 && abs(eta) < 2.4 '),
        lep2Selection = cms.string('pt > 2.0 && abs(eta) < 2.4 '),
        preVtxSelection  = cms.string('charge() == 0 && mass < 2.0'),
        postVtxSelection  = cms.string('userFloat("vtx_normChi2") > 0.001'),
        beamSpot = cms.InputTag('offlineBeamSpot')
#         preVtxSelection  = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1.'
#                                       '&& userFloat("lep_deltaR") > 0.03'),
#         postVtxSelection = cms.string('0 < userFloat("fitted_mass") && userFloat("fitted_mass") < 5.0'
#                                       '&& userFloat("sv_prob") > 0.001')
    )

    process.dimuonTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
        src = cms.InputTag("dimuonVtxProducer"),
        cut = cms.string(""),
        name = cms.string("dimuon"),
        doc  = cms.string("collection of dimuon vertices"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(CandVars,
#               sv_chi2 = Var("userFloat('sv_chi2')", float, doc="Vtx fit probability"),
              normChi2 = Var("userFloat('vtx_normChi2')", float, doc = "Vtx chi2 norm"),
              vtx_x = Var("userFloat('vtx_x')", float, doc = "Vtx position in x"),
              vtx_y = Var("userFloat('vtx_y')", float, doc = "Vtx position in y"),
              vtx_z = Var("userFloat('vtx_z')", float, doc = "Vtx position in y"),
              vtx_lxy = Var("userFloat('vtx_lxy')", float, doc = "Vtx 2d displacement from BS"),
              l1_pt  = Var("userCand('l1').pt", float, doc = "pt of lepton 1"),
              l1_eta = Var("userCand('l1').eta", float, doc = "eta of lepton 1"),
              l1_phi = Var("userCand('l1').phi", float, doc = "phi of lepton 1"),
              l1_charge = Var("userCand('l1').charge", int, doc = "charge of lepton 1"),
              l2_pt  = Var("userCand('l2').pt", float, doc = "pt of lepton 2"),
              l2_eta = Var("userCand('l2').eta", float, doc = "eta of lepton 2"),
              l2_phi = Var("userCand('l2').phi", float, doc = "phi of lepton 2"),
              l2_charge = Var("userCand('l2').charge", int, doc = "charge of lepton 2"),
        )
    )
    

    ## add vertices from displaced muons
    process.disDimuonVtxProducer = cms.EDProducer('DimuonVertex',
        src = cms.InputTag('slimmedDisplacedMuons'),
        lep1Selection = cms.string('pt > 2.0 && abs(eta) < 2.4 '),
        lep2Selection = cms.string('pt > 2.0 && abs(eta) < 2.4 '),
        preVtxSelection  = cms.string('charge() == 0 && mass < 2.0'),
        postVtxSelection  = cms.string('userFloat("vtx_normChi2") > 0.001'),
        beamSpot = cms.InputTag('offlineBeamSpot')
    )

    process.disDimuonTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
        src = cms.InputTag("disDimuonVtxProducer"),
        cut = cms.string(""),
        name = cms.string("disDimuon"),
        doc  = cms.string("collection of dimuon vertices from displaced muons"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(CandVars,
#               sv_chi2 = Var("userFloat('sv_chi2')", float, doc="Vtx fit probability"),
              normChi2 = Var("userFloat('vtx_normChi2')", float, doc = "Vtx chi2 norm"),
              vtx_x = Var("userFloat('vtx_x')", float, doc = "Vtx position in x"),
              vtx_y = Var("userFloat('vtx_y')", float, doc = "Vtx position in y"),
              vtx_z = Var("userFloat('vtx_z')", float, doc = "Vtx position in y"),
              vtx_lxy = Var("userFloat('vtx_lxy')", float, doc = "Vtx 2d displacement from BS"),
              l1_pt  = Var("userCand('l1').pt", float, doc = "pt of lepton 1"),
              l1_eta = Var("userCand('l1').eta", float, doc = "eta of lepton 1"),
              l1_phi = Var("userCand('l1').phi", float, doc = "phi of lepton 1"),
              l1_charge = Var("userCand('l1').charge", int, doc = "charge of lepton 1"),
              l2_pt  = Var("userCand('l2').pt", float, doc = "pt of lepton 2"),
              l2_eta = Var("userCand('l2').eta", float, doc = "eta of lepton 2"),
              l2_phi = Var("userCand('l2').phi", float, doc = "phi of lepton 2"),
              l2_charge = Var("userCand('l2').charge", int, doc = "charge of lepton 2"),
        )
    )


    ## add vertices from electrons
    process.dieleVtxProducer = cms.EDProducer('DielectronVertex',
        src = cms.InputTag('slimmedElectrons'),
        lep1Selection = cms.string('pt > 2.0 && abs(eta) < 2.4 '),
        lep2Selection = cms.string('pt > 2.0 && abs(eta) < 2.4 '),
        preVtxSelection  = cms.string('charge() == 0 && mass < 2.0'),
        postVtxSelection  = cms.string('userFloat("vtx_normChi2") > 0.001'),
        beamSpot = cms.InputTag('offlineBeamSpot')
    )


    process.dieleTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
        src = cms.InputTag("dieleVtxProducer"),
        cut = cms.string(""),
        name = cms.string("diele"),
        doc  = cms.string("collection of diele vertices"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(CandVars,
#               sv_chi2 = Var("userFloat('sv_chi2')", float, doc="Vtx fit probability"),
              normChi2 = Var("userFloat('vtx_normChi2')", float, doc = "Vtx chi2 norm"),
              vtx_x = Var("userFloat('vtx_x')", float, doc = "Vtx position in x"),
              vtx_y = Var("userFloat('vtx_y')", float, doc = "Vtx position in y"),
              vtx_z = Var("userFloat('vtx_z')", float, doc = "Vtx position in y"),
              vtx_lxy = Var("userFloat('vtx_lxy')", float, doc = "Vtx 2d displacement from BS"),
              l1_pt  = Var("userCand('l1').pt", float, doc = "pt of lepton 1"),
              l1_eta = Var("userCand('l1').eta", float, doc = "eta of lepton 1"),
              l1_phi = Var("userCand('l1').phi", float, doc = "phi of lepton 1"),
              l1_charge = Var("userCand('l1').charge", int, doc = "charge of lepton 1"),
              l2_pt  = Var("userCand('l2').pt", float, doc = "pt of lepton 2"),
              l2_eta = Var("userCand('l2').eta", float, doc = "eta of lepton 2"),
              l2_phi = Var("userCand('l2').phi", float, doc = "phi of lepton 2"),
              l2_charge = Var("userCand('l2').charge", int, doc = "charge of lepton 2"),
        )
    )


    process.custom_dilepton_task = cms.Task(
       process.dimuonVtxProducer,
       process.dimuonTable,
       process.disDimuonVtxProducer, 
       process.disDimuonTable,
       process.dieleVtxProducer,
       process.dieleTable,
    )

    # Associate the task 
    process.schedule.associate(process.custom_dilepton_task)
    return process




def customizeStau(process):

  isMC = True
  # customize stored objects

  ## for CHS
  process = customise_run3_jets(process)
  process = customize_process_and_associate(process, isMC, useCHSJets = True)
# #   ## for puppi tune v18
# # #   from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll as pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll
# # #   from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll as pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll
# # #   from RecoBTag.ONNXRuntime.pfUnifiedParticleTransformerAK4_cff import _pfUnifiedParticleTransformerAK4JetTagsAll as pfUnifiedParticleTransformerAK4JetTagsAll
# # #   
# # #   btagDiscriminatorsAK4 = cms.PSet(
# # #    names=cms.vstring(
# # #     'pfDeepFlavourJetTags:probb',
# # #     'pfDeepFlavourJetTags:probbb',
# # #     'pfDeepFlavourJetTags:problepb',
# # #     'pfDeepFlavourJetTags:probc',
# # #     'pfDeepFlavourJetTags:probuds',
# # #     'pfDeepFlavourJetTags:probg')
# # #     + pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll
# # #     + pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll
# # #     + pfUnifiedParticleTransformerAK4JetTagsAll
# # #   )
# # # 
# # #   process = puppiAK4METReclusterFromMiniAOD(process, runOnMC=True, useExistingWeights=False, btagDiscriminatorsAK4=btagDiscriminatorsAK4)
# #   ## end for puppi tune v18
# #   
# #   ## for any puppi
# # # #   process = customize_process_and_associate(process, isMC, useCHSJets = False)
# 
  ## btv custom
  process = BTVCustomNanoAODStaus(process, isMC)
  
  process.selectedFinalJetsConstituents = cms.EDFilter("PATPackedCandidatePtrSelector",
       src = cms.InputTag("finalJetsConstituentsTable"),
       cut = cms.string("pt > -1")
  )
  process.customizedPFCandsTask.add(process.selectedFinalJetsConstituents)
  process.customConstituentsExtTable.src="selectedFinalJetsConstituents"
  process.customAK8ConstituentsTable.candidates="selectedFinalJetsConstituents"
  process.customAK4ConstituentsTable.candidates="selectedFinalJetsConstituents"
  
  ## for CHS
  process.finalJetsAK4Constituents.src = cms.InputTag("finalJets")
  process.finalJetsAK4Constituents.cut = cms.string('(pt > 25) && (abs(eta) < 2.1)')
  process.customAK4ConstituentsTable.jets = cms.InputTag("finalJets")
  process.finalJets.cut = cms.string('(pt > 25) && (abs(eta) < 2.1)')
  
  ## add info on dilepton vertices, to study material interaction
  process = addDileptonVertices(process, isMC)

  
  return process
