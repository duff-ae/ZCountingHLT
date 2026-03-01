import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# То, чего *нет* по умолчанию — можно регистрировать
options.register(
    "inFile",
    "",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Path to the input file"
)

options.register(
    "outFile",
    "file:/tmp/jpsi_scout_mc_ntuple.root",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Path of the output file"
)

# maxEvents уже существует, просто при желании меняем default:
# options.maxEvents = -1  # all events
# или оставляем дефолт VarParsing-а

options.parseArguments()

if options.inFile == "":
    raise RuntimeError("Please provide --inFile=...")

process = cms.Process("SCOUTMC")

process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string(options.outFile)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inFile)
)

process.scout_mc = cms.EDAnalyzer(
    "ScoutingMCJpsiNtuplizer",
    genParticles      = cms.InputTag("prunedGenParticles"),
    pileup            = cms.InputTag("slimmedAddPileupInfo"),
    beamSpot          = cms.InputTag("offlineBeamSpot"),

    scoutingMuons     = cms.InputTag("hltScoutingMuonPacker"),
    scoutingTracks    = cms.InputTag("hltScoutingTrackPacker"),
    scoutingVertices  = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    scoutingElectrons = cms.InputTag("hltScoutingEgammaPacker", "electrons"),
    scoutingPhotons   = cms.InputTag("hltScoutingEgammaPacker", "photons"),
    scoutingPFJets    = cms.InputTag("hltScoutingPFPacker", "pfJets"),
    scoutingParticles = cms.InputTag("hltScoutingPFPacker", "pfCandidates"),

    patMuons          = cms.InputTag("slimmedMuons")
)

process.p = cms.Path(process.scout_mc)