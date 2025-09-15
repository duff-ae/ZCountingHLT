import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register("inFile", "file:", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Path to the input file")
options.register("outFile", "file:/tmp/out.root", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Path of the output file")
options.parseArguments()

process = cms.Process("Ntuple")

# --- L1 Unpacker for Scouting Data ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag("hltFEDSelectorL1")  # Typical for Scouting

process.load("FWCore.MessageService.MessageLogger_cfi")

from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = autoCond['run3_data']

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFile)
)

process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inFile)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

L1Info = ["L1_SingleMu0_BMTF", "L1_SingleMu3", "L1_SingleMu5", "L1_SingleMu5_BMTF"]

process.ntuplizer = cms.EDAnalyzer("ScoutingMuonNtuplizer",
    muonCollection = cms.InputTag("hltScoutingMuonPackerVtx"),
    primaryVertexCollection = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    muonVertexCollection = cms.InputTag("hltScoutingMuonPackerVtx", "displacedVtx"),
    pfMetPhiValue = cms.InputTag("hltScoutingPFPacker", "pfMetPhi"),
    pfMetPtValue = cms.InputTag("hltScoutingPFPacker", "pfMetPt"),
    l1Seeds           = cms.vstring(L1Info),
    AlgInputTag       = cms.InputTag("gtStage2Digis"),
    l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
    l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
    doL1              = cms.bool(True),
)

process.p = cms.Path(process.ntuplizer)

process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
