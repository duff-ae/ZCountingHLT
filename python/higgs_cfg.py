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

process.jpsi = cms.EDAnalyzer("ScoutingJpsiHiggs",
    muonCollection          = cms.InputTag("hltScoutingMuonPackerVtx"),
    primaryVertexCollection = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    muonVertexCollection    = cms.InputTag("hltScoutingMuonPackerVtx", "displacedVtx"),
    photonCollection        = cms.InputTag("hltScoutingEgammaPacker"),
    pfMetPhiValue           = cms.InputTag("hltScoutingPFPacker", "pfMetPhi"),
    pfMetPtValue            = cms.InputTag("hltScoutingPFPacker", "pfMetPt"),

    # окна/пороги
    jMin = cms.double(2.6), jMax = cms.double(3.6),
    hMin = cms.double(115.0), hMax = cms.double(135.0),

    minPtForJpsiMuon = cms.double(3.0),
    maxAbsEta        = cms.double(2.4),
    minDeltaR        = cms.double(0.02),
    sigmaJpsi        = cms.double(0.045),

    doPhotons        = cms.untracked.bool(True),
    minPtPhoton      = cms.double(25.0),
    maxAbsEtaPhoton  = cms.double(2.5),
    minDRMuGamma     = cms.double(0.2),
)

process.p = cms.Path(process.jpsi)

process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
