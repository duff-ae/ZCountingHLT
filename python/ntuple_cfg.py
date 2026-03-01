import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# --- CLI options ---
options = VarParsing.VarParsing('analysis')
options.register("inFile",  "file:",                 VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Path to the input file")
options.register("outFile", "file:/tmp/ntuple.root", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Path of the output file")
options.parseArguments()

process = cms.Process("Ntuple")

# --- L1 Unpacker for Scouting Data (можно и убрать, ntuplizerу не обязателен) ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag("hltFEDSelectorL1")

# --- Logging ---
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# --- Conditions / GlobalTag ---
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = autoCond['run3_data']

# --- Output ROOT file ---
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFile)
)

# --- Input ---
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inFile)
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

# --- Our ntuplizer: full scouting content ---
process.ntuplizer = cms.EDAnalyzer("ScoutingNtuplizer",
    # muons & tracks
    muonCollection  = cms.InputTag("hltScoutingMuonPackerVtx"),    # vector<Run3ScoutingMuon>
    trackCollection = cms.InputTag("hltScoutingTrackPacker"),      # vector<Run3ScoutingTrack>

    # electrons & photons (оба из EgammaPacker)
    electronCollection = cms.InputTag("hltScoutingEgammaPacker"),  # vector<Run3ScoutingElectron>
    photonCollection   = cms.InputTag("hltScoutingEgammaPacker"),  # vector<Run3ScoutingPhoton>

    # PF jets & PF particles
    pfJetCollection       = cms.InputTag("hltScoutingPFPacker"),   # vector<Run3ScoutingPFJet>
    pfCandidateCollection = cms.InputTag("hltScoutingPFPacker"),   # vector<Run3ScoutingParticle>

    # vertices
    primaryVertexCollection        = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    displacedVertexCollectionNoVtx = cms.InputTag("hltScoutingMuonPackerNoVtx", "displacedVtx"),
    displacedVertexCollectionVtx   = cms.InputTag("hltScoutingMuonPackerVtx",   "displacedVtx"),

    # flags
    storeMuons     = cms.untracked.bool(False),
    storeTracks    = cms.untracked.bool(False),
    storeElectrons = cms.untracked.bool(False),
    storePhotons   = cms.untracked.bool(False),
    storePFJets    = cms.untracked.bool(False),
    storePFCands   = cms.untracked.bool(True),
    storeVertices  = cms.untracked.bool(False),
    storeMET       = cms.untracked.bool(True),
    storeHLT       = cms.untracked.bool(True),
)

process.p = cms.Path(process.ntuplizer)

process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    wantSummary   = cms.untracked.bool(True)
)