import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# ---------------------------
# CLI options
# ---------------------------
options = VarParsing.VarParsing('analysis')
options.register("inFile",
                 "file:",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to the input file")
options.register("outFile",
                 "file:/tmp/out.root",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path of the output file")
options.register("useNoVtx",
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Use hltScoutingMuonPackerNoVtx instead of ...Vtx")
options.parseArguments()

process = cms.Process("ZCount")

# ---------------------------
# Conditions / logging
# ---------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = autoCond['run3_data']

# ---------------------------
# I/O
# ---------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inFile)
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFile)
)

# (optional) L1 unpacker often used with Scouting
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag("hltFEDSelectorL1")

# ---------------------------
# Collections (per your edmDump)
# ---------------------------
muon_tag = cms.InputTag("hltScoutingMuonPackerVtx")
if options.useNoVtx:
    muon_tag = cms.InputTag("hltScoutingMuonPackerNoVtx")

pv_tag    = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx")
track_tag = cms.InputTag("hltScoutingTrackPacker")
trg_tag   = cms.InputTag("TriggerResults","","HLT")

# ---------------------------
# Analyzer: ScoutingZCounting
# ---------------------------
MuonHLT = [
    "^HLT_IsoMu24_v\\d+$",
    "^HLT_Mu50_v\\d+$",
]

process.zscout = cms.EDAnalyzer("ScoutingZCounting",
    muonCollection          = muon_tag,
    primaryVertexCollection = pv_tag,
    trackCollection         = track_tag,                 # <-- добавлено
    triggerResults          = trg_tag,
    MuonTriggerNames        = cms.vstring(MuonHLT),

    # --- Kinematic cuts (tag/probe) ---
    PtCutL1 = cms.untracked.double(24.0),
    PtCutL2 = cms.untracked.double(15.0),
    EtaCutL1 = cms.untracked.double(2.4),
    EtaCutL2 = cms.untracked.double(2.4),

    # --- Mass window for Z ---
    MassBin = cms.untracked.int32(60),
    MassMin = cms.untracked.double(60.0),
    MassMax = cms.untracked.double(120.0),

    # --- PV histogram binning ---
    PVBin = cms.untracked.int32(100),
    PVMin = cms.untracked.double(0.0),
    PVMax = cms.untracked.double(100.0),

    # --- PV selection (как в ZCounting) ---
    VtxNTracksFitMin = cms.untracked.double(0.0),
    VtxNdofMin       = cms.untracked.double(4.0),
    VtxAbsZMax       = cms.untracked.double(24.0),
    VtxRhoMax        = cms.untracked.double(2.0),

    # --- ID / ISO ---
    IDType  = cms.untracked.string("CustomTight"),   # Loose/Medium/Tight/CustomTight/NoneID
    IsoType = cms.untracked.string("PF-based"),      # Tracker-based/PF-based/NoneIso
    IsoCut  = cms.untracked.double(10.0)
)

process.p = cms.Path(process.zscout)

process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    wantSummary   = cms.untracked.bool(True)
)