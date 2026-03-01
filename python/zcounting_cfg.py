import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# ---------------------------
# CLI
# ---------------------------
options = VarParsing.VarParsing('analysis')
options.register("inFile",
                 "file:",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to the input file")
options.register("outFile",
                 "file:/tmp/zscout.root",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path of the output file")
options.parseArguments()

process = cms.Process("ZCount")

# ---------------------------
# Services / Conditions / Logging
# ---------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = autoCond['run3_data']

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outFile)
)

# ---------------------------
# Source
# ---------------------------
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inFile)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)   # при тестах можно заменить на 1000 и т.п.
)

# ---------------------------
# Collections (проверь по edmDump имена инстансов под твой файл!)
# ---------------------------
muon_tag      = cms.InputTag("hltScoutingMuonPackerVtx")                     # Run3ScoutingMuon
pv_tag        = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx")  # Run3ScoutingVertex
track_tag     = cms.InputTag("hltScoutingTrackPacker")                       # Run3ScoutingTrack
muon_vtx_tag  = cms.InputTag("hltScoutingMuonPackerVtx", "displacedVtx")     # Run3ScoutingVertex

# ---------------------------
# Analyzer
# ---------------------------
process.zscout = cms.EDAnalyzer(
    "ScoutingZCounting",

    # --- входные коллекции ---
    muonCollection          = muon_tag,
    primaryVertexCollection = pv_tag,
    trackCollection         = track_tag,
    muonVertexCollection    = muon_vtx_tag,     # для J/ψ prompt/non-prompt по Lxy

    # --- включатели подсекций ---
    doZ      = cms.untracked.bool(True),
    doJpsi   = cms.untracked.bool(True),
    doW      = cms.untracked.bool(True),

    # --- Probe kinematics для Z T&P ---
    PtCutL2  = cms.untracked.double(15.0),
    EtaCutL2 = cms.untracked.double(2.4),

    # --- Z mass window ---
    MassBin  = cms.untracked.int32(60),
    MassMin  = cms.untracked.double(60.0),
    MassMax  = cms.untracked.double(120.0),

    # --- PV histogram binning ---
    PVBin    = cms.untracked.int32(100),
    PVMin    = cms.untracked.double(0.0),
    PVMax    = cms.untracked.double(100.0),

    # --- PV selection ---
    VtxNTracksFitMin = cms.untracked.double(0.0),
    VtxNdofMin       = cms.untracked.double(4.0),
    VtxAbsZMax       = cms.untracked.double(24.0),
    VtxRhoMax        = cms.untracked.double(2.0),

    # --- "HLT" proxy: IsoMu24-like (используется только как tag-прокси) ---
    #   hltIsoRelCut — максимум относительной изоляции (trk+ECAL+HCAL)/pt
    hltIsoRelCut     = cms.untracked.double(0.15),

    # --- Muon quality proxies (аналог passGlobalMuon / CustomTight) ---
    minTrackerLayers = cms.untracked.int32(5),
    minPixelLayers   = cms.untracked.int32(1),
    maxGlobalChi2    = cms.untracked.double(20.0),

    # --- MET (для W): проверь названия в edmDump ---
    pfMetPtValue     = cms.untracked.InputTag("hltScoutingPFPacker", "pfMetPt"),
    pfMetPhiValue    = cms.untracked.InputTag("hltScoutingPFPacker", "pfMetPhi"),

    # --- прочее ---
    runNumber        = cms.untracked.int32(-1),
    saveIntermediate = cms.untracked.bool(True)
)

# ---------------------------
# Path
# ---------------------------
process.p = cms.Path(process.zscout)

process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    wantSummary   = cms.untracked.bool(True)
)