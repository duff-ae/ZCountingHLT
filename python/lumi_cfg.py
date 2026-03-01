import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register("inFile", "file:", VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string, "Input file")
options.register("outFile", "file:/tmp/scouting_sc.root", VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string, "Output ROOT file")
options.parseArguments()

process = cms.Process("SCOUTSC")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = autoCond["run3_data"]

process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outFile))

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inFile))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

# Collections (switch here to compare both collections in separate runs)
muon_tag     = cms.InputTag("hltScoutingMuonPackerVtx")                      # Run3ScoutingMuon
#muon_tag     = cms.InputTag("hltScoutingMuonPackerNoVtx")                   # Run3ScoutingMuon

# For J/psi prompt/nonprompt classification (optional in code; recommended if available)
muon_vtx_tag = cms.InputTag("hltScoutingMuonPackerVtx","displacedVtx")       # Run3ScoutingVertex
#muon_vtx_tag = cms.InputTag("hltScoutingMuonPackerNoVtx","displacedVtx")    # Run3ScoutingVertex

pv_tag       = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx")   # Run3ScoutingVertex
track_tag    = cms.InputTag("hltScoutingTrackPacker")                        # Run3ScoutingTrack

process.sc = cms.EDAnalyzer(
    "ScoutingZ",
    muonCollection          = muon_tag,
    primaryVertexCollection = pv_tag,
    trackCollection         = track_tag,

    # Optional: only needed for J/psi prompt/nonprompt split; code handles absence safely
    muonVertexCollection    = muon_vtx_tag,

    doZ    = cms.untracked.bool(True),
    doJpsi = cms.untracked.bool(True),

    # Tag quality (GLOBAL NOT required; measured separately)
    tagEtaMax           = cms.untracked.double(2.4),
    tagTrkChi2NdofMax   = cms.untracked.double(3.0),
    tagMinTrackerLayers = cms.untracked.int32(4),
    tagMinPixelLayers   = cms.untracked.int32(1),

    # Probe base (inclusive)
    probePtMin = cms.untracked.double(3.0),
    probeEtaMax= cms.untracked.double(2.4),

    # Measured global component
    gloNormChi2Max = cms.untracked.double(999.0),

    # Z tag pT: main + scan
    zTagPtMain = cms.untracked.double(15.0),
    zTagPtScan = cms.untracked.vdouble(3.0, 10.0, 15.0, 25.0),

    # Track probe
    trkPtMin        = cms.untracked.double(3.0),
    trkEtaMax       = cms.untracked.double(2.4),
    trkChi2NdofMax  = cms.untracked.double(3.0),
    trkMinLayers    = cms.untracked.int32(4),
    dr2Match        = cms.untracked.double(0.09),

    # Z mass windows
    zMassBins     = cms.untracked.int32(60),
    zMassMin      = cms.untracked.double(60.0),
    zMassMax      = cms.untracked.double(120.0),
    zMassBinsCtrl = cms.untracked.int32(120),
    zMassMinCtrl  = cms.untracked.double(40.0),
    zMassMaxCtrl  = cms.untracked.double(160.0),

    # J/psi window + prompt split
    jpsiMassBins    = cms.untracked.int32(50),
    jpsiMassMin     = cms.untracked.double(2.6),
    jpsiMassMax     = cms.untracked.double(3.6),
    lxySigPromptMax = cms.untracked.double(3.0),

    # J/psi candidate limiting (to reduce combinatorics but allow >1 real pair)
    jpsiTopN              = cms.untracked.uint32(2),
    jpsiUseDisjointPairs  = cms.untracked.bool(True),

    # PV selection (for good-PV counting used in buffers; control only)
    vtxNdofMin  = cms.untracked.double(4.0),
    vtxAbsZMax  = cms.untracked.double(24.0),
    vtxRhoMax   = cms.untracked.double(2.0),
)

process.p = cms.Path(process.sc)

process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring("ProductNotFound"),
    wantSummary   = cms.untracked.bool(True),
)