import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register(
    "inFile",
    "file:/tmp/in.root",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Path to the input file"
)
options.register(
    "outFile",
    "/tmp/out.root",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Path of the output file"
)
options.parseArguments()

process = cms.Process("RHO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.cerr.FwkReport.reportEvery = 10

from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = autoCond["run3_data"]

process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outFile))

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inFile))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.rho = cms.EDAnalyzer(
    "ScoutingHiggsRho",

    # --- REQUIRED tags (do not rename collections) ---
    photonCollection = cms.InputTag("hltScoutingEgammaPacker"),
    pfCandCollection = cms.InputTag("hltScoutingPFPacker"),
    pfJetCollection  = cms.InputTag("hltScoutingPFPacker"),

    # --- Control ---
    storeHistograms = cms.untracked.bool(True),

    # MC/data switch for the module (keep even if you run data now)
    isMC = cms.untracked.bool(True),

    # --- GEN inputs (kept for MC; harmless for data with ProductNotFound handling) ---
    pileupSummary      = cms.untracked.InputTag("slimmedAddPileupInfo", "", "PAT"),
    prunedGenParticles = cms.untracked.InputTag("prunedGenParticles", "", "PAT"),
    packedGenParticles = cms.untracked.InputTag("packedGenParticles", "", "PAT"),

    # --- GEN fiducial knobs ---
    fidPhoMaxAbsEta = cms.untracked.double(2.5),
    fidPiMaxAbsEta  = cms.untracked.double(2.5),
    genPhoMinPt     = cms.untracked.double(0.0),
    genPiMinPt      = cms.untracked.double(0.0),

    # --- Scouting photon selection ---
    phoMinPt     = cms.untracked.double(0.0),
    phoMaxAbsEta = cms.untracked.double(0.0),

    # --- Scouting PF selection ---
    pfcMinPt     = cms.untracked.double(1.0),
    pfcMaxAbsEta = cms.untracked.double(2.4),
    pfcMaxAbsDz  = cms.untracked.double(0.2),
    pfcMaxAbsDxy = cms.untracked.double(0.1),
    pfPionOnly   = cms.untracked.bool(True),
    pfUseTrkVars = cms.untracked.bool(False),

    # --- “best-candidate” cuts (pf-level) ---
    pfRequireSameVtx = cms.untracked.bool(True),
    pfMaxAbsDzDiff   = cms.untracked.double(0.05),
    pfMaxAbsDxyDiff  = cms.untracked.double(0.05),
    pfMaxDR2PiPi     = cms.untracked.double(0.25),
    pfRhoMinPt       = cms.untracked.double(15.0),
    pfMinAbsDPhiGRho = cms.untracked.double(2.7),
    pfMinPtBal       = cms.untracked.double(0.5),
    pfMaxPtBal       = cms.untracked.double(2.0),
    pfRhoMassMin     = cms.untracked.double(0.0),
    pfRhoMassMax     = cms.untracked.double(0.0),

    # --- Truth matching ---
    doTruthMatching = cms.untracked.bool(True),
    matchPhoDR2     = cms.untracked.double(0.01),
    matchPionDR2    = cms.untracked.double(0.01),

    # --- Scouting jets (kept; GEN jets were removed earlier, but scouting jets are still used/kept) ---
    jetMinPt     = cms.untracked.double(0.0),
    htJetMinPt   = cms.untracked.double(0.0),
    jetMaxAbsEta = cms.untracked.double(4.7),

    # --- L1 Stage-2 (RECO) ---
    l1EGamma = cms.untracked.InputTag("caloStage2Digis", "EGamma"),
    l1Jet    = cms.untracked.InputTag("caloStage2Digis", "Jet"),

    doL1 = cms.untracked.bool(True),
    l1UseBX0Only = cms.untracked.bool(True),

    l1EGMinPt = cms.untracked.double(0.0),
    l1JetMinPt = cms.untracked.double(0.0),

    matchL1EG_DR2  = cms.untracked.double(0.04),
    matchL1Jet_DR2 = cms.untracked.double(0.09),
)

process.p = cms.Path(process.rho)

process.options = cms.untracked.PSet(
    TryToContinue=cms.untracked.vstring('ProductNotFound'),
    wantSummary=cms.untracked.bool(True)
)