import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register("inFile",  "", VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string, "Input file")
options.register("outFile", "ntuple.root", VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string, "Output ROOT file")
options.parseArguments()

process = cms.Process("ANA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inFile)
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFile)
)

process.scoutingW = cms.EDAnalyzer("ScoutingWAnalyzer",
    fid_mu_pt_min = cms.double(25.0),
    fid_mu_abs_eta_max = cms.double(2.4),

    pv_min_ndof = cms.int32(4),
    pv_max_abs_z = cms.double(24.0),
    pv_max_rho = cms.double(2.0),

    trk_max_normchi2 = cms.double(10.0),
    trk_max_abs_dz = cms.double(0.5),
    trk_max_abs_dxy = cms.double(0.2),
    trk_max_lostInnerHits = cms.int32(2),

    puppi_drmax = cms.double(0.4),
    puppi_eps = cms.double(0.02),
    puppi_max_neighbors = cms.int32(64),
    puppi_beta_pow = cms.double(1.0),

    scoutingMuons = cms.InputTag("hltScoutingMuonPacker"),
    scoutingVertices = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    scoutingParticles = cms.InputTag("hltScoutingPFPacker"),

    pfMetPtValue  = cms.InputTag("hltScoutingPFPacker", "pfMetPt"),
    pfMetPhiValue = cms.InputTag("hltScoutingPFPacker", "pfMetPhi"),
    rhoValue      = cms.InputTag("hltScoutingPFPacker", "rho"),

    pileup       = cms.InputTag("slimmedAddPileupInfo"),
    genParticles = cms.InputTag("prunedGenParticles"),
    genMetTrue   = cms.InputTag("genMetTrue", "", "HLT"),
)

process.p = cms.Path(process.scoutingW)