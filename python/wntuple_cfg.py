# run_scoutingW_cfg.py
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register("inFile",  "", VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string, "Input file")
options.register("outFile", "ntuple.root", VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string, "Output ROOT file")
options.parseArguments()

process = cms.Process("Ntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inFile)
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFile)
)

process.ntuplizer = cms.EDAnalyzer(
    "ScoutingWNtuplizer",
    fid_mu_pt_min      = cms.double(3.0),
    fid_mu_abs_eta_max = cms.double(2.4),
    dr_match_tight     = cms.double(0.10),
    dr_match_loose     = cms.double(0.20),

    scoutingMuons = cms.InputTag("hltScoutingMuonPacker"),
    scoutingPrimaryVertices = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    pfMetPtValue  = cms.InputTag("hltScoutingPFPacker", "pfMetPt"),
    pfMetPhiValue = cms.InputTag("hltScoutingPFPacker", "pfMetPhi"),
    rhoValue      = cms.InputTag("hltScoutingPFPacker", "rho"),

    genEventInfo = cms.InputTag("generator"),
    pileup       = cms.InputTag("slimmedAddPileupInfo"),
    genParticles = cms.InputTag("prunedGenParticles"),
    genMetTrue   = cms.InputTag("genMetTrue"),
)
process.p = cms.Path(process.ntuplizer)