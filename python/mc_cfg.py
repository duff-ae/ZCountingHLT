import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register("inFile", "file:", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Path to the input file")
options.register("outFile", "file:/tmp/out.root", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Path of the output file")
options.parseArguments()

process = cms.Process("SCOUTMC")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFile)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inFile)
)

process.demo = cms.EDAnalyzer('ScoutingMCJpsiAnalyzer')
process.p = cms.Path(process.demo)