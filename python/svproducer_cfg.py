import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import *
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

options.register('processName','Tree',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');
options.parseArguments()

process = cms.Process(options.processName,eras.Run2_2018)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff") 
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#verify which one to use
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'#'GR_E_V48'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:python/DYJets_M-50_2018.root'
        #'/store/mc/RunIISummer20UL18RECO/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/130000/016B6C98-203F-1447-B9B8-D2122D01D94A.root'
        #'/store/mc/RunIISummer20UL18RECO/GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/120000/372C073C-49E4-2141-B025-C4C1C98CB799.root'
    ))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

outfilename = 'testOutputFile.root'
options.register('outputFileName',outfilename,VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFileName))

process.timedVertices = cms.EDProducer("SVsTimingProducer",
                                       vertexLabel = cms.string("timedVertices"),
                                       electronSrc = cms.InputTag("gedGsfElectrons"),#"lowPtGsfElectrons"),
                                       timedVertices = cms.InputTag("timedVertices"),
)

process.p = cms.Path(process.timedVertices)
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(outfilename),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      "keep *_timedVertices_*_*", )
)


process.e = cms.EndPath(process.out)

