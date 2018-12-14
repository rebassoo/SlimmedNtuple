import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:step3_fpmc_MiniAOD.root'
        'file:/home/users/rebassoo/work/2018_12_13_TestingUpdatedProtonReco/CMSSW_9_4_11/src/Validation/CTPPS/test_2017/acceptance_test/ctppsSim.root'
        )
)

process.load("SlimmedNtuple.Ntupler.CfiFile_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
MC=True

#Global tags from here:https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
if MC:
    process.GlobalTag.globaltag ='94X_mc2017_realistic_v14'
else:
    process.GlobalTag.globaltag ='94X_dataRun2_v6'

process.demo = cms.EDAnalyzer('Ntupler')

process.p = cms.Path(process.demo)












