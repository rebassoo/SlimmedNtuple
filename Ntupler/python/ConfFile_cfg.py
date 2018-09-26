import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #        'file:/tmp/jjhollar/0C08DB9E-3543-E811-BB4D-0CC47A4D75F4.root'
        #        'file:/tmp/jjhollar/BulkGravToWW_narrow_M-1000_3A7945FC-1704-E811-9B66-008CFA110C5C.root'
        #        'file:/tmp/jjhollar/QCD_HT700to1000_F6F5131E-CCF9-E711-B15F-0242AC130002.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/BulkGravToWW_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/B8CED1AB-4D46-E811-B3F9-24BE05BDAE61.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/70000/FED622C0-7975-E811-AA65-0242AC130002.root'
        #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/36FDBA37-1C56-E811-B75E-00266CF2E388.root'
        #'root://xrootd.t2.ucsd.edu:2040//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B4005649-E955-E811-BE7B-0CC47A7C353E.root'
        'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B4005649-E955-E811-BE7B-0CC47A7C353E.root'
        #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/00469E05-E055-E811-807F-008CFAF7174A.root'
        #'file:step3_fpmc_MiniAOD.root'
        )
)

process.load("SlimmedNtuple.Ntupler.CfiFile_cfi")

process.load("SlimmedNtuple.Ntupler.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
MC=True
if MC:
    process.GlobalTag.globaltag ='94X_mc2017_realistic_v14'
else:
    process.GlobalTag.globaltag ='94X_dataRun2_v6'

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsAK8'),
   labelName = 'UpdatedJECAK8',
   jetCorrections = ('AK8PFchs', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual']), 'NONE')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual']), 'NONE')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)


from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.slimmedJetsAK8JetId = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                           filterParams = pfJetIDSelector.clone(),
                                           #src = cms.InputTag("slimmedJetsAK8"),
                                           src = cms.InputTag("updatedPatJetsUpdatedJECAK8"),
                                           filter = cms.bool(True)
                                           )
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.slimmedJetsJetId = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                           filterParams = pfJetIDSelector.clone(),
                                           #src = cms.InputTag("slimmedJets"),
                                           src = cms.InputTag("updatedPatJetsUpdatedJEC"), 
                                           filter = cms.bool(True)

                                           )


process.demo = cms.EDAnalyzer('Ntupler')
# Select data or MC - this controls which jet corrections are used and whether PU reweighting info is filled                                                          
# Currently the only year+era options are 2017 for MC, and 2017B for data.                                                                                            
process.demo.isMC = cms.bool(True)
#process.demo.isMC = cms.bool(False)
process.demo.year = cms.int32(2017)
process.demo.era = cms.string("C")

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("SlimmedNtuple.TotalEvents.CfiFile_cfi")

process.p = cms.Path(process.totalEvents*
                     process.hltFilter*
                     process.egmGsfElectronIDSequence*
                     process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC*
                     process.patJetCorrFactorsUpdatedJECAK8*
                     process.updatedPatJetsUpdatedJECAK8*
                     process.slimmedJetsAK8JetId*
                     process.slimmedJetsJetId*
                     process.demo)












