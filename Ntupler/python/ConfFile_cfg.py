import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #        'file:/tmp/jjhollar/0C08DB9E-3543-E811-BB4D-0CC47A4D75F4.root'
        #        'file:/tmp/jjhollar/BulkGravToWW_narrow_M-1000_3A7945FC-1704-E811-9B66-008CFA110C5C.root'
        #        'file:/tmp/jjhollar/QCD_HT700to1000_F6F5131E-CCF9-E711-B15F-0242AC130002.root'
        #'root://xrootd.t2.ucsd.edu:2040//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B4005649-E955-E811-BE7B-0CC47A7C353E.root'
        #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B4005649-E955-E811-BE7B-0CC47A7C353E.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2017C/SingleMuon/MINIAOD/17Nov2017-v1/40000/E664484F-7BDB-E711-A6A3-0025904C516C.root'
        #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/0299C20A-7736-E811-ABC8-008CFAE453D8.root'
        'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/CEEB885D-0354-E811-B17D-0CC47A4C8E16.root'
        #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/00469E05-E055-E811-807F-008CFAF7174A.root'
        #'file:step3_fpmc_MiniAOD.root'
        )
)

process.load('Configuration.StandardSequences.GeometryDB_cff')

process.load("SlimmedNtuple.Ntupler.CfiFile_cfi")

process.load("SlimmedNtuple.Ntupler.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.load("SlimmedNtuple.Ntupler.METFilter_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
MC=True

#Global tags from here:https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
#To print out global tag: conddb list 94X_mc2017_realistic_v16
#For v16 global tag uses v23 of Jet energy corrections. Newest is 32.
if MC:
    process.GlobalTag.globaltag ='94X_mc2017_realistic_v17'
else:
    process.GlobalTag.globaltag ='94X_dataRun2_v11'

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

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData= not MC
                           )


process.demo = cms.EDAnalyzer('Ntupler')
# Select data or MC - this controls which jet corrections are used and whether PU reweighting info is filled                                                          
# Currently the only year+era options are 2017 for MC, and 2017B for data.                                                                                            
process.demo.isMC = cms.bool(MC)
#process.demo.isMC = cms.bool(False)
process.demo.year = cms.int32(2017)
process.demo.era = cms.string("C")

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("SlimmedNtuple.TotalEvents.CfiFile_cfi")

process.load("RecoMET.METFilters.ecalBadCalibFilter_cfi")

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )


#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#Update for MET filter here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
if MC:
    process.p = cms.Path(
        process.totalEvents*
        process.metFilterMC*
        process.hltFilter*
        process.ecalBadCalibReducedMINIAODFilter*
        process.fullPatMetSequence*
        #process.dump*
        process.egmGsfElectronIDSequence*
        process.patJetCorrFactorsUpdatedJEC * 
        process.updatedPatJetsUpdatedJEC*
        process.patJetCorrFactorsUpdatedJECAK8*
        process.updatedPatJetsUpdatedJECAK8*
        process.slimmedJetsAK8JetId*
        process.slimmedJetsJetId*
        process.demo
        )
    
else:
    process.p = cms.Path(
        #process.dump*
        process.totalEvents*
        process.metFilter*
        process.hltFilter*
        process.ecalBadCalibReducedMINIAODFilter*
        process.fullPatMetSequence*
        process.egmGsfElectronIDSequence*
        process.patJetCorrFactorsUpdatedJEC * 
        process.updatedPatJetsUpdatedJEC*
        process.patJetCorrFactorsUpdatedJECAK8*
        process.updatedPatJetsUpdatedJECAK8*
        process.slimmedJetsAK8JetId*
        process.slimmedJetsJetId*
        process.demo
        )












