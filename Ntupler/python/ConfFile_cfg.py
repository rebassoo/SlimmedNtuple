import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #        'file:/tmp/jjhollar/0C08DB9E-3543-E811-BB4D-0CC47A4D75F4.root'
        #        'file:/tmp/jjhollar/BulkGravToWW_narrow_M-1000_3A7945FC-1704-E811-9B66-008CFA110C5C.root'
        #        'file:/tmp/jjhollar/QCD_HT700to1000_F6F5131E-CCF9-E711-B15F-0242AC130002.root'
        'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/BulkGravToWW_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/B8CED1AB-4D46-E811-B3F9-24BE05BDAE61.root'
    )
)

process.load("SlimmedNtuple.Ntupler.CfiFile_cfi")

process.load("SlimmedNtuple.Ntupler.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.slimmedJetsAK8JetId = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                           filterParams = pfJetIDSelector.clone(),
                                           src = cms.InputTag("slimmedJetsAK8"),
                                           filter = cms.bool(True)
                                           )
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.slimmedJetsJetId = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                           filterParams = pfJetIDSelector.clone(),
                                           src = cms.InputTag("slimmedJets"),
                                           filter = cms.bool(True)
                                           )

process.demo = cms.EDAnalyzer('Ntupler')
# Select data or MC - this controls which jet corrections are used and whether PU reweighting info is filled                                                          
# Currently the only year+era options are 2017 for MC, and 2017B for data.                                                                                            
process.demo.isMC = cms.bool(True)
#process.demo.isMC = cms.bool(False)
process.demo.year = cms.int32(2017)
process.demo.era = cms.string("C")

process.p = cms.Path(#process.hltFilter
                     process.slimmedJetsAK8JetId
                     process.slimmedJetsJetId
                     *process.demo)












