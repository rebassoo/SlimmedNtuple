import FWCore.ParameterSet.Config as cms
import os

TFileService = cms.Service("TFileService", fileName = cms.string("SlimmedNtuple.root") )

print os.getcwd()



demo = cms.EDAnalyzer('Ntupler',
                      #particleFile = cms.FileInPath('/home/users/rebassoo/work/2016_11_14_FinnNtupler/CMSSW_8_1_0_pre8/src/SlimmedNtuple/Ntupler/plugins/alignment_collection.out'),
                      #particleFile2 = cms.FileInPath('/home/users/rebassoo/work/2016_11_14_FinnNtupler/CMSSW_8_1_0_pre8/src/SlimmedNtuple/Ntupler/plugins/optics.root')#,
                      particleFile = cms.string("alignment_collection.out"),
                      particleFile2 = cms.string('optics.root')#,
                      #eleIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1")
                      #eleIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")
                      )
