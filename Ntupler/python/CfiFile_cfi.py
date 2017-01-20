import FWCore.ParameterSet.Config as cms
import os

TFileService = cms.Service("TFileService", fileName = cms.string("SlimmedNtuple.root") )

print os.getcwd()

demo = cms.EDAnalyzer('Ntupler',
                      ismc=cms.bool(False),
                      channel=cms.string("mue"),
                      particleFile = cms.string("alignment_collection.out"),
                      particleFile2 = cms.string('optics.root')
                      )
