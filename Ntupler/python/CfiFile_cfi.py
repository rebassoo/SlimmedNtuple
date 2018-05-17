import FWCore.ParameterSet.Config as cms
import os

TFileService = cms.Service("TFileService", fileName = cms.string("SlimmedNtuple.root") )

print os.getcwd()

demo = cms.EDAnalyzer('Ntupler',
                      ismc=cms.bool(False),
                      ispps=cms.bool(False),
                      channel=cms.string("mue"),
                      alignment = cms.string("alignment_collection.out"),
                      optics = cms.string('optics.root'),
                      hepmcCollection = cms.string('generatorSmeared')
                      )
