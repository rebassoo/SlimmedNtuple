import FWCore.ParameterSet.Config as cms
import copy

# Trigger                                                                                                                                                                     
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
metFilter = copy.deepcopy(hltHighLevel)
metFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","RECO")
metFilter.HLTPaths = ['Flag_goodVertices',
                      'Flag_globalSuperTightHalo2016Filter',
                      'Flag_HBHENoiseFilter',
                      'Flag_HBHENoiseIsoFilter',
                      'Flag_EcalDeadCellTriggerPrimitiveFilter',
                      'Flag_BadPFMuonFilter',
                      'Flag_BadChargedCandidateFilter',
                      'Flag_eeBadScFilter'
                      #This needs to be rerun 'ecalBadCalibReducedMINIAODFilter'
                      ]
metFilter.throw = cms.bool(False)


metFilterMC = copy.deepcopy(hltHighLevel)
metFilterMC.TriggerResultsTag = cms.InputTag("TriggerResults","","PAT")
metFilterMC.HLTPaths = ['Flag_goodVertices',
                        #'Flag_globalSuperTightHalo2016Filter',
                        'Flag_HBHENoiseFilter',
                        'Flag_HBHENoiseIsoFilter',
                        'Flag_EcalDeadCellTriggerPrimitiveFilter',
                        'Flag_BadPFMuonFilter',
                        'Flag_BadChargedCandidateFilter'#,
                        #'Flag_eeBadScFilter'
                      #This needs to be rerun 'ecalBadCalibReducedMINIAODFilter'
                      ]

metFilterMC.throw = cms.bool(False)
