import FWCore.ParameterSet.Config as cms
import copy

# Trigger                                                                                                                                                                     
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
hltFilter = copy.deepcopy(hltHighLevel)
hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
hltFilter.HLTPaths = ['HLT_Ele35_WPTight_Gsf_v*',
                      'HLT_Ele32_WPTight_Gsf_v*',
                      'HLT_IsoMu27_v*',
                      'HLT_IsoMu24_v*']
#hltFilter.HLTPaths = ['HLT_IsoMu27_v*']
#hltFilter.HLTPaths = ['HLT_PFHT1050_v*',
#                      'HLT_AK8PFJet500_v*',
#                      'HLT_AK8PFJet360_TrimMass30_v*',
#                      'HLT_AK8PFJet380_TrimMass30_v*',
#                      'HLT_AK8PFJet400_TrimMass30_v*',
#                      'HLT_AK8PFJet420_TrimMass30_v*',
#                      'HLT_AK8PFHT750_TrimMass50_v*',
#                      'HLT_AK8PFHT800_TrimMass50_v*',
#                      'HLT_AK8PFHT850_TrimMass50_v*',
#                      'HLT_AK8PFHT900_TrimMass50_v*']

hltFilter.throw = cms.bool(False)
