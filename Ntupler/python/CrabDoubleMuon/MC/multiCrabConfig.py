from CRABClient.UserUtilities import config
config = config()

config.section_('General')
#config.General.instance = 'You need to set the CRAB3 instance type: e.g. "private" '
#config.General.serverUrl = 'You need to set the CRAB3 server URL'
#config.General.instance = 'private'
#config.General.serverUrl     = 'crab3-gwms-1.cern.ch'
#config.General.requestName = 'DoubleMuonRun2016C_2ndTry'
config.General.requestName = 'REQ'
config.General.transferLogs = True
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'ConfFile_cfg.py'
#config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles  = ['alignment_collection.out','optics.root','shared_alignment.h','shared_fill_info.h','shared_reconstruction.h','shared_track.h']
config.JobType.inputFiles  = ['alignment_collection.out','xi_as_a_function_of_x_graph_b1.root','xi_as_a_function_of_x_graph_b2.root','MCPileupHighStats.root','MyDataPileupHistogram0to75_MuonPhys.root']

config.section_('Data')
config.Data.inputDataset ='/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/AODSIM'
#config.Data.inputDataset ='DATAS'
#config.Data.inputDataset ='/DoubleMuon/Run2016C-PromptReco-v2/AOD'
#config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5
#config.Data.unitsPerJob = 50000
#config.Data.unitsPerJob = 1000000
#config.Data.unitsPerJob = 40000
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventBased'
config.Data.outLFNDirBase = '/store/user/rebassoo/'
#config.Data.publication = True
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'

config.section_('Site')
config.Site.storageSite = 'T2_US_UCSD'
config.Site.whitelist = ['T2_US_Florida']



if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects_2018_02_21-SingleMuTrigger'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################


#    config.General.requestName = 'DY-10-50'
#    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
#    config.Data.unitsPerJob = 900
#    #config.Data.totalUnits = 8
#    submit(config)

    config.General.requestName = 'DY-50-up-94x-mumu'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM'
#    config.Data.unitsPerJob = 900
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'DY-50-up-94x-Madgraph-mumu'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-RECOSIMstep_94X_mc2017_realistic_v10-v1/AODSIM'
#    config.Data.unitsPerJob = 900
#    config.Data.totalUnits = 8
#    submit(config)


    config.General.requestName = 'DY-50-up-94x-Madgraph-ext1-mumu'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-RECOSIMstep_94X_mc2017_realistic_v10_ext1-v1/AODSIM'
#    config.Data.unitsPerJob = 900
#    config.Data.totalUnits = 8
#    submit(config)


#6105137 events
    config.General.requestName = 'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-mumu'
    config.Data.inputDataset = '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
#    config.Data.unitsPerJob = 900
#    #config.Data.totalUnits = 8
#    submit(config)


#1999000 events
    config.General.requestName = 'WWTo2L2Nu_13TeV-powheg-mumu'
    config.Data.inputDataset = '/WWTo2L2Nu_13TeV-powheg/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
##    config.Data.unitsPerJob = 900
##    #config.Data.totalUnits = 8
#    submit(config)
