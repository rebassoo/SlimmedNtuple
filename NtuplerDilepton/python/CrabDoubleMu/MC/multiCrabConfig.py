from CRABClient.UserUtilities import config
import glob
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
config.JobType.inputFiles  = glob.glob("2017-JEC-JER/*.txt")+["PUHistos_data.root","PUHistos_mc.root"]

config.section_('Data')
#config.Data.inputDataset ='/DoubleMuon/Run2016B-01Jul2016-v1/AOD'
config.Data.inputDataset ='DATAS'
#config.Data.inputDataset ='/DoubleMuon/Run2016C-PromptReco-v2/AOD'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 50
#config.Data.unitsPerJob = 100
#config.Data.unitsPerJob = 40000
#config.Data.splitting = 'EventBased'
config.Data.outLFNDirBase = '/store/user/rebassoo/'
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#config.Data.lumiMask = 'combined_RPIN_CMS.json'
#config.Data.ignoreLocality=True

config.section_('Site')
config.Site.storageSite = 'T2_US_UCSD'
#config.Site.blacklist = ['T2_ES_IFCA']
#config.Site.blacklist = ['T2_US_Purdue']
#config.Site.whitelist = ['T2_US_Florida','T2_US_MIT','T2_US_UCSD']
#config.Site.whitelist = ['T2_US_*']
#config.Site.whitelist = ['T2_IT_*','T2_FR_*','T2_DE_*','T2_CH_*','T2_UK_*']



if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects_2019_06_19-MC-Fix'

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

    config.General.requestName = 'DY1JetsToLL_M-50_LHEZpT_50-150'
    config.Data.inputDataset = '/DY1JetsToLL_M-50_LHEZpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams=['interactive=0','pileupName=DY1JetsToLL_M_50_LHEZpT_50_150_TuneCP5_13TeV_amcnloFXFX_pythia8']
    submit(config)


    config.General.requestName = 'DY1JetsToLL_M-50_LHEZpT_150-250'
    config.Data.inputDataset = '/DY1JetsToLL_M-50_LHEZpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams=['interactive=0','pileupName=DY1JetsToLL_M_50_LHEZpT_150_250_TuneCP5_13TeV_amcnloFXFX_pythia8']
    submit(config)

    config.General.requestName = 'DY1JetsToLL_M-50_LHEZpT_250-400'
    config.Data.inputDataset = '/DY1JetsToLL_M-50_LHEZpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams=['interactive=0','pileupName=DY1JetsToLL_M_50_LHEZpT_250_400_TuneCP5_13TeV_amcnloFXFX_pythia8']
    submit(config)

    config.General.requestName = 'DY1JetsToLL_M-50_LHEZpT_400-inf'
    config.Data.inputDataset = '/DY1JetsToLL_M-50_LHEZpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams=['interactive=0','pileupName=DY1JetsToLL_M_50_LHEZpT_400_inf_TuneCP5_13TeV_amcnloFXFX_pythia8']
    submit(config)

    config.General.requestName = 'DY2JetsToLL_M-50_LHEZpT_50-150'
    config.Data.inputDataset = '/DY2JetsToLL_M-50_LHEZpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams=['interactive=0','pileupName=DY2JetsToLL_M_50_LHEZpT_50_150_TuneCP5_13TeV_amcnloFXFX_pythia8']
    submit(config)


    config.General.requestName = 'DY2JetsToLL_M-50_LHEZpT_150-250'
    config.Data.inputDataset = '/DY2JetsToLL_M-50_LHEZpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams=['interactive=0','pileupName=WZ_TuneCP5_13TeV_pythia8']
    submit(config)

    config.General.requestName = 'DY2JetsToLL_M-50_LHEZpT_250-400'
    config.Data.inputDataset = '/DY2JetsToLL_M-50_LHEZpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams=['interactive=0','pileupName=DY2JetsToLL_M_50_LHEZpT_250_400_TuneCP5_13TeV_amcnloFXFX_pythia8']
    submit(config)

    config.General.requestName = 'DY2JetsToLL_M-50_LHEZpT_400-inf'
    config.Data.inputDataset = '/DY2JetsToLL_M-50_LHEZpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams=['interactive=0','pileupName=DY2JetsToLL_M_50_LHEZpT_400_inf_TuneCP5_13TeV_amcnloFXFX_pythia8']
    submit(config)



