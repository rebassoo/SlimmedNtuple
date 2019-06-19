from CRABClient.UserUtilities import config
import glob
#import os
config = config()

config.section_('General')
config.General.requestName = 'REQ'
config.General.transferLogs = True
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'ConfFile_cfg.py'
#fileList = [os.path.basename(x) for x in glob.glob('2017-JEC-JER/*.txt')]
config.JobType.inputFiles  = glob.glob("2017-JEC-JER/*.txt")

config.section_('Data')
#config.Data.inputDataset ='/DoubleMuon/Run2016B-01Jul2016-v1/AOD'
config.Data.inputDataset ='DATAS'
#config.Data.inputDataset ='/DoubleMuon/Run2016C-PromptReco-v2/AOD'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 50
#config.Data.unitsPerJob = 100
#config.Data.unitsPerJob = 40000
#config.Data.splitting = 'EventBased'
config.Data.outLFNDirBase = '/store/user/rebassoo/'
config.Data.lumiMask = 'combined_RPIN_CMS.json'
#config.Data.lumiMask = 'notFinishedLumis.json'
#config.Data.lumiMask = 'notFinishedLumis-RunC.json'
#config.Data.ignoreLocality=True

config.section_('Site')
config.Site.storageSite = 'T2_US_UCSD'
#config.Site.blacklist = ['T2_ES_IFCA']
#config.Site.blacklist = ['T2_US_Caltech','T2_US_UCSD']
#config.Site.whitelist = ['T2_US_Florida','T2_US_MIT','T2_US_UCSD']
#config.Site.whitelist = ['T2_US_UCSD']
#config.Site.whitelist = ['T2_IT_*','T2_FR_*','T2_DE_*','T2_CH_*','T2_UK_*']



if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects_2019_06_19_Fix'

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

    config.General.requestName = 'Run2017B-Dilepton'
    config.Data.inputDataset = '/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD'
    config.JobType.pyCfgParams=['era=B','interactive=0']
    #config.Data.unitsPerJob = 20
    submit(config)
    
    config.General.requestName = 'Run2017C-Dilepton'
    config.Data.inputDataset = '/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD'
    config.JobType.pyCfgParams=['era=C','interactive=0']
    #config.Data.unitsPerJob = 50
    submit(config)
    
    config.General.requestName = 'Run2017D-Dilepton'
    config.Data.inputDataset = '/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD'
    #config.Data.unitsPerJob = 50
    config.JobType.pyCfgParams=['era=D','interactive=0']
    submit(config)
    

    config.General.requestName = 'Run2017E-Dilepton'
    config.Data.inputDataset = '/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD'
    #config.Data.unitsPerJob = 50
    config.JobType.pyCfgParams=['era=E','interactive=0']
    submit(config)


    config.General.requestName = 'Run2017F-Dilepton'
    config.Data.inputDataset = '/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD'
    config.JobType.pyCfgParams=['era=F','interactive=0']
    #config.Data.unitsPerJob = 50
    submit(config)


