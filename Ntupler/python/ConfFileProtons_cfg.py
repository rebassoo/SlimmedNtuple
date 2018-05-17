import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1001) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/home/users/rebassoo/work/2017_10_23_ProtonTransport/CMSSW_8_0_27/src/SimCTPPS/OpticsParameterisation/test/ctppsSim.root'
        #'file:/home/users/rebassoo/work/2017_10_23_ProtonTransport/CMSSW_8_0_27/src/readHepMC_cff_py_GEN_200k.root'
        #'file:/home/users/rebassoo/work/CMSSW_7_1_26/src/FSQ-RunIISummer15wmLHEGS-00002.root'
        #'file:/home/users/rebassoo/work/2017_12_18_Testing2/CMSSW_8_0_28/src/SimCTPPS/OpticsParameterisation/test/ctppsSim.root'
        'file:/home/users/rebassoo/work/2018_04_27_JansNewRecipeUpdated/CMSSW_8_0_28/src/Validation/CTPPS/test/fast_simu_validation/ctppsSim.root'
        #'file:/home/users/rebassoo/work/2017_10_23_ProtonTransport/CMSSW_8_0_27/src/SubmittingCondorJobs/readHepMC_cff_py_GEN-FPMC-a0w-5e-6.root'
#        'file:/home/users/rebassoo/work/2017_10_23_ProtonTransport/CMSSW_8_0_27/src/readHepMC_cff_py_GEN-FPMCgit_200k_leptons.root'
        #'file:/home/users/rebassoo/work/2017_10_23_ProtonTransport/CMSSW_8_0_27/src/readHepMC_cff_py_GEN_SIM_FPMC.root'
        #'file:/home/users/rebassoo/work/2017_10_23_ProtonTransport/CMSSW_8_0_27/src/SimCTPPS/OpticsParameterisation/test/ctppsSim.root'
        #'file:/home/users/rebassoo/work/2017_10_23_ProtonTransport/CMSSW_8_0_27/src/cepWW_step2.root'
        #'file:LaurentsFileSim.root'
    )#,
   #skipEvents=cms.untracked.uint32(4)
)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.METFilters.badGlobalMuonTaggersAOD_cff')


#ISMC = False
ISMC = True
from Configuration.AlCa.GlobalTag import GlobalTag
if ISMC:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
else:
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v4'  #

#CHANNEL="mue"
process.load("SlimmedNtuple.Ntupler.CfiFile_cfi") 
process.demo.ismc=ISMC
process.demo.ispps=True
process.demo.channel="mue"
process.demo.hepmcCollection="generatorSmeared"

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
useAOD = True
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff'
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
    ]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
if ISMC:
    #process.p = cms.Path(process.egmGsfElectronIDSequence * process.demo)
    process.p = cms.Path(process.demo)
else:
    process.p = cms.Path(process.noBadGlobalMuons * process.egmGsfElectronIDSequence * process.demo)
    #process.p = cms.Path(process.egmGsfElectronIDSequence * process.demo)
#process.p = cms.Path(process.demo)
