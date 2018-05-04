import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/home/users/rebassoo/EED9C2D8-BF8F-E611-B8CE-A0369F7FE9FC.root'
        #'root://cmsxrootd.fnal.gov//store/data/Run2016F/DoubleMuon/AOD/23Sep2016-v1/50000/9E0F8266-0590-E611-A14E-0242AC130002.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2017B/DoubleMuon/AOD/17Nov2017-v1/30000/001137DB-26D6-E711-8E99-02163E014465.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/94X_mc2017_realistic_v10-v1/00000/003FC78E-80EC-E711-8B6A-008CFAC93D9C.root'
        'file:/hadoop/cms/store/user/rebassoo/TestFiles/DYJets2017-003FC78E-80EC-E711-8B6A-008CFAC93D9C.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2017C/DoubleMuon/AOD/17Nov2017-v1/30000/507F651E-CCD8-E711-8022-00259020084C.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2017F/DoubleMuon/AOD/17Nov2017-v1/70000/BAF4959F-7BDE-E711-A6CD-02163E0138BD.root'
        )#,
         #                   skipEvents=cms.untracked.uint32(169)
)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.METFilters.badGlobalMuonTaggersAOD_cff')


ISMC = False
#ISMC = True
from Configuration.AlCa.GlobalTag import GlobalTag
if ISMC:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
else:
    process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'  #
    #process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.load("SlimmedNtuple.Ntupler.CfiFile_cfi") 
process.demo.ismc=ISMC
process.demo.ispps=True
process.demo.channel="mumu"

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
useAOD = True
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

process.load("PhysicsTools.PatAlgos.patSequences_cff")
#from Configuration.EventContent.EventContent_cff import *

# define which IDs we want to produce
my_id_modules = [
#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff'
#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
    ]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.tightPatJetsPFlow = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                         filterParams = pfJetIDSelector.clone(quality=cms.string("TIGHT")),
                                         src = cms.InputTag("selectedPatJets")
                                         #src = cms.InputTag("cleanPatJets")
                                         )

if not ISMC:
    process.patJetCorrFactors.levels = cms.vstring('L1FastJet','L2Relative','L3Absolute','L2L3Residual')


process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
if ISMC:
    #process.p = cms.Path(process.egmGsfElectronIDSequence * process.patDefaultSequence* process.tightPatJetsPFlow*process.dump*process.demo)
    #process.p = cms.Path(process.egmGsfElectronIDSequence * process.patDefaultSequence* process.tightPatJetsPFlow*process.demo)
    process.p = cms.Path(process.egmGsfElectronIDSequence * process.makePatJets* process.selectedPatJets*process.tightPatJetsPFlow*process.demo)
    #process.p = cms.Path(process.egmGsfElectronIDSequence * process.demo)
else:
    #process.p = cms.Path(process.egmGsfElectronIDSequence * process.makePatJets* process.selectedPatJets*process.tightPatJetsPFlow*process.demo)
    process.p = cms.Path(process.egmGsfElectronIDSequence * process.patJetCorrections+process.patJetCharge*process.patJets*process.selectedPatJets*process.tightPatJetsPFlow*process.demo)
    #process.p = cms.Path(process.noBadGlobalMuons * process.egmGsfElectronIDSequence * process.patDefaultSequence*process.demo)
    #process.p = cms.Path(process.egmGsfElectronIDSequence * process.demo)
#process.p = cms.Path(process.demo)
