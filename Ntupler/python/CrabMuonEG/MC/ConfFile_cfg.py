import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1001) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        ##'file:/tmp/rebassoo/miniAOD_PAT_1.root'
        ##'file:/tmp/rebassoo/B2846F07-B046-E611-8CC0-0CC47A13CDA0.root'
        ##'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90001/B2846F07-B046-E611-8CC0-0CC47A13CDA0.root'
        ##'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90001/E49F4520-B046-E611-B7DF-003048F5B2F0.root'
        ##'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90000/1C4F8054-D845-E611-8251-FA163E149DE8.root'
        ##'root://cmsxrootd.fnal.gov//store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/20000/004938DD-26FC-E511-A80F-02163E017620.root'
        ##'root://cms-xrd-global.cern.ch//store/data/Run2016C/MuonEG/AOD/23Sep2016-v1/70000/0009D8A4-DC87-E611-9090-0CC47A044888.root'
        ##/SingleElectron/Run2016B-01Jul2016-v1/AOD
        ##'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90002/8A922D2A-5C47-E611-B3E7-001A648F1BFA.root'
        ##'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90000/A656B0C9-BF45-E611-9921-002590200A28.root'
        ##'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90000/CCB0A6F2-B145-E611-A8AB-B083FEC76567.root'
        ##'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/MINIAOD/01Jul2016-v1/90000/FEF32273-1647-E611-9AA3-848F69FD2D6F.root'
        ##'file:/hadoop/cms/store/user/rebassoo/58D3DB31-B9B1-E611-8F7A-001EC94BA173.root'
        ##'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16DR80Premix/WWTo2L2Nu_13TeV-powheg/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/04F11512-AEB0-E611-8949-0242AC130008.root'
        ##'root://cms-xrd-global.cern.ch//store/data/Run2016D/MuonEG/AOD/23Sep2016-v1/100000/2649EE95-C388-E611-AAD1-3417EBE64402.root'
        #'root://xrootd.pic.es:1094/pnfs/pic.es/data/cms/disk/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/58D3DB31-B9B1-E611-8F7A-001EC94BA173.root'
        #'file:pickevents_2016B.root',
        #'file:pickevents_2016C.root',
        #'file:pickevents_2016G.root'
        #'file:readHepMC_cff_py_GEN_SIM_RECOBEFMIX_DIGI_RECO.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/MuonEG-RunG-D43CFA2B-C192-E611-AC1E-0090FAA57BA0.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/WWTo2L2Nu_13TeV-powheg-2E768F23-D8B0-E611-B3D3-0242AC130004.root'
        #'file:SUS-RunIISummer16DR80Premix-00036.root',
        #'file:readHepMC_cff_py_GEN_SIM_RECOBEFMIX_DIGI_RECO.root',
        #'file:/home/users/rebassoo/work/CMSSW_8_0_21/src/SUS-RunIISummer16DR80Premix-00036.root',
        #'file:FSQ-RunIISummer16DR80Premix-00036.root'
        #'file:/home/users/rebassoo/work/CMSSW_8_0_21/src/2017-08-10-SMExclusiveWW-25k/FSQ-RunIISummer16DR80Premix-00036.root'
        #'file:/home/users/rebassoo/work/CMSSW_8_0_21/src/2017-08-10-SMExclusiveWW-2nd-25k/FSQ-RunIISummer16DR80Premix-00036.root'
        #'file:/hadoop/cms/store/user/rebassoo/GammaGammaMuMuFiles/GammaGammaMuMu-Elastic-RunIIWinter15R-Summer16Moriond17-0002.root'
        #'root://se01.indiacms.res.in//store/mc/RunIISummer16DR80Premix/WWTo2L2Nu_13TeV-powheg/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/6C709675-36B2-E611-8E68-FA163E8DA59B.root',
        #'root://se01.indiacms.res.in//store/mc/RunIISummer16DR80Premix/WWTo2L2Nu_13TeV-powheg/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/746CEFB5-35B1-E611-A5B3-001E674FAF23.root',
        #'root://se01.indiacms.res.in//store/mc/RunIISummer16DR80Premix/WWTo2L2Nu_13TeV-powheg/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/7496F3A6-A8B0-E611-BA28-001E67456F68.root',
        #'root://se01.indiacms.res.in//store/mc/RunIISummer16DR80Premix/WWTo2L2Nu_13TeV-powheg/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/786ED34A-DAB0-E611-9FA4-001E674FB063.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/DoubleEG-RunB-A6C55E4E-E997-E611-A869-008CFA1111F4.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/MuonEG-PromptReco-v3_0C2205B1-A29F-E611-8E3F-02163E014472.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/MuonEG-RunD_2649EE95-C388-E611-AAD1-3417EBE64402.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/MuonEG-Run2016H-v3-087BD90F-6B9F-E611-912F-02163E01448A.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/MuonEG-RunD-CE007154-EB88-E611-A13B-0025905C54DA.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/WWTo2L2Nu_13TeV-powheg-herwigpp_02B67C54-09B2-E611-866F-6C3BE5B51168.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-129B1D42-5AB0-E611-B1B0-A4BADB22B643.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/MuonEG-AOD-Run2016G-Aug17-081723E8-AD93-E711-B74E-24BE05C68671.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/MuonEG-AOD-Run2016G-Sep23-EED9C2D8-BF8F-E611-B8CE-A0369F7FE9FC.root'
        #'file:/hadoop/cms/store/user/rebassoo/TestFiles/MuonEG-AOD-Run2016-Sep16-485BAA4D-FA8B-E611-BAB7-0025905C54C6.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2016G/DoubleMuon/AOD/23Sep2016-v1/80001/EED9C2D8-BF8F-E611-B8CE-A0369F7FE9FC.root'
        #'file:/home/users/rebassoo/EED9C2D8-BF8F-E611-B8CE-A0369F7FE9FC.root'
        #'root://cmsxrootd.fnal.gov//store/data/Run2016F/DoubleMuon/AOD/23Sep2016-v1/50000/9E0F8266-0590-E611-A14E-0242AC130002.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2017B/DoubleMuon/AOD/17Nov2017-v1/30000/001137DB-26D6-E711-8E99-02163E014465.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/94X_mc2017_realistic_v10-v1/00000/003FC78E-80EC-E711-8B6A-008CFAC93D9C.root'
        'file:/hadoop/cms/store/user/rebassoo/TestFiles/DYJets2017-003FC78E-80EC-E711-8B6A-008CFAC93D9C.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/94X_mc2017_realistic_v10-v1/00000/003FC78E-80EC-E711-8B6A-008CFAC93D9C.root',
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/94X_mc2017_realistic_v10-v1/00000/020D88F7-E5F1-E711-8BE1-003048F35234.root',
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/94X_mc2017_realistic_v10-v1/00000/022AD7D0-87EC-E711-8C10-008CFAE4504C.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/RECOSIMstep_94X_mc2017_realistic_v10_ext1-v1/00000/06E07522-61F0-E711-9843-02163E0119C8.root'
        )#,
         #                   skipEvents=cms.untracked.uint32(169)
)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('RecoMET.METFilters.badGlobalMuonTaggersAOD_cff')


#ISMC = False
ISMC = True
from Configuration.AlCa.GlobalTag import GlobalTag
if ISMC:
    #process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
    process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v10')
else:
    process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'  #
    #process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
    #process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

#CHANNEL="mue"
process.load("SlimmedNtuple.Ntupler.CfiFile_cfi") 
process.demo.ismc=ISMC
process.demo.ispps=False
process.demo.channel="mue"

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

#process.patJets.addGenJetMatch = cms.bool(False)
#process.patJets.addGenPartonMatch = cms.bool(False)
#process.patJets.getJetMCFlavour = cms.bool(False)
#process.patJets.embedGenJetMatch = cms.bool(False)
#process.patJets.embedGenPartonMatch = cms.bool(False)

#process.cleanPatJets.checkOverlaps.electrons.requireNoOverlaps = cms.bool(True)
#process.cleanPatJets.checkOverlaps.muons.requireNoOverlaps = cms.bool(True)
#process.cleanPatJets.checkOverlaps.photons.requireNoOverlaps = cms.bool(True)
#process.cleanPatJets.checkOverlaps.taus.requireNoOverlaps = cms.bool(True)

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
