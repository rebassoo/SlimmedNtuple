import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
#process = cms.Process("CTPPSTestProtonReconstruction", eras.ctpps_2016)

#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)

#Details for this here https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('pileupName',
                 'WJetsToLNu_2J_TuneCP5_13TeV_amcatnloFXFX_pythia8',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Name for pileup sample")
options.register('era',
                 'C',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Era")
options.register('interactive',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Determines whether it is an interactive or crab run")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch//store/data/Run2017C/SingleMuon/MINIAOD/31Mar2018-v1/30000/EA6122AF-1137-E811-B552-FA163E28D344.root'
        'file:EA6122AF-1137-E811-B552-FA163E28D344.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/0299C20A-7736-E811-ABC8-008CFAE453D8.root'
        #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/B4005649-E955-E811-BE7B-0CC47A7C353E.root'
        #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/00469E05-E055-E811-807F-008CFAF7174A.root'
        )
)

process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('Configuration.StandardSequences.GeometryDB_cff')

process.load("SlimmedNtuple.Ntupler.CfiFile_cfi")

process.load("SlimmedNtuple.Ntupler.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.load("SlimmedNtuple.Ntupler.METFilter_cfi")


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
MC=False


# If using full proton re-reco (legacy) - local RP reconstruction chain with standard settings                                                                                            
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includeStrips = cms.bool(True)
process.ctppsProtons.doSingleRPReconstruction = cms.bool(True)
process.ctppsProtons.doMultiRPReconstruction = cms.bool(True)

#Global tags from here:https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
#To print out global tag: conddb list 94X_mc2017_realistic_v16
#For v16 global tag uses v23 of Jet energy corrections. Newest is 32.
if MC:
    process.GlobalTag.globaltag ='94X_mc2017_realistic_v17'
else:
    #process.GlobalTag.globaltag ='94X_dataRun2_v11'
    process.GlobalTag.globaltag ='106X_dataRun2_testPPS_v1'


### ADD SOME NEW JET COLLECTIONS                                                                                                              
# New (March 8, 2019) - to recover ak8 CHS jets with 2017 MiniAOD
#################################
  ###  JET TOOLBOX FOR CHS ###
#################################
# AK R=0.8 jets from CHS inputs with basic grooming, W tagging, and top tagging                                                            
from JMEAnalysis.JetToolbox.jetToolbox_cff import *
jetToolbox( process, 'ak8', 'ak8JetSubs', 'noOutput',
#jetToolbox( process, 'ak8', 'jetSequence', 'noOutput',
                PUMethod='CHS',runOnMC=MC,
                addPruning=True, addSoftDrop=False ,           # add basic grooming                                                            
                addTrimming=False, addFiltering=False,
                addSoftDropSubjets=False,
                addNsub=True, maxTau=4,                       # add Nsubjettiness tau1, tau2, tau3, tau4                                      
                Cut='pt > 10.0',
                bTagDiscriminators=['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
                #bTagDiscriminators=['pfCombinedSecondaryVertexV2BJetTags'],
                # added L1FastJet on top of the example config file
                JETCorrPayload = 'AK8PFchs', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
                )

jetToolbox( process, 'ak4', 'ak4JetSubs', 'noOutput',
#jetToolbox( process, 'ak4', 'jetSequence', 'noOutput',
                PUMethod='CHS',runOnMC=MC,
                addPruning=True, addSoftDrop=False ,           # add basic grooming                                                            
                addTrimming=False, addFiltering=False,
                addSoftDropSubjets=False,
                addNsub=True, maxTau=4,                       # add Nsubjettiness tau1, tau2, tau3, tau4                                      
                Cut='pt > 10.0',
                bTagDiscriminators=['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
                #bTagDiscriminators=['pfCombinedSecondaryVertexV2BJetTags'],
                # added L1FastJet on top of the example config file
                JETCorrPayload = 'AK4PFchs', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
                )


if MC:
    #################################
    ###  JER SMEARING AK8###
    #################################
    COLLECTIONTOSMEARAK8="selectedPatJetsAK8PFCHS"
    from RecoMET.METProducers.METSigParams_cfi import *
    process.slimmedAK8JetsSmeared = cms.EDProducer('SmearedPATJetProducer',
        src = cms.InputTag(COLLECTIONTOSMEARAK8),
        enabled = cms.bool(True),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        algo = cms.string('AK8PFchs'),
        algopt = cms.string('AK8PFchs_pt'),

        genJets = cms.InputTag('slimmedGenJets'),
        dRMax = cms.double(0.2),
        dPtMaxFactor = cms.double(3),
        seed = cms.uint32(37428479),
        debug = cms.untracked.bool(False),
    # Systematic variation
    # 0: Nominal
    # -1: -1 sigma (down variation)
    # 1: +1 sigma (up variation)
     variation = cms.int32(0)  # If not specified, default to 0
       )
    #################################
    ###  JER SMEARING AK4###
    #################################
    COLLECTIONTOSMEARAK4="selectedPatJetsAK4PFCHS"
    from RecoMET.METProducers.METSigParams_cfi import *
    process.slimmedAK4JetsSmeared = cms.EDProducer('SmearedPATJetProducer',
        src = cms.InputTag(COLLECTIONTOSMEARAK4),
        enabled = cms.bool(True),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        algo = cms.string('AK4PFchs'),
        algopt = cms.string('AK4PFchs_pt'),

        genJets = cms.InputTag('slimmedGenJets'),
        dRMax = cms.double(0.2),
        dPtMaxFactor = cms.double(3),
        seed = cms.uint32(37424479),
        debug = cms.untracked.bool(False),
    # Systematic variation
    # 0: Nominal
    # -1: -1 sigma (down variation)
    # 1: +1 sigma (up variation)
    variation = cms.int32(0)  # If not specified, default to 0
        )


COLLECTIONFORIDAK8="slimmedAK8JetsSmeared"
COLLECTIONFORIDAK4="slimmedAK4JetsSmeared"
if not MC:
    COLLECTIONFORIDAK8="selectedPatJetsAK8PFCHS"
    COLLECTIONFORIDAK4="selectedPatJetsAK4PFCHS"

    

from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.slimmedJetsAK8JetId = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                           filterParams = pfJetIDSelector.clone(),
                                           #src = cms.InputTag("slimmedJetsAK8"),
                                           #src = cms.InputTag("updatedPatJetsUpdatedJECAK8"),
                                           src = cms.InputTag(COLLECTIONFORIDAK8),
                                           #src = cms.InputTag("selectedPatJetsAK8PFCHS"),
                                           filter = cms.bool(True)
                                           )
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.slimmedJetsJetId = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                           filterParams = pfJetIDSelector.clone(),
                                           #src = cms.InputTag("slimmedJets"),
                                           #src = cms.InputTag("updatedPatJetsUpdatedJEC"), 
                                           src = cms.InputTag(COLLECTIONFORIDAK4),  
                                           #src = cms.InputTag("selectedPatJetsAK4PFCHS"),  
                                           filter = cms.bool(True)
                                           )


from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData= not MC
                           )


process.demo = cms.EDAnalyzer('Ntupler')
# Select data or MC - this controls which jet corrections are used and whether PU reweighting info is filled                           
process.demo.isMC = cms.bool(MC)
process.demo.isSignalMC = cms.bool(False)
process.demo.year = cms.int32(2017)
process.demo.era = cms.string(options.era)
process.demo.isInteractive = cms.bool(bool(options.interactive))
#pileupName="WJetsToLNu_2J_TuneCP5_13TeV_amcatnloFXFX_pythia8"
process.demo.mcName=cms.string("h_pileup_"+options.pileupName)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("SlimmedNtuple.TotalEvents.CfiFile_cfi")

process.load("RecoMET.METFilters.ecalBadCalibFilter_cfi")

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )


process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#Update for MET filter here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
if MC:
    process.p = cms.Path(
        process.totalEvents*
        process.hltFilter*
        process.metFilterMC*
        process.ecalBadCalibReducedMINIAODFilter*
        process.fullPatMetSequence*
        process.egmGsfElectronIDSequence*
        process.slimmedAK8JetsSmeared*
        process.slimmedJetsAK8JetId*
        process.slimmedAK4JetsSmeared*
        process.slimmedJetsJetId*
        #process.ctppsProtonReconstructionOFDB*
        process.demo
        )
    
else:
    process.p = cms.Path(
        #process.hltFilter*
        #process.metFilter*
        #process.totalEvents*
        process.ecalBadCalibReducedMINIAODFilter*
        process.fullPatMetSequence*
        process.egmGsfElectronIDSequence*
        #process.dump*
        process.slimmedJetsAK8JetId*
        process.slimmedJetsJetId*
        #process.ctppsProtonReconstructionOFDB*

        #process.totemRPUVPatternFinder *
        #process.totemRPLocalTrackFitter *
        #process.ctppsDiamondRecHits *
        #process.ctppsDiamondLocalTracks *
        #process.ctppsPixelLocalReconstruction *
        #process.ctppsLocalTrackLiteProducer *
        process.ctppsProtons *
        #process.dump#*
        process.demo
        )


#print process.dumpPython()









