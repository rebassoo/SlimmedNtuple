#include "SlimmedNtuple/Ntupler/plugins/Ntupler.h"

//
// constants, enums and typedefs
//
//
// static data member definitions
//

//
// constructors and destructor
//
Ntupler::Ntupler(const edm::ParameterSet& iConfig):
  jet_token_(consumes<edm::View<pat::Jet>>(edm::InputTag(("slimmedJetsJetId")))),
  jetAK8_token_(consumes<edm::View<pat::Jet>>(edm::InputTag(("slimmedJetsAK8JetId")))),
  muon_token_(consumes<edm::View<pat::Muon>>(edm::InputTag("slimmedMuons"))),
  pps_token_(consumes<std::vector<CTPPSLocalTrackLite>>(edm::InputTag(("ctppsLocalTrackLiteProducer")))),
  recoProtonsSingleRPToken_   ( consumes<std::vector<reco::ForwardProton> >      (edm::InputTag("ctppsProtons","singleRP" ) ) ),
  recoProtonsMultiRPToken_ ( consumes<std::vector<reco::ForwardProton> > ( edm::InputTag("ctppsProtons","multiRP" ) ) ),
  vertex_token_(consumes<std::vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"))),
  rho_token_(consumes<double>(edm::InputTag(("fixedGridRhoAll")))),
  hlt_token_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))),
  pu_token_(consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"))),
  gen_part_token_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  gen_jet_token_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJetsAK8"))),
  met_token_(consumes<pat::METCollection>(edm::InputTag("slimmedMETsPuppi"))),
  //met_token_(consumes<pat::METCollection>(edm::InputTag("slimmedMETs"))),
  ecalBadCalibFilterUpdate_token(consumes< bool >(edm::InputTag("ecalBadCalibReducedMINIAODFilter"))),
  electron_token_(consumes<edm::View<pat::Electron>>(edm::InputTag("slimmedElectrons"))),
  //electron_token_(consumes<edm::View<reco::GsfElectron>>(edm::InputTag("slimmedElectrons"))),
  pfcand_token_(consumes<edm::View<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"))),
  mcweight_token_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
  
  usesResource("TFileService");
  isMC = iConfig.getParameter<bool>("isMC");
  isSignalMC = iConfig.getParameter<bool>("isSignalMC");
  year = iConfig.getParameter<int>("year");
  era = iConfig.getParameter<std::string>("era");
  mcName = iConfig.getParameter<std::string>("mcName");
  isInteractive = iConfig.getParameter<bool>("isInteractive");

  cout<<"Year: "<<year<<endl;
  cout<<"era: "<<era<<endl;
  cout<<"mcName: "<<mcName<<endl;
  cout<<"isMC: "<<isMC<<endl;

  //eleIdMapToken_=consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"));
  //eleIdMapToken_veto_=consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"));

  if(isMC == true)
    {
      //reco_protons_token_=consumes<std::vector<reco::ProtonTrack>>(edm::InputTag("ctppsProtonReconstructionOFDB"));
      //LumiWeights = new edm::LumiReWeighting("PUHistos_mc.root", "PUHistos_data.root", "pileup", "pileup");
      if(year==2017){
	LumiWeights = new edm::LumiReWeighting("PUHistos_mc.root", "PUHistos_data.root", mcName, "pileup");}
      if(year==2018){
	LumiWeights = new edm::LumiReWeighting("PUHistos_mc_2018.root", "PUHistos_data_2018.root","pileup", "pileup");}

    }
  
  std::vector<std::string> jecAK8PayloadNames_;
  std::vector<std::string> jecAK8PayloadNames_withL1_;
  std::string prefix;
  if (isInteractive==true){
    if(year==2017){
      prefix="2017-JEC-JER/";}
    if(year==2018){
      prefix="2018-JEC-JER/";}
  }
  else{
    prefix="";
  }

  if(isMC==false && year==2017 && era == "B")
    {
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017B_V32_DATA_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017B_V32_DATA_L3Absolute_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK8PFchs.txt");
    }
  if(isMC==false && year==2017 && era == "C")
    {
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017C_V32_DATA_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017C_V32_DATA_L3Absolute_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017C_V32_DATA_L2L3Residual_AK8PFchs.txt");

    }
  if(isMC==false && year==2017 && era == "D")
    {
       jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017DE_V32_DATA_L2Relative_AK8PFchs.txt");
       jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017DE_V32_DATA_L3Absolute_AK8PFchs.txt");
       jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017DE_V32_DATA_L2L3Residual_AK8PFchs.txt");
    }
  if(isMC==false && year==2017 && era == "E")
    {
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017DE_V32_DATA_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017DE_V32_DATA_L3Absolute_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017DE_V32_DATA_L2L3Residual_AK8PFchs.txt");
    }
  if(isMC==false && year==2017 && era == "F")
    {
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017F_V32_DATA_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017F_V32_DATA_L3Absolute_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017F_V32_DATA_L2L3Residual_AK8PFchs.txt");
    }
  
  if(isMC==true && year==2017)
    {
       jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017_V32_MC_L2Relative_AK8PFchs.txt");
       jecAK8PayloadNames_.push_back(prefix+"Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFchs.txt");
       
       //jecAK8PayloadNames_withL1_.push_back(prefix+"Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFchs.txt");
       //jecAK8PayloadNames_withL1_.push_back(prefix+"Fall17_17Nov2017_V32_MC_L2Relative_AK8PFchs.txt");
       //jecAK8PayloadNames_withL1_.push_back(prefix+"Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFchs.txt");
    }

  if(isMC==true && year==2018)
    {
      jecAK8PayloadNames_.push_back("Autumn18_V19_MC_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_V19_MC_L3Absolute_AK8PFchs.txt");
    }
  if(isMC==false && year==2018 && era == "A")
    {
      jecAK8PayloadNames_.push_back("Autumn18_RunA_V19_DATA_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_RunA_V19_DATA_L3Absolute_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_RunA_V19_DATA_L2L3Residual_AK8PFchs.txt");
    }
  if(isMC==false && year==2018 && era == "B")
    {
      jecAK8PayloadNames_.push_back("Autumn18_RunB_V19_DATA_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_RunB_V19_DATA_L3Absolute_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_RunB_V19_DATA_L2L3Residual_AK8PFchs.txt");
    }
  if(isMC==false && year==2018 && era == "C")
    {
      jecAK8PayloadNames_.push_back("Autumn18_RunC_V19_DATA_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_RunC_V19_DATA_L3Absolute_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_RunC_V19_DATA_L2L3Residual_AK8PFchs.txt");
    }
  if(isMC==false && year==2018 && era == "D")
    {
      jecAK8PayloadNames_.push_back("Autumn18_RunD_V19_DATA_L2Relative_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_RunD_V19_DATA_L3Absolute_AK8PFchs.txt");
      jecAK8PayloadNames_.push_back("Autumn18_RunD_V19_DATA_L2L3Residual_AK8PFchs.txt");
    }
  /*
   std::vector<JetCorrectorParameters> vPar_withL1;
   for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_withL1_.begin(),
	   payloadEnd = jecAK8PayloadNames_withL1_.end(),
	   ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar_withL1.push_back(pars);
   }
  */
   
   std::vector<JetCorrectorParameters> vPar;
   for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(),
	   payloadEnd = jecAK8PayloadNames_.end(),
	   ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
   }


   // Make the FactorizedJetCorrector                                                                                                                                                      
   jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
   //jecAK8_withL1_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar_withL1) );
   

   // Get JER smearing                                                                                                                                               
   if(isMC==true && year==2017) // Note - here we're using Summer16 for 2017 MC, until the 2017 version is ready
     {                                     
       jerAK8chsName_res_ = prefix+"Fall17_V3_MC_PtResolution_AK8PFchs.txt";
       jerAK8chsName_sf_ = prefix+"Fall17_V3_MC_SF_AK8PFchs.txt";

       jerAK4chsName_res_ = prefix+"Fall17_V3_MC_PtResolution_AK4PFchs.txt";
       jerAK4chsName_sf_ = prefix+"Fall17_V3_MC_SF_AK8PFchs.txt";


     }


   if(isMC==true && year==2018)
     {
       jerAK8chsName_res_ = prefix+"Autumn18_V7b_MC_PtResolution_AK8PFchs.txt";
       jerAK8chsName_sf_ = prefix+"Autumn18_V7b_MC_SF_AK8PFchs.txt";
       jerAK4chsName_res_ = prefix+"Autumn18_V7b_MC_PtResolution_AK4PFchs.txt";
       jerAK4chsName_sf_ = prefix+"Autumn18_V7b_MC_SF_AK8PFchs.txt";

     }



}


Ntupler::~Ntupler()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< bool > passecalBadCalibFilterUpdate ;
  iEvent.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate);
  bool    _passecalBadCalibFilterUpdate =  (*passecalBadCalibFilterUpdate );
  //*ecalBadCalFilter_=_passecalBadCalibFilterUpdate;
  if(!_passecalBadCalibFilterUpdate){return;}
  using namespace edm;
  using namespace std;


  // HLT
  //cout<<"I get here"<<endl;
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByToken(hlt_token_,hltResults);
  const edm::TriggerNames & trigNames = iEvent.triggerNames(*hltResults);
  std::string TriggerPrefix_mu = "HLT_IsoMu27_";
  std::string TriggerPrefix_e = "HLT_Ele35_WPTight_Gsf_v";
  bool passTrigger_mu=false;
  bool passTrigger_e=false;
  for(unsigned int i=0; i<trigNames.size();i++)
    {
      //cout<<"I get here 1"<<endl;
      int prescale_value=hltPrescaleProvider_.prescaleValue(iEvent, iSetup,trigNames.triggerName(i));
      //       int result = hltResults->accept(i);
      std::string TriggerName = trigNames.triggerName(i);
      //cout<<"TriggerName: "<<TriggerName<<endl;
      if(TriggerName.find(TriggerPrefix_mu) != std::string::npos){
	if((hltResults->accept(i)>0)&&(prescale_value==1)){
	  passTrigger_mu=true;
	}
      }
      
      if(TriggerName.find(TriggerPrefix_e) != std::string::npos){
	if((hltResults->accept(i)>0)&&(prescale_value==1)){
	  passTrigger_e=true;
	 }
      }
    }
  // Run and vertex multiplicity info
  *run_ = iEvent.id().run();
  *ev_ = iEvent.id().event();
  *lumiblock_ = iEvent.luminosityBlock();

  //get beam-crossing angle and LHC conditions
  if(isMC==false){
    edm::ESHandle<LHCInfo> hLHCInfo;
    std::string lhcInfoLabel("");
    iSetup.get<LHCInfoRcd>().get(lhcInfoLabel, hLHCInfo);
    if(hLHCInfo.isValid()){
      *crossingAngle_=hLHCInfo->crossingAngle();
      *betaStar_=hLHCInfo->betaStar();
      *instLumi_=hLHCInfo->instLumi();
    }
  }
  edm::Handle< std::vector<reco::Vertex> > vertices_;
  iEvent.getByToken(vertex_token_, vertices_);
  reco::VertexRef vtx(vertices_, 0);
   //Get Muons
  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByToken(muon_token_,muonHandle);
  int numMuLoose=0;
  for (const pat::Muon &MuonIt : *muonHandle) { 
    bool tightId=MuonIt.isTightMuon(*vtx);
    double iso = (MuonIt.pfIsolationR04().sumChargedHadronPt + max(0., MuonIt.pfIsolationR04().sumNeutralHadronEt + MuonIt.pfIsolationR04().sumPhotonEt - 0.5*MuonIt.pfIsolationR04().sumPUPt))/MuonIt.pt();
    if(tightId&&MuonIt.pt()>35&&fabs(MuonIt.eta())<2.4&&iso<0.1){
      (*muon_px_).push_back(MuonIt.px());
      (*muon_py_).push_back(MuonIt.py());
      (*muon_pz_).push_back(MuonIt.pz());
      (*muon_e_).push_back(MuonIt.energy());
      (*muon_charge_).push_back(MuonIt.charge());
      (*muon_pt_).push_back(MuonIt.pt());
      (*muon_eta_).push_back(MuonIt.eta());
      (*muon_phi_).push_back(MuonIt.phi());
      (*muon_iso_).push_back(iso);
      (*muon_dxy_).push_back(fabs(MuonIt.muonBestTrack()->dxy(vtx->position())));
      (*muon_dz_).push_back(fabs(MuonIt.muonBestTrack()->dz(vtx->position())));
    }
    bool looseId=MuonIt.isLooseMuon();
    if(looseId&&MuonIt.pt()>20&&fabs(MuonIt.eta())<2.4&&iso<0.1){numMuLoose++;}
    
  }//end of looping over muons

   //Get Electrons
   /*
   //Tight id
   edm::Handle<edm::ValueMap<bool> > ele_id_decisions;
   iEvent.getByToken(eleIdMapToken_ ,ele_id_decisions);
   //Veto id
   edm::Handle<edm::ValueMap<bool> > ele_id_decisions_veto;
   iEvent.getByToken(eleIdMapToken_veto_ ,ele_id_decisions_veto);
   */
   //edm::Handle<edm::View<reco::GsfElectron> > electronHandle;
   edm::Handle<edm::View<pat::Electron> > electronHandle;
   iEvent.getByToken(electron_token_,electronHandle);
   unsigned int n = electronHandle->size();
   //cout<<"I get before electrons: "<<endl;
   int num_ele_veto=0;
   for(unsigned int i = 0; i < n; ++i) {
     const auto el = electronHandle->ptrAt(i);
     //reco::GsfElectronRef ele(electronHandle, i);
     //reco::GsfTrackRef theTrack = el->gsfTrack();
     //dz_.push_back( theTrack->dz( firstGoodVertex->position() ) );
     //dxy_.push_back( theTrack->dxy( firstGoodVertex->position() ) );
     //bool passConvVeto = !ConversionTools::hasMatchedConversion(*ele,conversions_h,theBeamSpot->position());

     //bool isPassEleId = (*ele_id_decisions)[el];
     bool isPassEleId = el->electronID("cutBasedElectronID-Fall17-94X-V2-tight");
     //bool isPassEleId = el->isElectron("cutBasedElectronID-Fall17-94X-V2-tight");
     //if(el->pt()>30){ 
       //       cout<<"electron pt is: "<<el->pt()<<endl;
     //}
     
     if(isPassEleId&&el->pt()>40&&fabs(el->superCluster()->eta())<2.4){ 
       (*electron_pt_).push_back(el->pt());
       (*electron_dxy_).push_back(fabs(el->gsfTrack()->dxy(vtx->position())));
       (*electron_dz_).push_back(fabs(el->gsfTrack()->dz(vtx->position())));
       (*electron_eta_).push_back(el->superCluster()->eta());
       (*electron_phi_).push_back(el->superCluster()->phi());
       (*electron_px_).push_back(el->px());
       (*electron_py_).push_back(el->py());
       (*electron_pz_).push_back(el->pz());
       (*electron_e_).push_back(el->energy());
       (*electron_charge_).push_back(el->charge());
     }

     //bool isPassEleId_veto = (*ele_id_decisions_veto)[el];
     bool isPassEleId_veto = el->electronID("cutBasedElectronID-Fall17-94X-V2-veto");
     if(isPassEleId_veto&&el->pt()>35&&fabs(el->superCluster()->eta())<2.4){ 
       num_ele_veto++;
     }

   }
   //Get Gen Jets
   if(isMC == true){
     edm::Handle<reco::GenJetCollection> genJet;
     iEvent.getByLabel("slimmedGenJetsAK8",genJet);
     
     //TLorentzVector genjet1, genjet2, genjj;
     for (reco::GenJetCollection::const_iterator genJetIter=genJet->begin(); genJetIter != genJet->end(); genJetIter++ ) {
       (*gen_jet_pt_).push_back(genJetIter->pt());
       (*gen_jet_phi_).push_back(genJetIter->phi());
       (*gen_jet_eta_).push_back(genJetIter->eta());
       (*gen_jet_energy_).push_back(genJetIter->energy());
       //cout<<"genJetIter->pt(): "<<genJetIter->pt()<<endl;
       //cout<<"genJetIter->eta(): "<<genJetIter->eta()<<endl;
     }
     //genjet1.SetPtEtaPhiE((*gen_jet_pt_)[0],(*gen_jet_eta_)[0],(*gen_jet_phi_)[0],(*gen_jet_energy_)[0]);
     //genjet2.SetPtEtaPhiE((*gen_jet_pt_)[1],(*gen_jet_eta_)[1],(*gen_jet_phi_)[1],(*gen_jet_energy_)[1]);
     //genjj = genjet1+genjet2;
     //(*gen_dijet_mass_).push_back(genjj.M());
     //(*gen_dijet_y_).push_back(genjj.Rapidity());
   }

   //Get Vertices and rho                                                                                    
   edm::Handle<double> rhoHandle;
   iEvent.getByToken(rho_token_,rhoHandle);
   double rho = *rhoHandle; 

   // Get ak8Jets
   edm::Handle<edm::View<pat::Jet>> jets;
   iEvent.getByToken(jetAK8_token_, jets);
   unsigned int collSize=jets->size();
   TLorentzVector jet1, jet2, jj;
   int numbJetsAk8_id=0;
   for (unsigned int ijet=0; ijet<collSize; ijet++) 
     {
       //reco::Jet jet = (*jets)[ijet];;
       const edm::Ptr<pat::Jet> jet = jets->ptrAt(ijet);
       //cout<<jet->pt()<<endl;
       // CMSSW_8_X samples
       //       double pruned_mass = (*jets)[ijet].userFloat("ak8PFJetsCHSPrunedMass");
       //       double tau1         = (*jets)[ijet].userFloat("NjettinessAK8:tau1");
       //       double tau2         = (*jets)[ijet].userFloat("NjettinessAK8:tau2");

       /*
       double NHF  = jet.neutralHadronEnergyFraction();
       double NEMF = jet.neutralEmEnergyFraction();
       double CHF  = jet.chargedHadronEnergyFraction();
       double MUF  = jet.muonEnergyFraction();
       double CEMF = jet.chargedEmEnergyFraction();
       double NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
       double NumNeutralParticles =jet.neutralMultiplicity();
       double CHM      = jet.chargedMultiplicity(); 
       bool looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(jet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jet.eta())>2.4) && abs(jet.eta())<=2.7; 
       */

       // CMSSW_9_4_X
       //double pruned_mass       = (*jets)[ijet].userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass");
       //double tau1         = (*jets)[ijet].userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1");
       //double tau2         = (*jets)[ijet].userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2");
       //From jet toolbox
       double pruned_mass       = (*jets)[ijet].userFloat("ak8PFJetsCHSPrunedMass");
       double tau1         = (*jets)[ijet].userFloat("NjettinessAK8CHS:tau1");
       double tau2         = (*jets)[ijet].userFloat("NjettinessAK8CHS:tau2");
       
       double C_JER=1.0;
       if(isMC == true)
	 {
	   C_JER=calculateJERFactor(jet->pt(),jet->eta(),jet->phi(),jet->energy(),rho,true);
	 }//end of if MC
       

       if( (C_JER*jet->pt())>200&&fabs(jet->eta())<2.4){
       //if(jet->pt()>200&&fabs(jet->eta())<2.4){
	 bool isLepton=isJetLeptonAK8(jet->eta(),jet->phi());
	 if(!isLepton){
	   //(*jet_pt_).push_back(jet->pt());
	   (*jet_pt_).push_back(C_JER*jet->pt());
	   (*jet_phi_).push_back(jet->phi());
	   (*jet_eta_).push_back(jet->eta());
	   (*jet_px_).push_back(jet->px());
	   (*jet_py_).push_back(jet->py());
	   (*jet_pz_).push_back(jet->pz());
	   (*jet_energy_).push_back(C_JER*jet->energy());
	   (*jet_mass_).push_back(pruned_mass);
	   (*jet_tau1_).push_back(tau1);
	   (*jet_tau2_).push_back(tau2);
	   (*jet_vertexz_).push_back(jet->vz());

	   std::vector<std::pair<std::string, float> > btag=jet->getPairDiscri();
	   int size_v=btag.size();
	   for(int i=0;i<size_v;i++){
	     if(btag[i].first=="pfCombinedInclusiveSecondaryVertexV2BJetTags"){
	       //cout<<"btag and result: "<<btag[i].first<<", "<<btag[i].second<<endl;
	       //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
	       if(btag[i].second>0.8838){
		 numbJetsAk8_id++;
	       }//end of btag requirement
	     }//end of requireing looking at pfCombinedInclusiveSecondaryVertexV2BJetTags
	   }//end of looping over btags
	   
	   
	   // Now apply corrections!
	   double pruned_masscorr = 0;
	   double corr = 0;
	   //cout<<"jet pt Hollar: "<<jet->correctedP4(0).pt()<<endl;
	   //cout<<"jet pt Finn: "<<jet->correctedJet(0).pt()<<endl;
	   jecAK8_->setJetEta( jet->correctedJet(0).eta() );
	   jecAK8_->setJetPt ( jet->correctedJet(0).pt() );
	   jecAK8_->setJetE  ( jet->correctedJet(0).energy() );
	   jecAK8_->setJetA  ( jet->correctedJet(0).jetArea() );
	   //jecAK8_->setJetA  ( jet->jetArea() );
	   jecAK8_->setRho   ( rho );
	   jecAK8_->setNPV   ( vertices_->size() );
	   corr = jecAK8_->getCorrection();

	   /*
	   jecAK8_withL1_->setJetEta( jet->correctedJet(0).eta() );
	   jecAK8_withL1_->setJetPt ( jet->correctedJet(0).pt() );
	   jecAK8_withL1_->setJetE  ( jet->correctedJet(0).energy() );
	   jecAK8_withL1_->setJetA  ( jet->correctedJet(0).jetArea() );
	   jecAK8_withL1_->setRho   ( rho );
	   jecAK8_withL1_->setNPV   ( vertices_->size() );
	   double corr_with_L1 = jecAK8_withL1_->getCorrection();
	   
	   cout<<"Correction from Jet Corrector is: "<<corr_with_L1<<endl;
	   cout<<"Uncorrected jet pt: "<<jet->correctedP4(0).pt()<<endl;
	   cout<<"Uncorrected jet pt2: "<<jet->correctedJet(0).pt()<<endl;
	   cout<<"Corrected jet pt: "<<jet->pt()<<endl;
	   */
	   pruned_masscorr = corr*pruned_mass;
	   (*jet_corrmass_).push_back(pruned_masscorr);



	 }}//endof looking at if by lepton and if pt>200 GeV
     }
   *num_bjets_ak8_=numbJetsAk8_id;

  
   //Get MET
   edm::Handle<pat::METCollection> patMET; // PAT      
   iEvent.getByToken(met_token_,patMET);
   const pat::METRef MET(patMET, 0);
   *met_=MET->caloMETPt();
   //*met_=patMET.corPt();
   *met_x_=MET->corPx();
   *met_y_=MET->corPy();
   double METPt=MET->corPt();
   double METphi=MET->phi();
   *met_phi_=MET->phi();

   //Look at proton tracks
   //if(isMC==false||isSignalMC==true||isMC==true){
   // Proton lite tracks
   edm::Handle<std::vector<CTPPSLocalTrackLite> > ppsTracks;
   iEvent.getByToken( pps_token_, ppsTracks );
     
   for ( const auto& trk : *ppsTracks ) 
     {
       const CTPPSDetId detid( trk.getRPId() );
       // transform the raw, 32-bit unsigned integer detId into the TOTEM "decimal" notation
       const unsigned short raw_id = 100*detid.arm()+10*detid.station()+detid.rp();
       //cout<<"pps track x: "<<trk.getX();
       (*pps_track_x_).push_back(trk.getX());
       (*pps_track_y_).push_back(trk.getY());
       (*pps_track_rpid_).push_back(raw_id);
     }


     // Full reco protons
   edm::Handle<vector<reco::ForwardProton>> recoMultiRPProtons;
   iEvent.getByToken(recoProtonsMultiRPToken_, recoMultiRPProtons);
   edm::Handle<vector<reco::ForwardProton>> recoSingleRPProtons;
   iEvent.getByToken(recoProtonsSingleRPToken_, recoSingleRPProtons);
       
   int ismultirp = -999;
   unsigned int decRPId = -999;
   unsigned int armId = -999;
   float th_y = -999;
   float th_x = -999;
   float t = -999;
   float xi = -999;
   float trackx1 = -999.;
   float tracky1 = -999.;
   float trackx2 = -999.;
   float tracky2 = -999.;
   unsigned int trackrpid1 = -999;
   unsigned int trackrpid2 = -999;
   int pixshift1 = -999;
   int pixshift2 = -999;
   float time = -999.;

   // Single-RP
   for (const auto & proton : *recoSingleRPProtons)
     {
       if (proton.validFit())
	 {
	   th_y = proton.thetaY();
	   th_x = proton.thetaX();
	   xi = proton.xi();
	   t = proton.t();
	   time = proton.time(); 

	   trackx1 = (*proton.contributingLocalTracks().begin())->getX();
	   tracky1 = (*proton.contributingLocalTracks().begin())->getY();

	   CTPPSpixelLocalTrackReconstructionInfo pixtrackinfo1 = (*proton.contributingLocalTracks().begin())->getPixelTrackRecoInfo();
	   if(pixtrackinfo1 == CTPPSpixelLocalTrackReconstructionInfo::notShiftedRun || pixtrackinfo1 == CTPPSpixelLocalTrackReconstructionInfo::noShiftedPlanes ||
	      pixtrackinfo1 == CTPPSpixelLocalTrackReconstructionInfo::invalid)
	     pixshift1 = 0;
	   else
	     pixshift1 = 1;

	   CTPPSDetId rpId((*proton.contributingLocalTracks().begin())->getRPId());
	   decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();
	   ismultirp = 0;

	   (*proton_xi_).push_back(xi);
	   //cout<<"Single Pot proton_xi_"<<xi<<endl;
	   (*proton_thy_).push_back(th_y);
	   (*proton_thx_).push_back(th_x);
	   (*proton_t_).push_back(t);
	   (*proton_ismultirp_).push_back(ismultirp);
	   (*proton_rpid_).push_back(decRPId);
	   (*proton_arm_).push_back(armId);
	   
	   (*proton_time_).push_back(time);
	   (*proton_trackx1_).push_back(trackx1);
	   (*proton_tracky1_).push_back(tracky1);
	   (*proton_trackpixshift1_).push_back(pixshift1);
	   (*proton_rpid1_).push_back(decRPId);
	 }
     }

     // Multi-RP
   for (const auto & proton : *recoMultiRPProtons)
     {
       if (proton.validFit())
	 {
	   th_y = proton.thetaY();
	   th_x = proton.thetaX();
	   xi = proton.xi();
	   t = proton.t();
	   time = proton.time();

	   int ij=0;
	   for (const auto &tr : proton.contributingLocalTracks())
	     {
	       CTPPSDetId rpIdJ(tr->getRPId());
	       unsigned int rpDecIdJ = rpIdJ.arm()*100 + rpIdJ.station()*10 + rpIdJ.rp();
	       
	       CTPPSpixelLocalTrackReconstructionInfo pixtrackinfo = (*proton.contributingLocalTracks().begin())->getPixelTrackRecoInfo();

	       if(ij == 0)
		 {
		   trackx1 = tr->getX();
		   tracky1 = tr->getY();
		   trackrpid1 = rpDecIdJ;
		   armId = rpIdJ.arm();
		   if(pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::notShiftedRun || pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::noShiftedPlanes ||
		      pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::invalid)
		     pixshift1 = 0;
		   else
		     pixshift1 = 1;
		 }
	       if(ij == 1)
		 {
		   trackx2 = tr->getX();
		   tracky2 = tr->getY();
		   trackrpid2 = rpDecIdJ;
		   if(pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::notShiftedRun || pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::noShiftedPlanes ||
		      pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::invalid)
		     pixshift2 = 0;
		   else
		     pixshift2 = 1;
		 }
	       ij++;
	     }

	   ismultirp = 1;

	   (*proton_xi_).push_back(xi);
	   //cout<<"Multi Pot proton_xi_"<<xi<<endl;
	   (*proton_thy_).push_back(th_y);
	   (*proton_thx_).push_back(th_x);
	   (*proton_t_).push_back(t);
	   (*proton_ismultirp_).push_back(ismultirp);
	   (*proton_rpid_).push_back(decRPId);
	   (*proton_arm_).push_back(armId);
	   
	   (*proton_time_).push_back(time);
	   (*proton_trackx1_).push_back(trackx1);
	   (*proton_tracky1_).push_back(tracky1);
	   (*proton_trackx2_).push_back(trackx2);
	   (*proton_tracky2_).push_back(tracky2);
	   (*proton_trackpixshift1_).push_back(pixshift1);
	   (*proton_trackpixshift2_).push_back(pixshift2);
	   (*proton_rpid1_).push_back(trackrpid1);
	   (*proton_rpid2_).push_back(trackrpid2);
	 }
     }
   //}

   *nVertices_=-1;
   *nVertices_=vertices_->size();

   // Fill pileup reweighting info if running on MC
   *pileupWeight_=1;
   *mc_pu_trueinteractions_=-999;
   if(isMC == true)
   {
     float trueInteractions=0;
     edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
     iEvent.getByToken(pu_token_, PupInfo);
     std::vector<PileupSummaryInfo>::const_iterator PVI;
     //cout<<"True num interactions is: "<<PupInfo->begin()->getTrueNumInteractions()<<endl;                                                              
     trueInteractions=PupInfo->begin()->getTrueNumInteractions();
     *pileupWeight_ = LumiWeights->weight( trueInteractions );
     *mc_pu_trueinteractions_ = trueInteractions;
   }


   /*
   if(isMC == true){
     
     edm::Handle<GenEventInfoProduct> genEvtInfo; 
     iEvent.getByLabel( "generator", genEvtInfo );
     
     //std::vector<double>& evtWeights = genEvtInfo->weights();
     double theWeight = genEvtInfo->weight();
     //cout<<"theWeight: "<<theWeight<<endl;
     *mcWeight_=theWeight;
     
   }
   */

   // Fill GEN infor if running on MC
   if(isMC == true)
     {
       edm::Handle<reco::GenParticleCollection> genP;
       iEvent.getByLabel("prunedGenParticles",genP);
       TLorentzVector Wplus;
       TLorentzVector Wminus;
       for (reco::GenParticleCollection::const_iterator mcIter=genP->begin(); mcIter != genP->end(); mcIter++ ) {

	 if((fabs(mcIter->pdgId()) == 24) && (mcIter->status() == 22))
	 //if((fabs(mcIter->pdgId()) == 24))
	   {
	     (*gen_W_pt_).push_back(mcIter->pt());
	     (*gen_W_charge_).push_back(mcIter->charge());
	     //cout<<"W pt is:"<<mcIter->pt()<<endl;
	     //cout<<"W status is:"<<mcIter->status()<<endl;
	     //cout<<"W mother is:"<<mcIter->mother()->pdgId()<<endl;
	     // cout<<"W charge is:"<<mcIter->charge()<<endl;
	     if (mcIter->charge() > 0){
	       Wplus=TLorentzVector(mcIter->px(),mcIter->py(),mcIter->pz(),mcIter->energy());
	     }
	     else{
	       Wminus=TLorentzVector(mcIter->px(),mcIter->py(),mcIter->pz(),mcIter->energy());
	     }
	   }

	 if((mcIter->pdgId() == 2212) && (fabs(mcIter->pz()) > 650) && (mcIter->status() == 1))
	 //if((mcIter->pdgId() == 2212) && (mcIter->status() == 1))
	   {
	     double thexi = 1 - ((mcIter->energy())/(13000.0/2.0));
	     double thet = -(std::pow(mcIter->pt(), 2));
	     double thepz = mcIter->pz();
	     double thepx = mcIter->px();
	     double thepy = mcIter->py();
	     double theenergy= mcIter->energy();
	     //double vx = mcIter->vx(), vy = mcIter->vy(), vz = mcIter->vz();
	     //cout<<"vx: "<<vx<<", "<<"vy: "<<vy<<", "<<"vz: "<<vz<<endl;
	     (*gen_proton_xi_).push_back(thexi);
	     (*gen_proton_t_).push_back(thet);
	     (*gen_proton_pz_).push_back(thepz);
	     (*gen_proton_py_).push_back(thepy);
	     (*gen_proton_px_).push_back(thepx);
	     (*gen_proton_energy_).push_back(theenergy);
	   }
       }//end of looping over gen particles
       *genWWpt_=(Wplus+Wminus).Pt();
     }//is MC

   bool passes_muon=false;
   bool passes_electron=false;
   bool passes_dimuon=false;
   bool passes_dielectron=false;
   if((*muon_pt_).size()==1&&(*jet_eta_).size()==1&&numMuLoose==1&&num_ele_veto==0&&passTrigger_mu){ 
     if((*muon_pt_)[0]>50){       passes_muon=true;   }
   }
   if((*electron_pt_).size()==1&&(*jet_eta_).size()==1&&num_ele_veto==1&&numMuLoose==0&&passTrigger_e){
     if((*electron_pt_)[0]>50){       passes_electron=true;   }
   }
   if((*muon_pt_).size()==2&&(*jet_eta_).size()==1&&numMuLoose==2&&num_ele_veto==0&&passTrigger_mu){     passes_dimuon=true;   }
   if((*electron_pt_).size()==2&&(*jet_eta_).size()==1&&numMuLoose==0&&num_ele_veto==2&&passTrigger_e){     passes_dielectron=true;   }
   bool passFatJet=false;
   if(passes_muon||passes_electron||passes_dimuon||passes_dielectron){passFatJet=true;}

   //Only store samples that pass Fat Jet
   if(passFatJet){
     edm::Handle<edm::View<pat::Jet> > jetColl; // PAT      
     iEvent.getByLabel("slimmedJetsJetId", jetColl);
     int numbJets_id=0;
     int numJetsAK4_id=0;
     double C_JER_AK4=1.0;
     for (unsigned int i=0; i<jetColl->size(); i++) {
       const edm::Ptr<pat::Jet> jet_ak4 = jetColl->ptrAt(i);
       if((*jet_eta_).size()==1){
	 if(isMC==true){
	   C_JER_AK4=calculateJERFactor(jet_ak4->pt(),jet_ak4->eta(),jet_ak4->phi(),jet_ak4->energy(),rho,false);}
	 //cout<<"AK4 jet pt: "<<jet->pt()<<endl;
	 if((jet_ak4->pt()*C_JER_AK4)>30&&fabs(jet_ak4->eta())<2.4){
	   bool isLepton=isJetLepton(jet_ak4->eta(),jet_ak4->phi());
	   if(!isLepton){
	     //numJets_id++;  
	     //if(jet->phi(),jet->eta())
	     double deltaR=sqrt( ((*jet_eta_)[0]-jet_ak4->eta())*((*jet_eta_)[0]-jet_ak4->eta()) + ((*jet_phi_)[0]-jet_ak4->phi())*((*jet_phi_)[0]-jet_ak4->phi()) );
	     if(deltaR<0.8){
	      continue;
	     }
	     else{numJetsAK4_id++;}
	     std::vector<std::pair<std::string, float> > btag=jet_ak4->getPairDiscri();
	     int size_v=btag.size();
	     for(int i=0;i<size_v;i++){
	       if(btag[i].first=="pfCombinedInclusiveSecondaryVertexV2BJetTags"){
		 //cout<<"btag and result: "<<btag[i].first<<", "<<btag[i].second<<endl;
		 //medium is 0.8484, loose is 0.5426, tight is 0.9535
		 //if(btag[i].second>0.9535){
		 if(btag[i].second>0.9693){
		   numbJets_id++;
		 }//end of btag requirement
	     }//end of requireing looking at pfCombinedInclusiveSecondaryVertexV2BJetTags
	     }//end of looping over btags

	   }//end of looking at if it is a lepton
	   //cout<<"Jet energy, phi, eta: "<<jet->energy()<<", "<<jet->phi()<<", "<<jet->eta()<<endl;
	   //cout<<"Jet corrected energy: "<<jet->energy()<<endl;
	   //cout<<"Uncorrected energy: "<<jet->correctedP4(0).E()<<endl;
	   //cout<<"Jet corrected pt: "<<jet->pt()<<endl;
	   //cout<<"Uncorrected pt: "<<jet->correctedP4(0).Pt()<<endl;
	   //std::vector< std::string > jec_levels = jet_ak4->availableJECLevels();
	 //int size = jec_levels.size();
	 //for(int i=0;i<size;i++){
	 //  cout<<jec_levels[i]<<endl;
	 // }
	 }//end of requirement of 30 geV
       }//Need at least one fat jet to even look at this
     }//end of looping over jet collection
     *num_bjets_ak4_=numbJets_id;
     *num_jets_ak4_=numJetsAK4_id;

     edm::Handle<edm::View<pat::PackedCandidate> > pfCands; // PAT      
     iEvent.getByToken(pfcand_token_, pfCands);
     int count_pv=0;
     int count_pv_noDRl=0;
     for (const pat::PackedCandidate &pfcand : *pfCands) { 
       //if(fabs(pfcand.pdgId())==211&&pfcand.fromPV(0)==3){
       if((*jet_eta_).size()>0){
	 if(pfcand.fromPV(0)==3&&(fabs(pfcand.pdgId())==13||fabs(pfcand.pdgId())==211||fabs(pfcand.pdgId()) == 11 )){
	   double deltaR=sqrt(  (pfcand.eta()-(*jet_eta_)[0])*(pfcand.eta()-(*jet_eta_)[0]) + deltaPhi(pfcand.phiAtVtx(),(*jet_phi_)[0])*deltaPhi(pfcand.phiAtVtx(),(*jet_phi_)[0]));
	   if (deltaR>0.8){
	     //cout<<"z distance: "<<pfcand.vertex().z()-vtx->position().z()-vtx->position().z()<<endl;
	     //cout<<"vtx->position().z(): "<<vtx->position().z()<<endl;
	     //if ((abs(pfcand.vertex().z()-vtx->position().z())<0.02) && (abs(pfcand.vertex().x()-vtx->position().x())<0.006) && (abs(pfcand.vertex().y()-vtx->position().y())<0.006)){
	     double deltaR_mu=1000;
	     double deltaR_mu_2=1000;
	     double deltaR_e=1000;
	     double deltaR_e_2=1000;
	     if(passes_muon||passes_dimuon){
	       deltaR_mu=sqrt(  (pfcand.eta()-(*muon_eta_)[0])*(pfcand.eta()-(*muon_eta_)[0]) + deltaPhi(pfcand.phiAtVtx(),(*muon_phi_)[0])*deltaPhi(pfcand.phiAtVtx(),(*muon_phi_)[0]));
	       if(passes_dimuon) {deltaR_mu_2=sqrt(  (pfcand.eta()-(*muon_eta_)[1])*(pfcand.eta()-(*muon_eta_)[1]) + deltaPhi(pfcand.phiAtVtx(),(*muon_phi_)[1])*deltaPhi(pfcand.phiAtVtx(),(*muon_phi_)[1]));}
	       if(deltaR_mu>0.3&&passes_muon){count_pv=count_pv+1;}
	       if(deltaR_mu>0.3&&deltaR_mu_2>0.3&&passes_dimuon){count_pv=count_pv+1;}
	       if(!(deltaR_mu<0.3&&fabs(pfcand.pdgId())==13)&&passes_muon){count_pv_noDRl=count_pv_noDRl+1;}
	       if(!(deltaR_mu<0.3&&fabs(pfcand.pdgId())==13)&&!(deltaR_mu_2<0.3&&fabs(pfcand.pdgId())==13)&&passes_dimuon){count_pv_noDRl=count_pv_noDRl+1;}
	     }
	     if(passes_electron||passes_dielectron){
	       deltaR_e=sqrt(  (pfcand.eta()-(*electron_eta_)[0])*(pfcand.eta()-(*electron_eta_)[0]) + deltaPhi(pfcand.phiAtVtx(),(*electron_phi_)[0])*deltaPhi(pfcand.phiAtVtx(),(*electron_phi_)[0]));
	       if(passes_dielectron){deltaR_e_2=sqrt(  (pfcand.eta()-(*electron_eta_)[1])*(pfcand.eta()-(*electron_eta_)[1]) + deltaPhi(pfcand.phiAtVtx(),(*electron_phi_)[1])*deltaPhi(pfcand.phiAtVtx(),(*electron_phi_)[1]));}
	       if(deltaR_e>0.3&&passes_electron){count_pv=count_pv+1;}
	       if(deltaR_e>0.3&&deltaR_e_2>0.3&&passes_dielectron){count_pv=count_pv+1;}
	       if(!(deltaR_e<0.3&&fabs(pfcand.pdgId())==11)&&passes_electron){count_pv_noDRl=count_pv_noDRl+1;}
	       if(!(deltaR_e<0.3&&fabs(pfcand.pdgId())==11)&&!(deltaR_e_2<0.3&&fabs(pfcand.pdgId())==11)&&passes_dielectron){count_pv_noDRl=count_pv_noDRl+1;}
	     }
	     //}//end of requiring pfcand vertex close to primary vertex
	   }
	 }
       }
     }//end over loop over pf candidates
     *pfcand_nextracks_=count_pv;
     *pfcand_nextracks_noDRl_=count_pv_noDRl;
     

     TLorentzVector LeptonP4;
     if(passes_muon||passes_dimuon){
       LeptonP4=TLorentzVector((*muon_px_)[0],(*muon_py_)[0],(*muon_pz_)[0],(*muon_e_)[0]);}
     if(passes_electron||passes_dielectron){
       LeptonP4=TLorentzVector((*electron_px_)[0],(*electron_py_)[0],(*electron_pz_)[0],(*electron_e_)[0]);}
     //
     //Neutrino 
     TLorentzVector SecondLeptonP4;
     if(passes_muon||passes_electron){
       TLorentzVector Neutrino_nominal = Nu4Momentum(LeptonP4, METPt, METphi);   // Neutrino Method   
       double E_nu=sqrt(MET->corPx()*MET->corPx()+MET->corPy()*MET->corPy()+Neutrino_nominal.Pz()*Neutrino_nominal.Pz());
       SecondLeptonP4 = TLorentzVector(*met_x_,*met_y_,Neutrino_nominal.Pz(),E_nu); 
     }
     if(passes_dimuon){
       SecondLeptonP4 = TLorentzVector((*muon_px_)[1],(*muon_py_)[1],(*muon_pz_)[1],(*muon_e_)[1]); 
     }
     if(passes_dielectron){
       SecondLeptonP4 = TLorentzVector((*electron_px_)[1],(*electron_py_)[1],(*electron_pz_)[1],(*electron_e_)[1]); 
     }

     //W lep with neutrino - W leptonic
     TLorentzVector WLeptonic = LeptonP4 + SecondLeptonP4;
     double WLeptonic_Pt  = WLeptonic.Pt();
     double WLeptonic_M   = WLeptonic.M();
     double WLeptonic_phi = WLeptonic.Phi();
     TLorentzVector p4SumJets;
     //p4SumJets = TLorentzVector((*jet_px_)[0],(*jet_py_)[0],(*jet_pz_)[0], (*jet_energy_)[0]);
     p4SumJets.SetPtEtaPhiE((*jet_pt_)[0],(*jet_eta_)[0],(*jet_phi_)[0], (*jet_energy_)[0]);
     double WJM = p4SumJets.M();
     double WJphi = p4SumJets.Phi();
       //
       // WW
       //
     TLorentzVector WW_p4_nu  = WLeptonic + p4SumJets;               // Wlep with Nu
     double MWW_Nu  = WW_p4_nu.M();
     *recoMWhad_=WJM;

     *recoMWlep_=WLeptonic_M;
     double dphi=deltaPhi(WLeptonic_phi,WJphi);
     *dphiWW_=dphi;
     *WLeptonicPhi_=WLeptonic_phi;
     *WLeptonicPt_=WLeptonic_Pt;
     *WLeptonicEta_=WLeptonic.Rapidity();
     *recoMWW_=MWW_Nu;
     *recoWWpt_=WW_p4_nu.Pt();
     *recoRapidityWW_=WW_p4_nu.Rapidity();

     tree_->Fill();

   }
   else{
     //if(passes_dimuon||passes_dielectron){ tree_->Fill();}
     //if((*muon_pt_).size()==2&&numMuLoose==2&&num_ele_veto==0&&passTrigger_mu){ tree_->Fill();}
     //if((*electron_pt_).size()==2&&numMuLoose==0&&num_ele_veto==1&&passTrigger_e){tree_->Fill();}
   }
   

   (*muon_pt_).clear();
   (*muon_eta_).clear();
   (*muon_phi_).clear();
   (*muon_px_).clear();
   (*muon_py_).clear();
   (*muon_pz_).clear();
   (*muon_e_).clear();
   (*muon_charge_).clear();
   (*muon_iso_).clear();
   (*muon_dxy_).clear();
   (*muon_dz_).clear();

   (*electron_pt_).clear();
   (*electron_eta_).clear();
   (*electron_phi_).clear();
   (*electron_dxy_).clear();
   (*electron_dz_).clear();
   (*electron_px_).clear();
   (*electron_py_).clear();
   (*electron_pz_).clear();
   (*electron_e_).clear();
   (*electron_charge_).clear();


   (*jet_pt_).clear();
   (*jet_px_).clear();
   (*jet_py_).clear();
   (*jet_pz_).clear();
   (*jet_energy_).clear();
   (*jet_phi_).clear();
   (*jet_eta_).clear();
   (*jet_mass_).clear();
   (*jet_tau1_).clear();
   (*jet_tau2_).clear();
   (*jet_corrmass_).clear();
   (*jet_vertexz_).clear();
   (*jet_jer_res_).clear();
   (*jet_jer_sf_).clear();
   (*jet_jer_sfup_).clear();
   (*jet_jer_sfdown_).clear();

   (*pps_track_x_).clear();
   (*pps_track_y_).clear();
   (*pps_track_rpid_).clear();

   (*gen_W_pt_).clear();
   (*gen_W_charge_).clear();

   (*gen_proton_px_).clear();
   (*gen_proton_py_).clear();
   (*gen_proton_pz_).clear();
   (*gen_proton_energy_).clear();
   (*gen_proton_xi_).clear();
   (*gen_proton_t_).clear();
   
   (*proton_xi_).clear();
   (*proton_thy_).clear();
   (*proton_thx_).clear();
   (*proton_t_).clear();
   (*proton_ismultirp_).clear();
   (*proton_rpid_).clear();
   (*proton_arm_).clear();
   (*proton_time_).clear();
   (*proton_trackx1_).clear();
   (*proton_tracky1_).clear();
   (*proton_trackx2_).clear();
   (*proton_tracky2_).clear();
   (*proton_trackpixshift1_).clear();
   (*proton_trackpixshift2_).clear();
   (*proton_rpid1_).clear();
   (*proton_rpid2_).clear();

   (*gen_jet_pt_).clear();
   (*gen_jet_phi_).clear();
   (*gen_jet_eta_).clear();
   (*gen_jet_energy_).clear();

   (*hlt_).clear();



}


bool Ntupler::isJetLeptonAK8(double jet_eta, double jet_phi)
{
  
  for(uint i=0;i<(*muon_eta_).size();i++){
    double eta=(*muon_eta_)[i];
    double phi=(*muon_phi_)[i];
    double deltaR=sqrt((eta-jet_eta)*(eta-jet_eta)+(phi-jet_phi)*(phi-jet_phi));
    if(deltaR<1.0){
      return true;
    }
  }
  
  for(uint i=0;i<(*electron_eta_).size();i++){
    double eta=(*electron_eta_)[i];
    double phi=(*electron_phi_)[i];
    double deltaR=sqrt((eta-jet_eta)*(eta-jet_eta)+(phi-jet_phi)*(phi-jet_phi));
    if(deltaR<1.0){
      return true;
    }
  }

  return false;

}

bool Ntupler::isJetLepton(double jet_eta, double jet_phi)
{
  
  for(uint i=0;i<(*muon_eta_).size();i++){
    double eta=(*muon_eta_)[i];
    double phi=(*muon_phi_)[i];
    double deltaR=sqrt((eta-jet_eta)*(eta-jet_eta)+(phi-jet_phi)*(phi-jet_phi));
    if(deltaR<0.3){
      return true;
    }
  }
  
  for(uint i=0;i<(*electron_eta_).size();i++){
    double eta=(*electron_eta_)[i];
    double phi=(*electron_phi_)[i];
    double deltaR=sqrt((eta-jet_eta)*(eta-jet_eta)+(phi-jet_phi)*(phi-jet_phi));
    if(deltaR<0.3){
      return true;
    }
  }

  return false;

}

double Ntupler::calculateJERFactor(double jet_pt,double jet_eta,double jet_phi,double jet_energy,double rho,bool isAK8)
{
  double C_JER=1.;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor resolution_sf;
  if (isAK8 ==true){
    resolution = JME::JetResolution(jerAK8chsName_res_);
    resolution_sf = JME::JetResolutionScaleFactor(jerAK8chsName_sf_);
  }
  else{
    resolution = JME::JetResolution(jerAK4chsName_res_);
    resolution_sf = JME::JetResolutionScaleFactor(jerAK4chsName_sf_);
  }
  JME::JetParameters parameters;
  parameters.setJetPt(jet_pt);
  parameters.setJetEta(jet_eta);
  parameters.setRho(rho);
      
  float jer_res= resolution.getResolution(parameters);
  float jer_sf = resolution_sf.getScaleFactor(parameters);
  float jer_sf_up = resolution_sf.getScaleFactor(parameters, Variation::UP);
  float jer_sf_down = resolution_sf.getScaleFactor(parameters, Variation::DOWN);
	   
  if (isAK8 ==true){
    (*jet_jer_res_).push_back(jer_res);
    (*jet_jer_sf_).push_back(jer_sf);
    (*jet_jer_sfup_).push_back(jer_sf_up);
    (*jet_jer_sfdown_).push_back(jer_sf_down);
  }
  TLorentzVector recojtmp, genjtmp;
  TRandom3 randomSrc;
  recojtmp.SetPtEtaPhiE(jet_pt,jet_eta,jet_phi,jet_energy);
  
  int matchedgen=0;
  int indmatchedgen=-1;
  //Jet smearing calculation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
  //loop over gen jets
  for(uint ig=0; ig<(*gen_jet_pt_).size(); ig++){
    genjtmp.SetPtEtaPhiE((*gen_jet_pt_).at(ig),(*gen_jet_eta_).at(ig),(*gen_jet_phi_).at(ig),(*gen_jet_energy_).at(ig));

    if (isAK8==true){
      if( (recojtmp.DeltaR(genjtmp) < (0.8/2.)) && (fabs(recojtmp.Pt() - genjtmp.Pt())<(3*jer_res*recojtmp.Pt()) ) ){
	matchedgen=1; 
	indmatchedgen=ig;
      } // 0.8 is cone radius
    }//end of looking for AK8 gen jet
    else{
      if( (recojtmp.DeltaR(genjtmp) < (0.4/2.)) && (fabs(recojtmp.Pt() - genjtmp.Pt())<(3*jer_res*recojtmp.Pt()) ) ){
	matchedgen=1; 
	indmatchedgen=ig;
      } // 0.4 is cone radius
    }//end of looking for AK4 gen jet

  }
  
  if(matchedgen == 1) { C_JER = 1 + (jer_sf -1 )*( (recojtmp.Pt() - (*gen_jet_pt_).at(indmatchedgen)) / recojtmp.Pt() );
    if(C_JER < 0) {C_JER = 0;}
  }
  else       
    {	   C_JER = 1 + randomSrc.Gaus(0, jer_res)*(sqrt(max(jer_sf*jer_sf - 1., 0.)));	 }
  
  return C_JER;


}


// ------------ method called once each job just before starting event loop  ------------
void 
Ntupler::beginJob()
{

  edm::Service<TFileService> fs; 
  tree_=fs->make<TTree>("SlimmedNtuple","SlimmedNtuple");

  muon_pt_ = new std::vector<float>;
  muon_eta_ = new std::vector<float>;
  muon_phi_ = new std::vector<float>;
  muon_px_ = new std::vector<float>;
  muon_py_ = new std::vector<float>;
  muon_pz_ = new std::vector<float>;
  muon_e_ = new std::vector<float>;
  muon_charge_ = new std::vector<float>;
  muon_iso_ = new std::vector<float>;
  muon_dxy_ = new std::vector<float>;
  muon_dz_ = new std::vector<float>;

  electron_pt_ = new std::vector<float>;
  electron_eta_ = new std::vector<float>;
  electron_phi_ = new std::vector<float>;
  electron_dxy_ = new std::vector<float>;
  electron_dz_ = new std::vector<float>;
  electron_px_ = new std::vector<float>;
  electron_py_ = new std::vector<float>;
  electron_pz_ = new std::vector<float>;
  electron_e_ = new std::vector<float>;
  electron_charge_ = new std::vector<float>;

  met_ = new float;
  met_x_ = new float;
  met_y_ = new float;
  met_phi_ = new float;
  pfcand_nextracks_ = new int;
  pfcand_nextracks_noDRl_ = new int;
  num_bjets_ak8_ = new int;
  num_bjets_ak4_ = new int;
  num_jets_ak4_ = new int;

  recoMWhad_ = new float;
  recoMWlep_ = new float;
  recoMWW_ = new float;  
  recoWWpt_ = new float;
  genWWpt_ = new float;
  recoRapidityWW_ = new float;
  dphiWW_ = new float;
  WLeptonicPt_ = new float;
  WLeptonicPhi_ = new float;
  WLeptonicEta_ = new float;

  //ecalBadCalFilter_ = new bool;

  jet_pt_ = new std::vector<float>;
  jet_px_ = new std::vector<float>;
  jet_py_ = new std::vector<float>;
  jet_pz_ = new std::vector<float>;
  jet_energy_ = new std::vector<float>;
  jet_phi_ = new std::vector<float>;
  jet_eta_ = new std::vector<float>;
  jet_mass_ = new std::vector<float>;
  jet_tau1_ = new std::vector<float>;
  jet_tau2_ = new std::vector<float>;
  jet_corrmass_ = new std::vector<float>;
  jet_vertexz_ = new std::vector<float>;
  jet_jer_res_ = new std::vector<float>;
  jet_jer_sf_ = new std::vector<float>;
  jet_jer_sfup_ = new std::vector<float>;
  jet_jer_sfdown_ = new std::vector<float>;

  pps_track_x_ = new std::vector<float>;
  pps_track_y_ = new std::vector<float>;
  pps_track_rpid_ = new std::vector<int>;

  gen_W_pt_ = new std::vector<float>;
  gen_W_charge_ = new std::vector<float>;

  gen_proton_px_ = new std::vector<float>;
  gen_proton_py_ = new std::vector<float>;
  gen_proton_pz_ = new std::vector<float>;
  gen_proton_energy_ = new std::vector<float>;
  gen_proton_xi_ = new std::vector<float>;
  gen_proton_t_ = new std::vector<float>;

  proton_xi_ = new std::vector<float>;
  proton_thy_ = new std::vector<float>;
  proton_thx_ = new std::vector<float>;
  proton_t_ = new std::vector<float>;
  proton_ismultirp_ = new std::vector<int>;
  proton_rpid_ = new std::vector<int>;
  proton_arm_ = new std::vector<int>;

  proton_time_ = new std::vector<float>;
  proton_trackx1_ = new std::vector<float>;
  proton_tracky1_ = new std::vector<float>;
  proton_trackx2_ = new std::vector<float>;
  proton_tracky2_ = new std::vector<float>;
  proton_trackpixshift1_ = new std::vector<int>;
  proton_trackpixshift2_ = new std::vector<int>;
  proton_rpid1_ = new std::vector<int>;
  proton_rpid2_ = new std::vector<int>;

  gen_jet_pt_ = new std::vector<float>;
  gen_jet_phi_ = new std::vector<float>;
  gen_jet_eta_ = new std::vector<float>;
  gen_jet_energy_ = new std::vector<float>;

  hlt_ = new std::vector<string>;

  ev_ = new long int;
  run_ = new int;
  lumiblock_ = new int;
  crossingAngle_ = new float;
  betaStar_ = new float;
  instLumi_ = new float;
  nVertices_ = new int;
  pileupWeight_ = new float;
  mc_pu_trueinteractions_ = new float;
  mcWeight_ = new float;

  tree_->Branch("muon_pt",&muon_pt_);
  tree_->Branch("muon_eta",&muon_eta_);
  tree_->Branch("muon_phi",&muon_phi_);
  tree_->Branch("muon_px",&muon_px_);
  tree_->Branch("muon_py",&muon_py_);
  tree_->Branch("muon_pz",&muon_pz_);
  tree_->Branch("muon_e",&muon_e_);
  tree_->Branch("muon_charge",&muon_charge_);
  tree_->Branch("muon_iso",&muon_iso_);
  tree_->Branch("muon_dxy",&muon_dxy_);
  tree_->Branch("muon_dz",&muon_dz_);
  tree_->Branch("electron_pt",&electron_pt_);
  tree_->Branch("electron_eta",&electron_eta_);
  tree_->Branch("electron_phi",&electron_phi_);
  tree_->Branch("electron_dxy",&electron_dxy_);
  tree_->Branch("electron_dz",&electron_dz_);
  tree_->Branch("electron_px",&electron_px_);
  tree_->Branch("electron_py",&electron_py_);
  tree_->Branch("electron_pz",&electron_pz_);
  tree_->Branch("electron_e",&electron_e_);
  tree_->Branch("electron_charge",&electron_charge_);
  tree_->Branch("met",met_,"met/f");
  tree_->Branch("met_x",met_x_,"met_x/f");
  tree_->Branch("met_y",met_y_,"met_y/f");
  tree_->Branch("met_phi",met_phi_,"met_phi/f");
  tree_->Branch("pfcand_nextracks",pfcand_nextracks_,"pfcand_nextracks/I");
  tree_->Branch("pfcand_nextracks_noDRl",pfcand_nextracks_noDRl_,"pfcand_nextracks_noDRl/I");
  tree_->Branch("num_bjets_ak8",num_bjets_ak8_,"num_bjets_ak8/I");
  tree_->Branch("num_bjets_ak4",num_bjets_ak4_,"num_bjets_ak4/I");
  tree_->Branch("num_jets_ak4",num_jets_ak4_,"num_jets_ak4/I");
  tree_->Branch("recoMWhad",recoMWhad_,"recoMWhad/f");
  tree_->Branch("recoMWlep",recoMWlep_,"recoMWlep/f");
  tree_->Branch("recoMWW",recoMWW_,"recoMWW/f");
  tree_->Branch("recoRapidityWW",recoRapidityWW_,"recoRapidityWW/f");
  tree_->Branch("recoWWpt",recoWWpt_,"recoWWpt_/f");
  tree_->Branch("genWWpt",genWWpt_,"genWWpt_/f");
  tree_->Branch("dphiWW",dphiWW_,"dphiWW/f");
  tree_->Branch("WLeptonicPt",WLeptonicPt_,"WLeptonicPt/f");
  tree_->Branch("WLeptonicPhi",WLeptonicPhi_,"WLeptonicPhi/f");
  tree_->Branch("WLeptonicEta",WLeptonicEta_,"WLeptonicEta/f");

  tree_->Branch("jet_pt",&jet_pt_);
  tree_->Branch("jet_px",&jet_px_);
  tree_->Branch("jet_py",&jet_py_);
  tree_->Branch("jet_pz",&jet_pz_);
  tree_->Branch("jet_energy",&jet_energy_);
  tree_->Branch("jet_phi",&jet_phi_);
  tree_->Branch("jet_eta",&jet_eta_);
  tree_->Branch("jet_mass",&jet_mass_);
  tree_->Branch("jet_tau1",&jet_tau1_);
  tree_->Branch("jet_tau2",&jet_tau2_);
  tree_->Branch("jet_corrmass",&jet_corrmass_);
  tree_->Branch("jet_vertexz",&jet_vertexz_);
  tree_->Branch("jet_jer_res",&jet_jer_res_);
  tree_->Branch("jet_jer_sf",&jet_jer_sf_);
  tree_->Branch("jet_jer_sfup",&jet_jer_sfup_);
  tree_->Branch("jet_jer_sfdown",&jet_jer_sfdown_);
  tree_->Branch("pps_track_x",&pps_track_x_);
  tree_->Branch("pps_track_y",&pps_track_y_);
  tree_->Branch("pps_track_rpid",&pps_track_rpid_);
  tree_->Branch("hlt",&hlt_);
  tree_->Branch("gen_W_pt",&gen_W_pt_);
  tree_->Branch("gen_W_charge",&gen_W_charge_);
  tree_->Branch("gen_proton_px",&gen_proton_px_);
  tree_->Branch("gen_proton_py",&gen_proton_py_);
  tree_->Branch("gen_proton_pz",&gen_proton_pz_);
  tree_->Branch("gen_proton_energy",&gen_proton_energy_);
  tree_->Branch("gen_proton_xi",&gen_proton_xi_);
  tree_->Branch("gen_proton_t",&gen_proton_t_);
  tree_->Branch("proton_xi",&proton_xi_);
  tree_->Branch("proton_thx",&proton_thx_);
  tree_->Branch("proton_thy",&proton_thy_);
  tree_->Branch("proton_t",&proton_t_);
  tree_->Branch("proton_ismultirp_",&proton_ismultirp_);
  tree_->Branch("proton_rpid",&proton_rpid_);
  tree_->Branch("proton_arm",&proton_arm_);
  tree_->Branch("proton_time",&proton_time_);
  tree_->Branch("proton_trackx1",&proton_trackx1_);
  tree_->Branch("proton_tracky1",&proton_tracky1_);
  tree_->Branch("proton_trackx2",&proton_trackx2_);
  tree_->Branch("proton_tracky2",&proton_tracky2_);
  tree_->Branch("proton_trackpixshift1",&proton_trackpixshift1_);
  tree_->Branch("proton_trackpixshift2",&proton_trackpixshift2_);
  tree_->Branch("proton_rpid1",&proton_rpid1_);
  tree_->Branch("proton_rpid2",&proton_rpid2_);
  tree_->Branch("gen_jet_pt",&gen_jet_pt_);
  tree_->Branch("gen_jet_phi",&gen_jet_phi_);
  tree_->Branch("gen_jet_eta",&gen_jet_eta_);
  tree_->Branch("gen_jet_energy",&gen_jet_energy_);
  tree_->Branch("nVertices",nVertices_,"nVertices/i");
  tree_->Branch("pileupWeight",pileupWeight_,"pileupWeight/f");
  tree_->Branch("mc_pu_trueinteractions",mc_pu_trueinteractions_,"mc_pu_trueinteractions/f");
  tree_->Branch("mcWeight",mcWeight_,"mcWeight/f");
  tree_->Branch("run",run_,"run/I");
  tree_->Branch("event",ev_,"event/L");
  tree_->Branch("lumiblock",lumiblock_,"lumiblock/I");
  tree_->Branch("crossingAngle",crossingAngle_,"crossingAngle/f");
  tree_->Branch("betaStar",betaStar_,"betaStar/f");
  tree_->Branch("instLumi",instLumi_,"instLumi/f");
  //tree_->Branch("ecalBadCalFilter",ecalBadCalFilter_,"ecalBadCalFilter/O");


}

// ------------ method called once each jo bjust after ending the event loop  ------------
void 
Ntupler::endJob() 
{

  delete muon_px_;
  delete muon_pt_;
  delete muon_eta_;
  delete muon_phi_;
  delete muon_py_;
  delete muon_pz_;
  delete muon_e_;
  delete muon_charge_;
  delete muon_iso_;
  delete muon_dxy_;
  delete muon_dz_;

  delete electron_pt_;
  delete electron_eta_;
  delete electron_phi_;
  delete electron_dxy_;
  delete electron_dz_;
  delete electron_px_;
  delete electron_py_;
  delete electron_pz_;
  delete electron_e_;
  delete electron_charge_;

  delete met_;
  delete met_x_;
  delete met_y_;
  delete met_phi_;
  delete num_bjets_ak8_;
  delete num_bjets_ak4_;
  delete num_jets_ak4_;

  delete run_;
  delete ev_;
  delete lumiblock_;
  delete crossingAngle_;
  delete betaStar_;
  delete instLumi_;
  delete nVertices_;
  delete pileupWeight_;
  delete mc_pu_trueinteractions_;
  delete mcWeight_;

  delete pfcand_nextracks_;
  delete pfcand_nextracks_noDRl_;
  delete recoMWhad_;
  delete recoMWlep_;
  delete dphiWW_;
  delete recoMWW_;
  delete recoRapidityWW_;
  delete recoWWpt_;
  delete genWWpt_;
  delete WLeptonicPt_;
  delete WLeptonicPhi_;
  delete WLeptonicEta_;

  delete jet_pt_;
  delete jet_px_;
  delete jet_py_;
  delete jet_pz_;
  delete jet_energy_;
  delete jet_phi_;
  delete jet_eta_;
  delete jet_mass_;
  delete jet_tau1_;
  delete jet_tau2_;
  delete jet_corrmass_;
  delete jet_vertexz_;
  delete jet_jer_res_;
  delete jet_jer_sf_;
  delete jet_jer_sfup_;
  delete jet_jer_sfdown_;
  delete pps_track_x_;
  delete pps_track_y_;
  delete pps_track_rpid_;
  
  delete proton_xi_;
  delete proton_thy_;
  delete proton_thx_;
  delete proton_t_;
  delete proton_ismultirp_;
  delete proton_rpid_;
  delete proton_arm_;
  delete proton_time_;
  delete proton_trackx1_;
  delete proton_tracky1_;
  delete proton_trackx2_;
  delete proton_tracky2_;
  delete proton_trackpixshift1_;
  delete proton_trackpixshift2_;
  delete proton_rpid1_;
  delete proton_rpid2_;
  
  delete gen_W_pt_;
  delete gen_W_charge_;
  delete gen_proton_px_;
  delete gen_proton_py_;
  delete gen_proton_pz_;
  delete gen_proton_xi_;
  delete gen_proton_energy_;
  delete gen_proton_t_;
  delete gen_jet_pt_;
  delete gen_jet_phi_;
  delete gen_jet_eta_;
  delete gen_jet_energy_;
  delete hlt_;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntupler);
