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
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
  
  cout<<"I get to beginning of constructor"<<endl;
  //usesResource("TFileService");
  isMC = iConfig.getParameter<bool>("ismc");
  isPPS = iConfig.getParameter<bool>("ispps");
  channel = iConfig.getParameter<string>("channel");

  cout<<"channel is: "<<channel<<endl;
  consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
  consumes<std::vector<reco::Muon>>(edm::InputTag("muons"));
  consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));
  consumes<edm::DetSetVector<TotemRPLocalTrack>>(edm::InputTag("totemRPLocalTrackFitter"));
  consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>(edm::InputTag("ctppsPixelLocalTracks"));
  consumes<edm::DetSetVector<CTPPSDiamondLocalTrack> >(edm::InputTag("ctppsDiamondLocalTracks")); 
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
  consumes<reco::GsfElectronCollection >(edm::InputTag("gedGsfElectrons"));
  consumes<reco::ConversionCollection>(edm::InputTag("allConversions"));
  consumes<edm::View<pat::Jet>>(edm::InputTag("selectedPatJets"));
  consumes<edm::View<pat::Jet>>(edm::InputTag("tightPatJetsPFlow"));
  consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("addPileupInfo"));
  beamSpotToken_= consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  eleIdMapToken_=consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"));
  eleIdFullInfoMapToken_=consumes<edm::ValueMap<vid::CutFlowResult> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"));
  //eleIdMapToken_=consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"));
  //eleIdFullInfoMapToken_=consumes<edm::ValueMap<vid::CutFlowResult> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"));


  if(isMC){
    //LumiWeights = new edm::LumiReWeighting("MCPileup.root","MyDataPileupHistogramNotMu.root","h1","pileup");
    LumiWeights = new edm::LumiReWeighting("MCPileupHighStats.root","MyDataPileupHistogram0to75_MuonPhys.root","h_trueNumInter","pileup");
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

   using namespace edm;
   using namespace std;
   
   //This is not optimal if storing a lot of MC variables. In that case need to put this after trigger requirement. 
   //But need to be careful because h_trueNumInteractions in GetMC() needs to be before passTrigger
   //Right now it is ok because not saving any MC variables
   *pileupWeight_=1;
   if(isMC){
     GetMC(iEvent);
   }//end of looking at MC

   bool passTrigger=GetTrigger(iEvent,iSetup);
   //passTrigger=true;
   if(passTrigger){
     //cout<<"I pass trigger"<<endl;
     h_passTrigger->Fill(1);
     int numMuTight=0;
     int numETight = 0;

     if(isPPS){
       GetProtons(iEvent);
     }//end of if statement making sure that we want to look at these runs
     else{*ispps_=false;}

       
     //Get Vertices
     edm::Handle< std::vector<reco::Vertex> > vtxs;
     iEvent.getByLabel("offlinePrimaryVertices", vtxs);
     std::vector<reco::Vertex>::const_iterator vtxIt ;
     *vertex_nvtxs_ = vtxs.product()->size();
     for (vtxIt = vtxs->begin(); vtxIt != vtxs->end(); ++vtxIt) {
       (*allvertices_z_).push_back(vtxIt->position().z());
     }
     reco::VertexRef vtx(vtxs, 0);
     
     // get RECO tracks from the event
     edm::Handle<reco::TrackCollection> tks;
     iEvent.getByLabel("generalTracks", tks);
     //get the builder:
     edm::ESHandle<TransientTrackBuilder> theB;
     iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
     //do the conversion:
     vector<reco::TransientTrack> t_tks = (*theB).build(tks);
     std::vector<reco::TransientTrack> ttrkC_mu;
     std::vector<reco::TransientTrack> ttrkC_e_gsf;
     std::vector<reco::TransientTrack> ttrkC_e_ctf;
     std::vector<reco::TransientTrack> ttrkC;
     std::vector<uint> ttrkC_mu_it;
     std::vector<uint> ttrkC_e_ctf_it;

     GetMuons(iEvent,vtx,theB,ttrkC_mu,ttrkC_mu_it,t_tks,numMuTight);
     GetElectrons(iEvent,vtx,theB,ttrkC_e_gsf,ttrkC_e_ctf,ttrkC_e_ctf_it,t_tks,numETight);
     GetTracksPrimaryVertex(vtx,ttrkC_mu,ttrkC_e_ctf);
     
     TransientVertex myVertex;
     //If there is a good dilepton fit for this channel get lepton and track distances and track counting.
     if(FitLeptonVertex(myVertex,ttrkC,ttrkC_mu,ttrkC_e_gsf,channel)){
       
       GetMuonDistance(myVertex,ttrkC_mu);
       GetElectronDistance(myVertex,ttrkC_e_gsf);
       //Need to pass ttrkC_mu_it and ttrkC_e_ctf_it so not to count electron and muon ctf tracks
       GetTrackDistance(myVertex,t_tks,ttrkC_mu_it,ttrkC_e_ctf_it);
       //This counts electron and muon ctf tracks within 0.05 cm
     }//end of requiring valid vertex
     else{
       *fvertex_x_=999.;
       *fvertex_y_=999.;
       *fvertex_z_=999.;
       *fvertex_chi2ndof_=999.;
       *fvertex_nltracks_p5mm_=999.;
       *fvertex_ntracks_cms_p5mm_=999.;
       *fvertex_ntracks_ts_p5mm_=999.;
       *fvertex_ntracks_cms_p3mm_=999.;
       *fvertex_ntracks_ts_p3mm_=999.;
       *fvertex_closest_trk_cms_=999.;
       *fvertex_closest_trk_ts_=999.;
     }

     if(iEvent.id().event()==2974946501){
       cout<<"Num Mu: "<<numMuTight;
       cout<<"Num E: "<<numETight;
     }
     
     if((numMuTight+numETight)>1){
       //if(numMuTight==1&&numETight==1){
       GetJets(iEvent);
       *run_ = iEvent.id().run();
       *ev_ = iEvent.id().event();
       *lumiblock_ = iEvent.luminosityBlock();
       tree_->Fill();
     }
   }//end of looking at passing trigger
   
   
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
   (*electron_px_).clear();
   (*electron_py_).clear();
   (*electron_pz_).clear();
   (*electron_e_).clear();
   (*electron_charge_).clear();
   (*electron_passip_).clear();
   (*electron_dxy_).clear();
   (*electron_dz_).clear();

   (*jet_pt_).clear();
   (*jet_energy_).clear();
   (*jet_phi_).clear();
   (*jet_eta_).clear();
   
   (*allvertices_z_).clear();
   
   (*fvertex_tkdist_ts_p3mm_to_1p5mm_).clear();
   (*fvertex_tkdist_cms_p3mm_to_1p5mm_).clear();
   (*fvertex_tkpt_).clear();
   (*fvertex_tketa_).clear();
   (*muon_tkdist_).clear();
   (*muon_tkpt_).clear();
   (*electron_tkdist_).clear();
   (*electron_tkpt_).clear();


   if(isPPS){
     //(*rp_tracks_xraw_).clear();
     (*rp_tracks_y_).clear();
     (*rp_tracks_x_).clear();
     (*rp_tracks_xi_).clear();
     (*rp_tracks_xi_unc_).clear();
     (*rp_tracks_detId_).clear();
     (*rp_tracks_time_).clear();
   }

}

bool Ntupler::GetTrigger(const edm::Event& iEvent,const edm::EventSetup& iSetup)
{

   bool passTrigger=false;
   edm::Handle<edm::TriggerResults> hltResults;
   iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"),hltResults);
   const edm::TriggerNames & trigNames = iEvent.triggerNames(*hltResults);
   for(unsigned int i=0; i<trigNames.size();i++){
     int prescale_value=hltPrescaleProvider_.prescaleValue(iEvent, iSetup,trigNames.triggerName(i));
     //cout<<"trigNames.triggerName(i)"<<trigNames.triggerName(i)<<endl;
     //These triggers are for 2017 data
     if(channel=="mumu"){
       if((trigNames.triggerName(i)=="HLT_IsoMu27_v10"
	   ||trigNames.triggerName(i)=="HLT_IsoMu27_v8"
	   ||trigNames.triggerName(i)=="HLT_IsoMu27_v9"
	   ||trigNames.triggerName(i)=="HLT_IsoMu27_v11"
	   ||trigNames.triggerName(i)=="HLT_IsoMu27_v12"
	   ||trigNames.triggerName(i)=="HLT_IsoMu27_v13"
	   ||trigNames.triggerName(i)=="HLT_IsoMu27_v14")
	  //&&hltResults->accept(i)>0&&prescale_value==1){
	  &&hltResults->accept(i)>0){
	 //cout<<"Prescale value: "<<prescale_value<<endl;
	 passTrigger=true;
       }
     }//end of requirement of mumu channel

     //These are for 2017 data
     //Other triggers that are possible are: HLT Ele35 WPTight ,     HLT_Ele23_Ele12_CaloIdL_TrkIdL_IsoVL 
     if(channel=="ee"){
       if((trigNames.triggerName(i)=="HLT_DoubleEle33_CaloIdL_MW_v9"
	   ||trigNames.triggerName(i)=="HLT_DoubleEle33_CaloIdL_MW_v10"
	   ||trigNames.triggerName(i)=="HLT_DoubleEle33_CaloIdL_MW_v11"
	   ||trigNames.triggerName(i)=="HLT_DoubleEle33_CaloIdL_MW_v12"
	   ||trigNames.triggerName(i)=="HLT_DoubleEle33_CaloIdL_MW_v13"
	   ||trigNames.triggerName(i)=="HLT_DoubleEle33_CaloIdL_MW_v14"
	   ||trigNames.triggerName(i)=="HLT_DoubleEle33_CaloIdL_MW_v15")
	  &&hltResults->accept(i)>0&&prescale_value==1){
	 passTrigger=true;
       }
     }//end of requirement of ee channel
     
     //These are for 2017 data
     if(channel=="mue"){
       if((trigNames.triggerName(i)=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"||
	   trigNames.triggerName(i)=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2"||
	   trigNames.triggerName(i)=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3"||
	   trigNames.triggerName(i)=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4"||
	   trigNames.triggerName(i)=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5")
	  &&hltResults->accept(i)>0&&prescale_value==1){
       passTrigger=true;
       }
     }
   }//end of looping over triggers

   return passTrigger;

}

void
Ntupler::GetProtons(const edm::Event& iEvent)
{



	 edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > aodTracks;
	 iEvent.getByLabel("totemRPLocalTrackFitter",aodTracks);
	 
	 // produce collection of lite tracks (in future this will be done in miniAOD)        
	 
	 //TrackDataCollection liteTracks;
	 for (const auto &ds : *aodTracks)
	   {
	     for (const auto &tr : ds)
	       {
		 //liteTracks[rpId] = tr;
		 //cout<<"I get to the RP tracks, output X0: "<<tr.getX0()<<endl;
		 (*rp_tracks_x_).push_back(tr.getX0());   
		 (*rp_tracks_y_).push_back(tr.getY0()); 
		 //(*rp_tracks_tx_).push_back(tr.getTx());   
		 //(*rp_tracks_ty_).push_back(tr.getTy()); 
		 //liteTracks[rpId]=tr;
		 //cout<<"ds.detID(): "<<ds.detId()<<endl;
		 (*rp_tracks_detId_).push_back(ds.detId());   
	       }
	   }

	 
	 edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > pixLocalTracks;
	 iEvent.getByLabel("ctppsPixelLocalTracks",pixLocalTracks);
	 for (const auto &ds : *pixLocalTracks)
	   {
	     for (const auto &tr : ds)
	       {
		 //liteTracks[rpId] = tr;
		 //cout<<"I get to the RP tracks, output X0: "<<tr.getX0()<<endl;
		 (*rp_tracks_x_).push_back(tr.getX0());   
		 (*rp_tracks_y_).push_back(tr.getY0()); 
		 //(*rp_tracks_tx_).push_back(tr.getTx());   
		 //(*rp_tracks_ty_).push_back(tr.getTy()); 
		 //liteTracks[rpId]=tr;
		 //cout<<"ds.detID(): "<<ds.detId()<<endl;
		 (*rp_tracks_detId_).push_back(ds.detId());   
	       }
	   }
	 
	 edm::Handle< edm::DetSetVector<CTPPSDiamondLocalTrack> > diamondLocalTracks;
	 iEvent.getByLabel("ctppsDiamondLocalTracks",diamondLocalTracks);
	 for (const auto &ds : *diamondLocalTracks)
	   {
	     for (const auto &tr : ds)
	       {
		 //liteTracks[rpId] = tr;
		 //cout<<"I get to the RP tracks, output X0: "<<tr.getX0()<<endl;
		 (*rp_tracks_x_).push_back(tr.getX0());   
		 (*rp_tracks_y_).push_back(tr.getY0()); 
		 //(*rp_tracks_tx_).push_back(tr.getTx());   
		 //(*rp_tracks_ty_).push_back(tr.getTy()); 
		 //liteTracks[rpId]=tr;
		 //cout<<"ds.detID(): "<<ds.detId()<<endl;
		 (*rp_tracks_detId_).push_back(ds.detId());   
		 (*rp_tracks_time_).push_back(tr.getT());   
	       }
	   }

	 




}

void 
Ntupler::GetMC(const edm::Event& iEvent)
{
  
  //Need to fill this so can know total number of events in MC (since I don't keep every event  
  float trueInteractions=0;
  edm::InputTag PileupSrc_("addPileupInfo");
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(PileupSrc_, PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  //cout<<"True num interactions is: "<<PupInfo->begin()->getTrueNumInteractions()<<endl;
  trueInteractions=PupInfo->begin()->getTrueNumInteractions();
  h_trueNumInteractions->Fill(trueInteractions);
  *pileupWeight_ = LumiWeights->weight( trueInteractions );
     
  //cout<<" I get into MC"<<endl;
  /*
    Handle<reco::GenParticleCollection> genP;
    iEvent.getByLabel("genParticles",genP);
    for (reco::GenParticleCollection::const_iterator mcIter=genP->begin(); mcIter != genP->end(); mcIter++ ) {
    //cout<<"MC id is: "<<mcIter->pdgId()<<endl;
    //cout<<"MC status is: "<<mcIter->status()<<endl;
    //cout<<"Pz, energy, pt is: "<<mcIter->pz()<<", "<<mcIter->energy()<<", "<<mcIter->pt()<<endl;
    int n = mcIter->numberOfDaughters();
    for(int j = 0; j < n; ++ j) {
    const reco::Candidate * d = mcIter->daughter( j );
    int dauId = d->pdgId();
    //cout<<"Daughter pdg Id: "<<dauId<<endl;
    }
    
    //	     cout<<"Pt, eta, phi is: "<<mcIter->pt()<<", "<<mcIter->eta()<<", "<<mcIter->phi()<<endl;
    //if ( (fabs(mcIter->pdgId())==11|| fabs(mcIter->pdgId())==12 || fabs(mcIter->pdgId())==13 || fabs(mcIter->pdgId())==14 || fabs(mcIter->pdgId())==15 || fabs(mcIter->pdgId())==16 ) && mcIter->status() == 3 ){
    //if ( (fabs(mcIter->pdgId())==11|| fabs(mcIter->pdgId())==12 || fabs(mcIter->pdgId())==13 || fabs(mcIter->pdgId())==14 || fabs(mcIter->pdgId())==15 || fabs(mcIter->pdgId())==16 )){
    if ( (fabs(mcIter->pdgId())==11 || fabs(mcIter->pdgId())==13 )){
    h_lepton_pt->Fill(mcIter->pt());
    //cout<<", MC id is: "<<mcIter->pdgId()<<endl;
    //cout<<"Pt, eta, phi is: "<<mcIter->pt()<<", "<<mcIter->eta()<<", "<<mcIter->phi()<<endl;
    }
    }//end of looking at GEN
  */

}

void 
Ntupler::GetMuons(const edm::Event& iEvent,reco::VertexRef vtx,edm::ESHandle<TransientTrackBuilder>& theB,std::vector<reco::TransientTrack>& ttrkC_mu,std::vector<uint>& ttrkC_mu_it,vector<reco::TransientTrack> t_tks,int& numMuTight)
{
     edm::Handle< std::vector<reco::Muon> > muonHandle;
     iEvent.getByLabel("muons",muonHandle);
     std::vector<reco::Muon>::const_iterator MuonIt ;
     
     TLorentzVector mu1,mu2;

     for (MuonIt = muonHandle->begin(); MuonIt != muonHandle->end(); ++MuonIt) {
       //cout<<"Muon pt is: "<<MuonIt->pt()<<endl;
       
       //bool tightId = muon::isTightMuon(*MuonIt,*vtx);
       
       bool tightId_noIP=false;
       
       bool muID = muon::isGoodMuon(*MuonIt,muon::GlobalMuonPromptTight) && (MuonIt->numberOfMatchedStations() > 1);
       bool hits=false;
       if(MuonIt->isGlobalMuon()){
	 	 hits = MuonIt->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && MuonIt->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;}
       //bool ip = fabs(MuonIt->muonBestTrack()->dxy(vtx->position())) < 0.2 && fabs(MuonIt->muonBestTrack()->dz(vtx->position())) < 0.5;
       if(MuonIt->isPFMuon() && MuonIt->isGlobalMuon() && muID && hits ){ tightId_noIP=true; }
       
       //cout<<"Tight Muon Id is: "<<tightId<<endl;
       //if(tightId&&MuonIt->pt()>30&&fabs(MuonIt->eta())<2.4){
       if(tightId_noIP&&MuonIt->pt()>30&&fabs(MuonIt->eta())<2.4){


	 double iso = (MuonIt->pfIsolationR04().sumChargedHadronPt + max(0., MuonIt->pfIsolationR04().sumNeutralHadronEt + MuonIt->pfIsolationR04().sumPhotonEt - 0.5*MuonIt->pfIsolationR04().sumPUPt))/MuonIt->pt();
	 //cout<<"Muon Iso is: "<<iso<<endl;
	 //cout<<"Muon dz: "<<MuonIt->muonBestTrack()->dz(vtx->position())<<endl;
	 //cout<<"Muon dxy: "<<MuonIt->muonBestTrack()->dxy(vtx->position())<<endl;
	 (*muon_px_).push_back(MuonIt->px());
	 (*muon_py_).push_back(MuonIt->py());
	 (*muon_pz_).push_back(MuonIt->pz());
	 (*muon_e_).push_back(MuonIt->energy());
	 (*muon_charge_).push_back(MuonIt->charge());
	 (*muon_pt_).push_back(MuonIt->pt());
	 (*muon_eta_).push_back(MuonIt->eta());
	 (*muon_phi_).push_back(MuonIt->phi());
	 (*muon_iso_).push_back(iso);
	 (*muon_dxy_).push_back(fabs(MuonIt->muonBestTrack()->dxy(vtx->position())));
	 (*muon_dz_).push_back(fabs(MuonIt->muonBestTrack()->dz(vtx->position())));
	 //cout<<"Muon Pt: "<<MuonIt->pt()<<endl;
	 //cout<<"Muon charge: "<<MuonIt->charge()<<endl;
	 //cout<<"Vertex track size: "<<vtx->tracksSize()<<endl;
	 //cout<<"Get right before track muon"<<endl;
	 //ttrkC_mu.push_back((*theB).build(*MuonIt->innerTrack()));
	 const reco::TransientTrack b = (*theB).build(*MuonIt->innerTrack());
	 ttrkC_mu.push_back(b);
	 //cout<<"Get right after track muon"<<endl;
	 //ttrkC_mu_it.push_back(i);
	 //for(const auto at : t_tks){
	 for (uint i=0; i < t_tks.size();i++){
	   if(fabs(MuonIt->innerTrack()->pt()-t_tks[i].track().pt())<0.001&&fabs(MuonIt->innerTrack()->eta()-t_tks[i].track().eta())<0.001&&fabs(MuonIt->innerTrack()->phi()-t_tks[i].track().phi())<0.001){
	     //ttrkC_mu.push_back(t_tks[i]);
	     ttrkC_mu_it.push_back(i);
	     }
	   }
	   //reco::TrackRef mutrk = MuonIt->innerTrack();
	   
	 numMuTight++;
	 if(numMuTight==1){	   mu1.SetPx(MuonIt->px());mu1.SetPy(MuonIt->py());mu1.SetPz(MuonIt->pz());mu1.SetE(MuonIt->energy());	 }
	 if(numMuTight==2){	   mu2.SetPx(MuonIt->px());mu2.SetPy(MuonIt->py());mu2.SetPz(MuonIt->pz());mu2.SetE(MuonIt->energy());	 }
	 if(numMuTight>2){cout<<"There are more than 3 tight muons in the event"<<endl;}
       }//end of looking at tightId
     }//end of looking at muons
     if(numMuTight==2){
       TLorentzVector mumu = mu1+mu2;
       //*mumu_mass_=mumu.M();
       //*mumu_rapidity_=mumu.Rapidity();
       //cout<<"Invariant mass: "<<mumu.M()<<endl;
	 //cout<<"Rapidity: "<<mumu.Rapidity()<<endl;
     }
     



}


void 
Ntupler::GetElectrons(const edm::Event& iEvent,reco::VertexRef vtx,edm::ESHandle<TransientTrackBuilder>& theB,std::vector<reco::TransientTrack>& ttrkC_e_gsf,std::vector<reco::TransientTrack>& ttrkC_e_ctf,std::vector<uint>& ttrkC_e_ctf_it,vector<reco::TransientTrack> t_tks,int& numETight)
{

     //Electron information
     //https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_8.0.3/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc       
     edm::Handle<edm::ValueMap<bool> > ele_id_decisions;
     iEvent.getByToken(eleIdMapToken_ ,ele_id_decisions);
     // Full cut flow info for one of the working points:
     edm::Handle<edm::ValueMap<vid::CutFlowResult> > ele_id_cutflow_data;
     iEvent.getByToken(eleIdFullInfoMapToken_,ele_id_cutflow_data);
     
     
     edm::Handle<reco::BeamSpot> theBeamSpot;
     iEvent.getByToken(beamSpotToken_,theBeamSpot);
     
     // conversions                                                             
     edm::Handle<reco::ConversionCollection> conversions_h;
     iEvent.getByLabel("allConversions", conversions_h);
     
     //Loop over electrons in event
     edm::Handle<reco::GsfElectronCollection> els_h;
     iEvent.getByLabel("gedGsfElectrons", els_h);
     unsigned int n = els_h->size();

     for(unsigned int i = 0; i < n; ++i) {
       reco::GsfElectronRef ele(els_h, i);
       //reco::GsfTrackRef theTrack = el->gsfTrack();
       //dz_.push_back( theTrack->dz( firstGoodVertex->position() ) );
       //dxy_.push_back( theTrack->dxy( firstGoodVertex->position() ) );
       bool passConvVeto = !ConversionTools::hasMatchedConversion(*ele,conversions_h,theBeamSpot->position());
       bool isPassEleId = (*ele_id_decisions)[ele];
              
       
       //See if electron passes
       
       if(isPassEleId&&passConvVeto&&ele->pt()>30&&fabs(ele->superCluster()->eta())<2.4){
	 //reco::VertexRef vtx(vtxs, 0);
	 /*
	 cout<<"The privary vertex position: "<<vtx->position()<<endl;
	 cout<<"The beamspot position: "<<theBeamSpot->position()<<endl;
	 cout<<"electron d0: "<<ele->gsfTrack()->d0()<<endl;
	 cout<<"electron d0,bs: "<<ele->gsfTrack()->dxy(theBeamSpot->position())<<endl;
	 cout<<"electron d0,pv: "<<ele->gsfTrack()->dxy(vtx->position())<<endl;
	 cout<<"electron dz: "<<ele->gsfTrack()->dz()<<endl;
	 cout<<"electron dz,bs: "<<ele->gsfTrack()->dz(theBeamSpot->position())<<endl;
	 cout<<"electron dz,pv: "<<ele->gsfTrack()->dz(vtx->position())<<endl;
	 cout<<"electron vertex z: "<<ele->vertex().z();
	 */
	 bool passIPcuts=false;
	 if(fabs(ele->superCluster()->eta())<=1.479){
	   if(fabs(ele->gsfTrack()->dz(vtx->position()))<0.10&&fabs(ele->gsfTrack()->dxy(vtx->position()))<0.05){passIPcuts=true;}
	 }
	 if(fabs(ele->superCluster()->eta())>1.479){
	   if(fabs(ele->gsfTrack()->dz(vtx->position()))<0.20&&fabs(ele->gsfTrack()->dxy(vtx->position()))<0.10){passIPcuts=true;}
	 }
	 //cout<<"Does it pass IP cuts: "<<passIPcuts<<endl;
	 numETight++;
	 (*electron_pt_).push_back(ele->pt());
	 (*electron_passip_).push_back(passIPcuts);
	 (*electron_dxy_).push_back(fabs(ele->gsfTrack()->dxy(vtx->position())));
	 (*electron_dz_).push_back(fabs(ele->gsfTrack()->dz(vtx->position())));
	 //cout<<"Get Here 0:"<<endl;
	 //cout<<"Electron pt, eta, phi: "<<ele->pt()<<", "<<ele->eta()<<", "<<ele->phi()<<endl;
	 //cout<<"Electron charge: "<<ele->charge()<<endl;
	 //cout<<"Electron px, py, pz: "<<ele->px()<<", "<<ele->py()<<", "<<ele->pz()<<endl;
	 //cout<<"Electron pt from px, py: "<<sqrt(ele->px()*ele->px()+ele->py()*ele->py())<<endl;
	 //cout<<"Electron energy: "<<ele->energy()<<endl;
	 //if(ele->closestCtfTrackRef().isNonnull()){
	 //cout<<"Electron Ctf track pt, eta, phi: "<<ele->closestCtfTrackRef()->pt()<<", "<<ele->closestCtfTrackRef()->eta()<<", "<<ele->closestCtfTrackRef()->phi()<<endl;
	 //cout<<"Electron Ctf track px, py, pz: "<<ele->closestCtfTrackRef()->px()<<", "<<ele->closestCtfTrackRef()->py()<<", "<<ele->closestCtfTrackRef()->pz()<<endl;}
	 //cout<<"Electron GSF track pt, eta, phi: "<<ele->gsfTrack()->pt()<<", "<<ele->gsfTrack()->eta()<<", "<<ele->gsfTrack()->phi()<<endl;
	 (*electron_eta_).push_back(ele->superCluster()->eta());
	 (*electron_phi_).push_back(ele->superCluster()->phi());
	 (*electron_px_).push_back(ele->px());
	 (*electron_py_).push_back(ele->py());
	 (*electron_pz_).push_back(ele->pz());
	 (*electron_e_).push_back(ele->energy());
	 (*electron_charge_).push_back(ele->charge());
	 //cout<<"electron gsf pt"<<ele->gsfTrack()->pt()<<endl;
	 const reco::TransientTrack b = (*theB).build(*ele->gsfTrack());
	 //cout<<"Valid?: "<<b.isValid()<<endl;
	 ttrkC_e_gsf.push_back(b);
	 //ttrkC_e.push_back((*theB).build(*ele->gsfTrack()));
	 //ttrkC_e.push_back((*theB).build(ele->track()));
	 //cout<<(*theB).build(*ele->gsfTrack()).track().pt()<<endl;
	 //cout<<"Get right after track electron"<<endl;
	 if(ele->closestCtfTrackRef().isNonnull()){
	   //for(const auto at : t_tks){
	   for (uint i=0; i < t_tks.size();i++){
	     if(fabs(ele->closestCtfTrackRef()->pt()-t_tks[i].track().pt())<0.001&&fabs(ele->closestCtfTrackRef()->eta()-t_tks[i].track().eta())<0.001&&fabs(ele->closestCtfTrackRef()->phi()-t_tks[i].track().phi())<0.001){
	       //cout<<"This is the correct electron track, pt: "<<at.track().pt()<<endl;
	       //		   ttrkC_e.push_back(at);
	       //ttrkC_e.push_back(t_tks[i]);
	       ttrkC_e_ctf_it.push_back(i);
	       ttrkC_e_ctf.push_back(t_tks[i]);
	     }
	   }//end of looping over tracks to get track matching to electron
	 }//making sure closest Ctf Track is non-null
	 
     }//requirement that electron pass id and conv veto
       vid::CutFlowResult fullCutFlowData = (*ele_id_cutflow_data)[ele];
       bool verbose_electron=false;
       if(verbose_electron){
	 printf("\nDEBUG CutFlow, full info for cand with pt=%f:\n", ele->pt());
	 //printCutFlowResult(fullCutFlowData);
	 printf("    CutFlow name= %s    decision is %d\n", 
		fullCutFlowData.cutFlowName().c_str(),
		(int) fullCutFlowData.cutFlowPassed());
	 int ncuts = fullCutFlowData.cutFlowSize();
	 printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
	 for(int icut = 0; icut<ncuts; icut++){
	   printf("  %2d      %50s    %d        %f          %d\n", icut,
		  fullCutFlowData.getNameAtIndex(icut).c_str(),
		  (int)fullCutFlowData.isCutMasked(icut),
		  fullCutFlowData.getValueCutUpon(icut),
		  (int)fullCutFlowData.getCutResultByIndex(icut));
	 }
       }//end of looking at electrons cutflow
     }//end of loop over electrons



}

void
Ntupler::GetTracksPrimaryVertex(reco::VertexRef vtx,std::vector<reco::TransientTrack> ttrkC_mu,std::vector<reco::TransientTrack> ttrkC_e_ctf)
{

  *vertex_ntracks_ = vtx->tracksSize();
  //cout<<"vertex_ntracks_"<<vtx->tracksSize()<<endl;
  *vertex_x_ = vtx->position().x();
  *vertex_y_ = vtx->position().y();
  *vertex_z_ = vtx->position().z();
  

  int pass_muon_assoc=0;
  int pass_electron_assoc=0;
  //Look at primary vertex and loop to see if muon and electron match it
  for(reco::Vertex::trackRef_iterator vertex_Tracks = vtx->tracks_begin();vertex_Tracks<vtx->tracks_end(); vertex_Tracks++){
    for(uint it = 0; it<ttrkC_mu.size();it++){
      if( fabs((*vertex_Tracks)->pt()-ttrkC_mu[it].track().pt())<0.001 && 
	  fabs((*vertex_Tracks)->eta()-ttrkC_mu[it].track().eta())<0.001 && 
	  fabs((*vertex_Tracks)->phi()-ttrkC_mu[it].track().phi())<0.001){
	     pass_muon_assoc++;
      }
    }
    for(uint it = 0; it<ttrkC_e_ctf.size();it++){
      if (fabs((*vertex_Tracks)->pt()-ttrkC_e_ctf[it].track().pt())<0.001 && 
	  fabs((*vertex_Tracks)->eta()-ttrkC_e_ctf[it].track().eta())<0.001 && 
	  fabs((*vertex_Tracks)->phi()-ttrkC_e_ctf[it].track().phi())<0.001){
	pass_electron_assoc++;
      }
    }//end of for loop over electrons
  }
  //cout<<"pass_muon_assoc: "<<pass_muon_assoc<<endl;
  //cout<<"pass_electron_assoc: "<<pass_electron_assoc<<endl;
  if(channel=="mue"){ if(!(pass_muon_assoc==1&&pass_electron_assoc==1)){*vertex_ntracks_=1001;}      }
  if(channel=="mumu"){ if(pass_muon_assoc!=2){*vertex_ntracks_=1001;}       }
  if(channel=="ee"){ if(pass_electron_assoc!=2){*vertex_ntracks_=1001;}       }
  
}

void
Ntupler::GetMuonDistance(TransientVertex myVertex,std::vector<reco::TransientTrack> ttrkC_mu)
{

  //calculate muon distance
  for(uint att=0;att<ttrkC_mu.size();att++){
    //TrajectoryStateClosestToPoint tS_muon=ttrkC_mu[att].trajectoryStateClosestToPoint(myVertex.position());
    //if(tS_muon.isValid()){

      float closest_pos = sqrt( pow(myVertex.position().x()-ttrkC_mu[att].track().vertex().x(),2)+pow(myVertex.position().y()-ttrkC_mu[att].track().vertex().y(),2)+pow(myVertex.position().z()-ttrkC_mu[att].track().vertex().z(),2));  
      //float closest_pos = sqrt( pow(myVertex.position().x()-tS_muon.position().x(),2)+pow(myVertex.position().y()-tS_muon.position().y(),2)+pow(myVertex.position().z()-tS_muon.position().z(),2));
      h_mu_closest->Fill(closest_pos);
      (*muon_tkdist_).push_back(closest_pos);
      (*muon_tkpt_).push_back(ttrkC_mu[att].track().pt());
      if(myVertex.normalisedChiSquared()<10){h_mu_closest_chi2_10->Fill(closest_pos);}
      h_mu_chi2_vs_closest->Fill(closest_pos,myVertex.normalisedChiSquared());
      

    //cout<<"Closest muon distance: "<<closest_pos<<", Run: "<<iEvent.id().run()<<", Event: "<<iEvent.id().event()<<endl;
      //closest_pos = sqrt( pow(myVertex.position().x()-ttrkC_mu[att].track().vertex().x(),2)+pow(myVertex.position().y()-ttrkC_mu[att].track().vertex().y(),2)+pow(myVertex.position().z()-ttrkC_mu[att].track().vertex().z(),2));
      //cout<<"Closest muon distance: "<<closest_pos<<", Run: "<<iEvent.id().run()<<", Event: "<<iEvent.id().event()<<endl;
      //cout<<"muon track dist:"<<closest_pos<<endl;
      //cout<<"muon track pt:"<<ttrkC_mu[att].track().pt()<<endl;
      //if(closest_pos<0.05){ num_close_tracks++;}	   
      //    }
  }//end of loop over muons
  
}

void
Ntupler::GetElectronDistance(TransientVertex myVertex,std::vector<reco::TransientTrack> ttrkC_e_gsf)
{
  
  //calculate electron distance
  for(uint att=0;att<ttrkC_e_gsf.size();att++){
    //cout<<"Get in electrons"<<endl;
    float closest_pos = sqrt( pow(myVertex.position().x()-ttrkC_e_gsf[att].track().vertex().x(),2)+pow(myVertex.position().y()-ttrkC_e_gsf[att].track().vertex().y(),2)+pow(myVertex.position().z()-ttrkC_e_gsf[att].track().vertex().z(),2));
    h_e_closest->Fill(closest_pos);
    (*electron_tkdist_).push_back(closest_pos);
    (*electron_tkpt_).push_back(ttrkC_e_gsf[att].track().pt());
    if(myVertex.normalisedChiSquared()<10){h_e_closest_chi2_10->Fill(closest_pos);}
    h_e_chi2_vs_closest->Fill(closest_pos,myVertex.normalisedChiSquared());
    /*
      if(closest_pos > 0.5){
      cout<<"GSF Track vertex"<<ttrkC_e_gsf[att].track().vertex()<<endl;
      cout<<"My vertex: "<<myVertex.position()<<endl;
      cout<<"Closest electon distance: "<<closest_pos<<", Run: "<<iEvent.id().run()<<", Event: "<<iEvent.id().event()<<endl;
      }
    */
    //cout<<"electron track dist:"<<closest_pos<<endl;
    //TrajectoryStateClosestToPoint tS_e=ttrkC_e_gsf[att].trajectoryStateClosestToPoint(myVertex.position());
    //cout<<"Get after traj."<<endl;
    //if(tS_e.isValid()){
    //cout<<"tS is valid."<<endl;
    //closest_pos = sqrt( pow(myVertex.position().x()-tS_e.position().x(),2)+pow(myVertex.position().y()-tS_e.position().y(),2)+pow(myVertex.position().z()-tS_e.position().z(),2));
  }//end of electron distance
  

}

void
Ntupler::GetTrackDistance(TransientVertex myVertex,std::vector<reco::TransientTrack> t_tks,std::vector<uint> ttrkC_mu_it,std::vector<uint> ttrkC_e_ctf_it)
{
  
  //calculate track distance
  float closest_track_cms=1000;
  float closest_track_TrajState=1000;
  int num_lepton_tracks_cms_p5mm=0;
  int num_tracks_cms_p5mm=0;
  int num_tracks_cms_p3mm=0;
  int num_tracks_ts_p5mm=0;
  int num_tracks_ts_p3mm=0;
  for (uint i=0; i < t_tks.size();i++){

    float track_distance_cms = sqrt( pow(myVertex.position().x()-t_tks[i].track().vertex().x(),2)+pow(myVertex.position().y()-t_tks[i].track().vertex().y(),2)+pow(myVertex.position().z()-t_tks[i].track().vertex().z(),2));

    bool isEorMu=false;
    for(uint att=0;att<ttrkC_mu_it.size();att++){
      if(ttrkC_mu_it[att]==i){
	isEorMu=true;
	if(track_distance_cms<0.05){ num_lepton_tracks_cms_p5mm++;}
      }
    }
    for(uint att=0;att<ttrkC_e_ctf_it.size();att++){
      if(ttrkC_e_ctf_it[att]==i){
	isEorMu=true;
	if(track_distance_cms<0.05){ num_lepton_tracks_cms_p5mm++;}
      }
    }

    if(!isEorMu){
      if(track_distance_cms<closest_track_cms){
	closest_track_cms=track_distance_cms;
      }
      if(track_distance_cms<0.05){  num_tracks_cms_p5mm++; }
      if(track_distance_cms<0.03){  num_tracks_cms_p3mm++; }
      if(track_distance_cms<0.15&&track_distance_cms>0.03){
	(*fvertex_tkdist_cms_p3mm_to_1p5mm_).push_back(track_distance_cms);
      }
    }//end of requirement excluding e or mu tracks
    TrajectoryStateClosestToPoint tS=t_tks[i].trajectoryStateClosestToPoint(myVertex.position());
    //cout<<"Closest positivon on track: "<<tS.position().x()<<", "<<tS.position().y()<<", "<<tS.position().z()<<endl;
    
    //believe this is all in cm
    if(tS.isValid()){
      float track_distance_TrajState = sqrt( pow(myVertex.position().x()-tS.position().x(),2)+pow(myVertex.position().y()-tS.position().y(),2)+pow(myVertex.position().z()-tS.position().z(),2));
      //if(t_tks[i].track().pt()>6){cout<<"Closest position: "<<closest_pos<<endl;}

      if(!isEorMu){
	if(track_distance_TrajState<closest_track_TrajState){
	  closest_track_TrajState=track_distance_TrajState;
	}
	if(track_distance_TrajState<0.05){  num_tracks_ts_p5mm++; }
	if(track_distance_TrajState<=0.03){  num_tracks_ts_p3mm++; }
	if(track_distance_TrajState<0.15&&track_distance_TrajState>0.03){
	  (*fvertex_tkdist_ts_p3mm_to_1p5mm_).push_back(track_distance_TrajState);
	  //(*fvertex_tkpt_).push_back(t_tks[i].track().pt());
	  //(*fvertex_tketa_).push_back(t_tks[i].track().eta());
	}//fill ntuple with tracks within 0.2 cm
      }//end of requirement excluding e or mu tracks
    }//end of making sure Trajectory state is valid
    else{cout<<"TrajectoryStateClosestToPoint is not valid"<<endl;}
  }//end of looping over tracks

  *fvertex_closest_trk_cms_=closest_track_cms;
  *fvertex_closest_trk_ts_=closest_track_TrajState;
  *fvertex_ntracks_cms_p5mm_=num_tracks_cms_p5mm;
  *fvertex_ntracks_ts_p5mm_=num_tracks_ts_p5mm;
  *fvertex_ntracks_cms_p3mm_=num_tracks_cms_p3mm;
  *fvertex_ntracks_ts_p3mm_=num_tracks_ts_p3mm;
  *fvertex_nltracks_p5mm_=num_lepton_tracks_cms_p5mm;

}

bool
Ntupler::FitLeptonVertex(TransientVertex& myVertex,std::vector<reco::TransientTrack> ttrkC,std::vector<reco::TransientTrack> ttrkC_mu,std::vector<reco::TransientTrack> ttrkC_e_gsf,string channel)
{

  bool vertexExistAndValid=false;
  bool fitVertex = false;

  if(channel=="mue"){
    if(ttrkC_mu.size()==1&&ttrkC_e_gsf.size()==1){fitVertex=true;}
  }
  if(channel=="mumu"){
    if(ttrkC_mu.size()==2){fitVertex=true;}
  }
  if(channel=="ee"){
    if(ttrkC_e_gsf.size()==2){fitVertex=true;}
  }
  
  if(fitVertex){
    //Do fit to two lepton tracks
    //AdaptiveVertexFitter fitter;
    KalmanVertexFitter fitter;
    if(channel=="mue"){
      ttrkC.push_back(ttrkC_mu[0]);
      ttrkC.push_back(ttrkC_e_gsf[0]);
      myVertex = fitter.vertex(ttrkC);
    }
    if(channel=="mumu"){
      myVertex = fitter.vertex(ttrkC_mu);
    }
    if(channel=="ee"){
      myVertex = fitter.vertex(ttrkC_e_gsf);
    }
    //cout<<"Get Before valid vertex"<<endl;
    
    if(myVertex.isValid()){
      vertexExistAndValid=true;
      *fvertex_x_=myVertex.position().x();
      *fvertex_y_=myVertex.position().y();
      *fvertex_z_=myVertex.position().z();
      *fvertex_chi2ndof_=myVertex.normalisedChiSquared();
    }
    else{
      cout<<"Fitted vertex is not valid"<<endl;
    }
    
  }//end of requirement of fitting vertex to two tracks
  return vertexExistAndValid;
  
}


void
Ntupler::GetJets(const edm::Event& iEvent)
{

       //cout<<"Muon phi, eta: "<<(*muon_phi_)[0]<<", "<<(*muon_eta_)[0]<<endl;
       //cout<<"Electron phi, eta: "<<(*electron_phi_)[0]<<", "<<(*electron_eta_)[0]<<endl;
       //GetJets(iEvent);
       edm::Handle<edm::View<pat::Jet> > jetColl; // PAT      
       //iEvent.getByToken(jetToken_, jetColl);
       //iEvent.getByLabel("selectedPatJets", jetColl);
       iEvent.getByLabel("tightPatJetsPFlow", jetColl);
       int numJets_id=0;
       int numJets_id_25=0;
       int numJets_id_20=0;
       for (unsigned int i=0; i<jetColl->size(); i++) {
	 const edm::Ptr<pat::Jet> jet = jetColl->ptrAt(i);
	 //cout<<"chargedHadronEnergyFraction(): "<<jet->chargedHadronEnergyFraction()<<endl;
	 //cout<<"chargedHadronEnergyFraction(), 2nd method: "<<jet->chargedHadronEnergy()/jet->correctedP4(0).E()<<endl;
	 
	 if(jet->pt()>25){
	   bool isLepton=isJetLepton(jet->eta(),jet->phi());
	   if(!isLepton){
	     numJets_id++;  
	     (*jet_pt_).push_back(jet->pt());
	     (*jet_energy_).push_back(jet->energy());
	     (*jet_phi_).push_back(jet->phi());
	     (*jet_eta_).push_back(jet->eta());
	   }
	   //cout<<"Jet energy, phi, eta: "<<jet->energy()<<", "<<jet->phi()<<", "<<jet->eta()<<endl;
	   //cout<<"Jet corrected energy: "<<jet->energy()<<endl;
	   //cout<<"Uncorrected energy: "<<jet->correctedP4(0).E()<<endl;
	   //cout<<"Jet corrected pt: "<<jet->pt()<<endl;
	   //cout<<"Uncorrected pt: "<<jet->correctedP4(0).Pt()<<endl;
	   std::vector< std::string > jec_levels = jet->availableJECLevels();
	   int size = jec_levels.size();
	   //for(int i=0;i<size;i++){
	   //  cout<<jec_levels[i]<<endl;
	   //}
	 }//end of looping over jet pt > 25
	 if(jet->pt()>25){numJets_id_25++;       }
	 if(jet->pt()>20){numJets_id_20++;       }
	 //cout<<"Uncorrected energy: "<<jet->correctedP4(0).E()<<endl;
	 //cout<<"Uncorrected chargedHadronEnergyFraction(): "<<jet->correctedJet(0).chargedHadronEnergyFraction()<<endl;
       }

       /*       
       edm::Handle<edm::View<pat::Jet> > jetSColl; // PAT      
       iEvent.getByLabel("selectedPatJets", jetSColl);
       int numJets=0;
       for (unsigned int i=0; i<jetSColl->size(); i++) {
	 const edm::Ptr<pat::Jet> jet = jetSColl->ptrAt(i);
	 //cout<<"Jet energy, pt, eta: "<<jet->energy()<<", "<<jet->pt()<<", "<<jet->eta()<<endl;
	 //cout<<"chargedHadronEnergyFraction(): "<<jet->chargedHadronEnergyFraction()<<endl;
	 //cout<<"chargedHadronEnergyFraction(), 2nd method: "<<jet->chargedHadronEnergy()/jet->correctedP4(0).E()<<endl;
	 if(jet->pt()>30)
	   {numJets++;}
	 //cout<<"Uncorrected energy: "<<jet->correctedP4(0).E()<<endl;
	 //cout<<"Uncorrected chargedHadronEnergyFraction(): "<<jet->correctedJet(0).chargedHadronEnergyFraction()<<endl;
       }
       //cout<<"Num jets with no id, with id: "<<numJets<<", "<<numJets_id<<endl;
       */
       h_jets_30->Fill(numJets_id);
       h_jets_25->Fill(numJets_id_25);
       h_jets_20->Fill(numJets_id_20);


}


bool
Ntupler::isJetLepton(double jet_eta, double jet_phi)
{
  for(uint i=0;i<(*muon_eta_).size();i++){
    double eta=(*muon_eta_)[i];
    double phi=(*muon_phi_)[i];
    double deltaR=sqrt((eta-jet_eta)*(eta-jet_eta)+(phi-jet_phi)*(phi-jet_phi));
    if(deltaR<0.5){
      return true;
    }
  }

  for(uint i=0;i<(*electron_eta_).size();i++){
    double eta=(*electron_eta_)[i];
    double phi=(*electron_phi_)[i];
    double deltaR=sqrt((eta-jet_eta)*(eta-jet_eta)+(phi-jet_phi)*(phi-jet_phi));
    if(deltaR<0.5){
      return true;
    }
  }



  return false;

}


void
Ntupler::endRun(edm::Run const & iRun, edm::EventSetup const& iSetup) {}


void 
Ntupler::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

  //cout<<"I get to beginning of beginRun"<<endl;
  bool changed(true);
  hltPrescaleProvider_.init(iRun,iSetup,"HLT",changed);
  //  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
  if (hltConfig_.init(iRun,iSetup,"HLT",changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
      // The HLT config has actually changed wrt the previous Run, hence rebook your
      // histograms or do anything else dependent on the revised HLT config
     
    }
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    //    coutr") << " HLT config extraction failure with process name " << processName_<<endl;
    //    LogError("MyAnalyzer") << " HLT config extraction failure with process name " << "HLT";
    std::cout << " HLT config extraction failure with process name " << "HLT"<<std::endl;
    // In this case, all access methods will return empty values!
  }

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
  electron_px_ = new std::vector<float>;
  electron_py_ = new std::vector<float>;
  electron_pz_ = new std::vector<float>;
  electron_e_ = new std::vector<float>;
  electron_charge_ = new std::vector<float>;
  electron_passip_ = new std::vector<bool>;
  electron_dxy_ = new std::vector<float>;
  electron_dz_ = new std::vector<float>;
  allvertices_z_ = new std::vector<float>;

  jet_pt_ = new std::vector<float>;
  jet_energy_ = new std::vector<float>;
  jet_phi_ = new std::vector<float>;
  jet_eta_ = new std::vector<float>;

  vertex_ntracks_ = new int;
  vertex_x_ = new float;
  vertex_y_ = new float;
  vertex_z_ = new float;
  vertex_nvtxs_ = new int;

  fvertex_x_ = new float;
  fvertex_y_ = new float;
  fvertex_z_ = new float;
  fvertex_chi2ndof_ = new float;
  fvertex_nltracks_p5mm_ = new int;
  fvertex_ntracks_cms_p5mm_ = new int;
  fvertex_ntracks_ts_p5mm_ = new int;
  fvertex_ntracks_cms_p3mm_ = new int;
  fvertex_ntracks_ts_p3mm_ = new int;
  fvertex_closest_trk_cms_ = new float;
  fvertex_closest_trk_ts_ = new float;
  fvertex_tkdist_ts_p3mm_to_1p5mm_ = new std::vector<float>;
  fvertex_tkdist_cms_p3mm_to_1p5mm_ = new std::vector<float>;
  fvertex_tkpt_ = new std::vector<float>;
  fvertex_tketa_ = new std::vector<float>;
  muon_tkdist_ = new std::vector<float>;
  muon_tkpt_ = new std::vector<float>;
  electron_tkdist_ = new std::vector<float>;
  electron_tkpt_ = new std::vector<float>;

  if(isPPS){
    //rp_tracks_xraw_ = new std::vector<float>;
    rp_tracks_y_ = new std::vector<float>;
    rp_tracks_x_ = new std::vector<float>;
    rp_tracks_xi_ = new std::vector<float>;
    rp_tracks_xi_unc_ = new std::vector<float>;
    rp_tracks_detId_ = new std::vector<float>;
    rp_tracks_time_ = new std::vector<float>;
    //mumu_mass_ = new float;
    //mumu_rapidity_ = new float;
  }

  ev_ = new long int;
  run_ = new int;
  lumiblock_ = new int;
  ispps_ = new bool;

  //Tnpv_ = new float;
  pileupWeight_ = new float;


  
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
  tree_->Branch("electron_px",&electron_px_);
  tree_->Branch("electron_py",&electron_py_);
  tree_->Branch("electron_pz",&electron_pz_);
  tree_->Branch("electron_e",&electron_e_);
  tree_->Branch("electron_charge",&electron_charge_);
  tree_->Branch("electron_passip",&electron_passip_);
  tree_->Branch("electron_dxy",&electron_dxy_);
  tree_->Branch("electron_dz",&electron_dz_);
  tree_->Branch("allvertices_z",&allvertices_z_);
  tree_->Branch("jet_pt",&jet_pt_);
  tree_->Branch("jet_energy",&jet_energy_);
  tree_->Branch("jet_phi",&jet_phi_);
  tree_->Branch("jet_eta",&jet_eta_);
  tree_->Branch("vertex_ntracks",vertex_ntracks_,"vertex_ntracks/I");
  tree_->Branch("vertex_x",vertex_x_,"vertex_x/f");
  tree_->Branch("vertex_y",vertex_y_,"vertex_y/f");
  tree_->Branch("vertex_z",vertex_z_,"vertex_z/f");
  tree_->Branch("vertex_nvtxs",vertex_nvtxs_,"vertex_nvtxs/I");
  tree_->Branch("fvertex_x",fvertex_x_,"fvertex_x/f");
  tree_->Branch("fvertex_y",fvertex_y_,"fvertex_y/f");
  tree_->Branch("fvertex_z",fvertex_z_,"fvertex_z/f");
  tree_->Branch("fvertex_chi2ndof",fvertex_chi2ndof_,"fvertex_chi2ndof/f");
  tree_->Branch("fvertex_nltracks_p5mm",fvertex_nltracks_p5mm_,"fvertex_nltracks_p5mm/I");
  tree_->Branch("fvertex_ntracks_cms_p5mm",fvertex_ntracks_cms_p5mm_,"fvertex_ntracks_cms_p5mm/I");
  tree_->Branch("fvertex_ntracks_ts_p5mm",fvertex_ntracks_ts_p5mm_,"fvertex_ntracks_ts_p5mm/I");
  tree_->Branch("fvertex_ntracks_cms_p3mm",fvertex_ntracks_cms_p3mm_,"fvertex_ntracks_cms_p3mm/I");
  tree_->Branch("fvertex_ntracks_ts_p3mm",fvertex_ntracks_ts_p3mm_,"fvertex_ntracks_ts_p3mm/I");
  tree_->Branch("fvertex_closest_trk_cms",fvertex_closest_trk_cms_,"fvertex_closest_trk_cms/f");
  tree_->Branch("fvertex_closest_trk_ts",fvertex_closest_trk_ts_,"fvertex_closest_trk_ts/f");
  tree_->Branch("fvertex_tkdist_ts_p3mm_to_1p5mm",&fvertex_tkdist_ts_p3mm_to_1p5mm_);
  tree_->Branch("fvertex_tkdist_cms_p3mm_to_1p5mm",&fvertex_tkdist_cms_p3mm_to_1p5mm_);
  tree_->Branch("fvertex_tkpt",&fvertex_tkpt_);
  tree_->Branch("fvertex_tketa",&fvertex_tketa_);
  tree_->Branch("muon_tkdist",&muon_tkdist_);
  tree_->Branch("muon_tkpt",&muon_tkpt_);
  tree_->Branch("electron_tkdist",&electron_tkdist_);
  tree_->Branch("electron_tkpt",&electron_tkpt_);

  if(isPPS){
    //tree_->Branch("rp_tracks_xraw",&rp_tracks_xraw_);
    tree_->Branch("rp_tracks_y",&rp_tracks_y_);
    tree_->Branch("rp_tracks_x",&rp_tracks_x_);
    tree_->Branch("rp_tracks_xi",&rp_tracks_xi_);
    tree_->Branch("rp_tracks_xi_unc",&rp_tracks_xi_unc_);
    tree_->Branch("rp_tracks_detId",&rp_tracks_detId_);
    tree_->Branch("rp_tracks_time",&rp_tracks_time_);
    //tree_->Branch("mumu_mass",mumu_mass_,"mumu_mass/f");
    //tree_->Branch("mumu_rapidity",mumu_rapidity_,"mumu_rapidity/f");
  }
  tree_->Branch("run",run_,"run/I");
  tree_->Branch("event",ev_,"event/L");
  tree_->Branch("lumiblock",lumiblock_,"lumiblock/I");
  tree_->Branch("ispps",ispps_,"ispps/O");

  //tree_->Branch("Tnpv",Tnpv_,"Tnpv/f");
  tree_->Branch("pileupWeight",pileupWeight_,"pileupWeight/f");


  h_trueNumInteractions = fs->make<TH1F>("h_trueNumInteractions" , "PU" , 200 , 0 , 100 );
  h_passTrigger = fs->make<TH1F>("h_passTrigger" , "PU" , 1 , 0.5 , 1.5 );
  h_mu_closest = fs->make<TH1F>("h_mu_closest" , ";Track Distance (cm);" , 500 , 0 , 10 );
  h_mu_closest_chi2_10 = fs->make<TH1F>("h_mu_closest_chi2_10" , "For #chi^2 < 10;Track Distance (cm);" , 500 , 0 , 10 );
  h_mu_chi2_vs_closest = fs->make<TH2F>("h_mu_chi2_vs_closest" , ";Track Distance (cm);#chi^{2}" , 500 , 0 , 10 , 200, 0, 50);

  h_e_closest = fs->make<TH1F>("h_e_closest" , ";Track Distance (cm);" , 500 , 0 , 10 );
  h_e_closest_chi2_10 = fs->make<TH1F>("h_e_closest_chi2_10" , "For #chi^2 < 10;Track Distance (cm);" , 500 , 0 , 10 );
  h_e_chi2_vs_closest = fs->make<TH2F>("h_e_chi2_vs_closest" , ";Track Distance (cm);#chi^{2}" , 500 , 0 , 10 , 200, 0, 50);
  h_lepton_pt = fs->make<TH1F>("h_lepton_pt" , ";Lepton p_{T};" , 150, 0, 300);
  h_jets_30 = fs->make<TH1F>("h_jets_30" , ";Num jets;" , 15, -0.5, 14.5);
  h_jets_25 = fs->make<TH1F>("h_jets_25" , ";Num jets;" , 15, -0.5, 14.5);
  h_jets_20 = fs->make<TH1F>("h_jets_20" , ";Num jets;" , 15, -0.5, 14.5);



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

  delete electron_px_;
  delete electron_pt_;
  delete electron_eta_;
  delete electron_phi_;
  delete electron_py_;
  delete electron_pz_;
  delete electron_e_;
  delete electron_charge_;
  delete electron_passip_;
  delete electron_dxy_;
  delete electron_dz_;

  delete jet_pt_;
  delete jet_energy_;
  delete jet_phi_;
  delete jet_eta_;

  delete allvertices_z_;

  delete vertex_ntracks_;
  delete vertex_x_;
  delete vertex_y_;
  delete vertex_z_;
  delete vertex_nvtxs_;

  delete fvertex_x_;
  delete fvertex_y_;
  delete fvertex_z_;
  delete fvertex_chi2ndof_;
  delete fvertex_nltracks_p5mm_;
  delete fvertex_ntracks_cms_p5mm_;
  delete fvertex_ntracks_ts_p5mm_;
  delete fvertex_ntracks_cms_p3mm_;
  delete fvertex_ntracks_ts_p3mm_;
  delete fvertex_closest_trk_cms_;
  delete fvertex_closest_trk_ts_;
  delete fvertex_tkdist_ts_p3mm_to_1p5mm_;
  delete fvertex_tkdist_cms_p3mm_to_1p5mm_;
  delete fvertex_tkpt_;
  delete fvertex_tketa_;
  delete muon_tkdist_;
  delete muon_tkpt_;
  delete electron_tkdist_;
  delete electron_tkpt_;

  if(isPPS){
    //delete rp_tracks_xraw_;
    delete rp_tracks_y_;
    delete rp_tracks_x_;
    delete rp_tracks_xi_;
    delete rp_tracks_xi_unc_;
    delete rp_tracks_detId_;
    delete rp_tracks_time_;
    //delete mumu_mass_;
    //delete mumu_rapidity_;
  }

  delete run_;
  delete ev_;
  delete lumiblock_;
  //delete Tnpv_;
  delete pileupWeight_;
  delete ispps_;


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
