// -*- C++ -*-
//
// Package:    SlimmedNtuple/Ntupler
// Class:      Ntupler
// 
/**\class Ntupler Ntupler.cc SlimmedNtuple/Ntupler/plugins/Ntupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Finn O'Neill Rebassoo
//         Created:  Thu, 10 Nov 2016 01:41:43 GMT
//
//


// system include files
#include <memory>
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <iostream>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TLorentzVector.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"

//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//For Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//TOTEM reco
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/CTPPSDigi/interface/TotemRPDigi.h"
#include "DataFormats/CTPPSDigi/interface/TotemVFATStatus.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"

//For Conversions
//#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Conversion.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//For electrons
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"

//For electron id
//#include "DataFormats/Common/interface/ValueMap.h"
//#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "shared_track.h"
#include "shared_alignment.h"
#include "shared_reconstruction.h"
#include "shared_fill_info.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include <FWCore/Framework/interface/ESHandle.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Ntupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Ntupler(const edm::ParameterSet&);
      ~Ntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  // ----------member data ---------------------------
  
  edm::FileInPath fp;
  edm::FileInPath fp2;
  string fp0;
  string fp1;
  AlignmentResultsCollection alignment;
  TTree * tree_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  // ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapToken_;
  // One example of full information about the cut flow
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleIdFullInfoMapToken_;
  //HLTConfigProvider hltConfig_;

  std::vector<float> * muon_pt_;
  std::vector<float> * muon_eta_;
  std::vector<float> * muon_px_;
  std::vector<float> * muon_py_;
  std::vector<float> * muon_pz_;
  std::vector<float> * muon_e_;
  std::vector<float> * muon_charge_;

  uint * vertex_ntracks_;
  float * vertex_x_;
  float * vertex_y_;
  float * vertex_z_;
  uint * vertex_nvtxs_;

  float * fvertex_x_;
  float * fvertex_y_;
  float * fvertex_z_;
  float * fvertex_chi2ndof_;
  uint * fvertex_ntracks_;
  std::vector<float> * fvertex_tkdist_;
  std::vector<float> * fvertex_tkpt_;
  std::vector<float> * fvertex_tketa_;


  std::vector<float> * rp_tracks_xraw_;
  std::vector<float> * rp_tracks_y_;
  std::vector<float> * rp_tracks_x_;
  std::vector<float> * rp_tracks_xi_;
  std::vector<float> * rp_tracks_detId_;

  uint * run_;
  uint * ev_;
  uint * lumiblock_;

  float * mumu_mass_;
  float * mumu_rapidity_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ntupler::Ntupler(const edm::ParameterSet& iConfig)

{
  

  cout<<"I get to beginning of constructor"<<endl;
  //now do what ever initialization is needed
  usesResource("TFileService");
  consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
  consumes<std::vector<reco::Muon>>(edm::InputTag("muons"));
  consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));
  consumes< edm::DetSetVector<TotemRPLocalTrack> >(edm::InputTag("totemRPLocalTrackFitter"));
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  //consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));

  //consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMap"))
  consumes<reco::GsfElectronCollection >(edm::InputTag("gedGsfElectrons"));
  consumes<reco::ConversionCollection>(edm::InputTag("allConversions"));

  beamSpotToken_= consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  //eleIdMapToken_=consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"));
eleIdMapToken_=consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"));
  //eleIdFullInfoMapToken_=consumes<edm::ValueMap<vid::CutFlowResult> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"));
 eleIdFullInfoMapToken_=consumes<edm::ValueMap<vid::CutFlowResult> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"));

  fp0 = iConfig.getParameter<string>("particleFile");
  fp1 = iConfig.getParameter<string>("particleFile2");
  //fp = iConfig.getParameter<edm::FileInPath>("particleFile");
  //fp2 = iConfig.getParameter<edm::FileInPath>("particleFile2");
  //cout<<"file path alignment: "<<fp.fullPath()<<endl;
  //cout<<"file path optics: "<<fp.fullPath()<<endl;
  //cout<<"Is local?: "<<fp.isLocal()<<endl;

  InitReconstruction(fp1);
  InitFillInfoCollection();

  if (alignment.Load(fp0.c_str()) != 0)
    {  
      printf("ERROR: can't load alignment data.\n");
    }


  muon_pt_ = new std::vector<float>;
  muon_eta_ = new std::vector<float>;
  muon_px_ = new std::vector<float>;
  muon_py_ = new std::vector<float>;
  muon_pz_ = new std::vector<float>;
  muon_e_ = new std::vector<float>;
  muon_charge_ = new std::vector<float>;

  vertex_ntracks_ = new uint;
  vertex_x_ = new float;
  vertex_y_ = new float;
  vertex_z_ = new float;
  vertex_nvtxs_ = new uint;

  fvertex_x_ = new float;
  fvertex_y_ = new float;
  fvertex_z_ = new float;
  fvertex_chi2ndof_ = new float;
  fvertex_ntracks_ = new uint;
  fvertex_tkdist_ = new std::vector<float>;
  fvertex_tkpt_ = new std::vector<float>;
  fvertex_tketa_ = new std::vector<float>;

  rp_tracks_xraw_ = new std::vector<float>;
  rp_tracks_y_ = new std::vector<float>;
  rp_tracks_x_ = new std::vector<float>;
  rp_tracks_xi_ = new std::vector<float>;
  rp_tracks_detId_ = new std::vector<float>;

  ev_ = new uint;
  run_ = new uint;
  lumiblock_ = new uint;

  mumu_mass_ = new float;
  mumu_rapidity_ = new float;

  edm::Service<TFileService> fs; 
  tree_=fs->make<TTree>("SlimmedNtuple","SlimmedNtuple");


  tree_->Branch("muon_pt",&muon_pt_);
  tree_->Branch("muon_eta",&muon_eta_);
  tree_->Branch("muon_px",&muon_px_);
  tree_->Branch("muon_py",&muon_py_);
  tree_->Branch("muon_pz",&muon_pz_);
  tree_->Branch("muon_e",&muon_e_);
  tree_->Branch("muon_charge",&muon_charge_);
  tree_->Branch("vertex_ntracks",vertex_ntracks_,"vertex_ntracks/i");
  tree_->Branch("vertex_x",vertex_x_,"vertex_x/f");
  tree_->Branch("vertex_y",vertex_y_,"vertex_y/f");
  tree_->Branch("vertex_z",vertex_z_,"vertex_z/f");
  tree_->Branch("vertex_nvtxs",vertex_nvtxs_,"vertex_nvtxs/i");
  tree_->Branch("fvertex_x",fvertex_x_,"fvertex_x/i");
  tree_->Branch("fvertex_y",fvertex_y_,"fvertex_y/i");
  tree_->Branch("fvertex_z",fvertex_z_,"fvertex_z/i");
  tree_->Branch("fvertex_chi2ndof",fvertex_chi2ndof_,"fvertex_chi2ndof/f");
  tree_->Branch("fvertex_ntracks",fvertex_ntracks_,"fvertex_ntracks/i");
  tree_->Branch("fvertex_tkdist",&fvertex_tkdist_);
  tree_->Branch("fvertex_tkpt",&fvertex_tkpt_);
  tree_->Branch("fvertex_tketa",&fvertex_tketa_);
  tree_->Branch("rp_tracks_xraw",&rp_tracks_xraw_);
  tree_->Branch("rp_tracks_y",&rp_tracks_y_);
  tree_->Branch("rp_tracks_x",&rp_tracks_x_);
  tree_->Branch("rp_tracks_xi",&rp_tracks_xi_);
  tree_->Branch("rp_tracks_detId",&rp_tracks_detId_);
  
  tree_->Branch("run",run_,"run/i");
  tree_->Branch("event",ev_,"event/i");
  tree_->Branch("lumiblock",lumiblock_,"lumiblock/i");

  tree_->Branch("mumu_mass",mumu_mass_,"mumu_mass/f");
  tree_->Branch("mumu_rapidity",mumu_rapidity_,"mumu_rapidity/f");

}


Ntupler::~Ntupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)


  delete muon_px_;
  delete muon_pt_;
  delete muon_eta_;
  delete muon_py_;
  delete muon_pz_;
  delete muon_e_;
  delete muon_charge_;

  delete vertex_ntracks_;
  delete vertex_x_;
  delete vertex_y_;
  delete vertex_z_;
  delete vertex_nvtxs_;

  delete fvertex_x_;
  delete fvertex_y_;
  delete fvertex_z_;
  delete fvertex_chi2ndof_;
  delete fvertex_tkdist_;
  delete fvertex_tkpt_;
  delete fvertex_tketa_;

  delete rp_tracks_xraw_;
  delete rp_tracks_y_;
  delete rp_tracks_x_;
  delete rp_tracks_xi_;
  delete rp_tracks_detId_;
  

  delete run_;
  delete ev_;
  delete lumiblock_;

  delete mumu_mass_;
  delete mumu_rapidity_;

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
   //if(iEvent.id().event()==883153054||iEvent.id().event()==814088551){
     //cout<<"iEvent.id().run()"<<iEvent.id().run()<<endl;
     //cout<<"iEvent.id().luminosityBlock"<<iEvent.luminosityBlock()<<endl;
     //cout<<"iEvent.id().event()"<<iEvent.id().event()<<endl;

     bool passTrigger=false;
     edm::Handle<TriggerResults> hltResults;
     iEvent.getByLabel(InputTag("TriggerResults","","HLT"),hltResults);
     const TriggerNames & trigNames = iEvent.triggerNames(*hltResults);
     for(unsigned int i=0; i<trigNames.size();i++){
       //cout<<"Trigger_name: "<<trigNames.triggerName(i)<<endl;
       //cout<<"Trigger decision: "<<hltResults->accept(i)<<endl;
       //cout<<"Prescale: "<<hltConfig_.prescaleValue(iEvent,iSetup,trigNames.triggerName(i))<<endl;
       if(trigNames.triggerName(i)=="HLT_DoubleMu38NoFiltersNoVtx_v3"&&hltResults->accept(i)>0){
	 passTrigger=true;
       }
     }
     passTrigger=true;
     if(passTrigger){
       //cout<<"I pass trigger"<<endl;
       /*
       //////////////////Start of reconstruction of RP si strips////////////////////////////////
       map<unsigned int, bool> tr; 
       tr[2] = false;
       tr[3] = false;
       tr[102] = false;
       tr[103] = false;
       Handle< edm::DetSetVector<TotemRPLocalTrack> > tracks;     
       iEvent.getByLabel("totemRPLocalTrackFitter",tracks);
                                                   
       TrackDataCollection trackData_raw;
       for (const auto &ds : *tracks)
	 {
	   const auto &rpId = ds.detId();
	   for (const auto &t : ds)
	     {
	       if (rpId == 3 || rpId == 2 || rpId == 102 || rpId == 103){
		 trackData_raw[rpId] = t;
		 //cout<<"I get raw track"<<endl;
		 //cout<<"x: "<<trackData_raw[rpId].x<<endl;
		 //cout<<"y: "<<trackData_raw[rpId].y<<endl;
		 (*rp_tracks_xraw_).push_back(trackData_raw[rpId].x);
	       }
	     }
	 }

       cout<<"I get here 3:"<<endl;                                                   
       const auto &fillInfo = fillInfoCollection.FindByRun(iEvent.id().run());
       cout<<"I get here 4:"<<endl;                                                   
       const auto alignment_it = alignment.find(fillInfo.alignmentTag);
       if (alignment_it == alignment.end())
	 {
	   printf("ERROR: no alignment for tag '%s'.\n", fillInfo.alignmentTag.c_str());
	   //return 1;
	 }
       //       else{
       // printf("ERROR: alignment for tag '%s'.\n", fillInfo.alignmentTag.c_str());
       //}
       TrackDataCollection trackData_al = alignment_it->second.Apply(trackData_raw);

       // split track collection per arm                                                                                                     
       TrackDataCollection trackData_L, trackData_R;
       for (const auto &p : trackData_al)
	 {
	   //int arm = p.first / 100;
	   TSpline3 *s_x_to_xi = m_s_x_to_xi[p.first];
	   
	   (*rp_tracks_x_).push_back(p.second.x);
	   (*rp_tracks_y_).push_back(p.second.y);
	   (*rp_tracks_detId_).push_back(p.first);
	   (*rp_tracks_xi_).push_back(s_x_to_xi->Eval(p.second.x*1E-3));
	 }
       //////////////////End of reconstruction of RP si strips////////////////////////////////
       */

       edm::Handle< std::vector<reco::Vertex> > vtxs;
       iEvent.getByLabel("offlinePrimaryVertices", vtxs);
       std::vector<reco::Vertex>::const_iterator vtxIt ;
       //cout<<"Number of vertices: "<<vtxs.product()->size();
       *vertex_nvtxs_ = vtxs.product()->size();
       //for (vtxIt = vtxs->begin(); vtxIt != vtxs->end(); ++vtxIt) {
       // cout<<"Vertex track size: "<<vtxIt->tracksSize()<<endl;
       //}

       edm::Handle< std::vector<reco::Muon> > muonHandle;
       iEvent.getByLabel("muons",muonHandle);
       std::vector<reco::Muon>::const_iterator MuonIt ;

       // get RECO tracks from the event
       edm::Handle<reco::TrackCollection> tks;
       iEvent.getByLabel("generalTracks", tks);

       //get the builder:
       edm::ESHandle<TransientTrackBuilder> theB;
       iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
       //do the conversion:
       vector<reco::TransientTrack> t_tks = (*theB).build(tks);
       std::vector<reco::TransientTrack>::const_iterator ttrk_It;
       //t_tks.setBeamSpot(beamSpot)
       std::vector<reco::TransientTrack> ttrkC;

       TLorentzVector mu1,mu2;
       int numMuTight=0;
       reco::VertexRef vtx(vtxs, 0);
       *vertex_ntracks_ = vtx->tracksSize();
       *vertex_x_ = vtx->position().x();
       *vertex_y_ = vtx->position().y();
       *vertex_z_ = vtx->position().z();
       for (MuonIt = muonHandle->begin(); MuonIt != muonHandle->end(); ++MuonIt) {
	 bool tightId = muon::isTightMuon(*MuonIt,*vtx);
	 if(tightId&&MuonIt->pt()>20&&fabs(MuonIt->eta())<2.4){
	   (*muon_px_).push_back(MuonIt->px());
	   (*muon_py_).push_back(MuonIt->py());
	   (*muon_pz_).push_back(MuonIt->pz());
	   (*muon_e_).push_back(MuonIt->energy());
	   (*muon_charge_).push_back(MuonIt->charge());
	   (*muon_pt_).push_back(MuonIt->pt());
	   (*muon_eta_).push_back(MuonIt->eta());
	   //cout<<"Pt: "<<MuonIt->pt()<<endl;
	   //cout<<"Vertex track size: "<<vtx->tracksSize()<<endl;
	   for(const auto at : t_tks){
	     if(fabs(MuonIt->pt()-at.track().pt())<0.001&&fabs(MuonIt->eta()-at.track().eta())<0.001&&fabs(MuonIt->phi()-at.track().phi())<0.001){
	       //cout<<"This is the correct track, pt: "<<MuonIt->pt()<<endl;
	       ttrkC.push_back(at);
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
	 *mumu_mass_=mumu.M();
	 *mumu_rapidity_=mumu.Rapidity();
	 //cout<<"Invariant mass: "<<mumu.M()<<endl;
	 //cout<<"Rapidity: "<<mumu.Rapidity()<<endl;
       }

       //Electron information
       //https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_8.0.3/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc       
       edm::Handle<edm::ValueMap<bool> > ele_id_decisions;
       iEvent.getByToken(eleIdMapToken_ ,ele_id_decisions);
       // Full cut flow info for one of the working points:
       edm::Handle<edm::ValueMap<vid::CutFlowResult> > ele_id_cutflow_data;
       iEvent.getByToken(eleIdFullInfoMapToken_,ele_id_cutflow_data);
       
       
       // beam spot                                                                                     
       //edm::Handle<reco::BeamSpot> beamspot_h;
       //iEvent.getByLabel(beamSpotInputTag_, beamspot_h);                                           
       //iEvent.getByLabel("offlineBeamSpot", beamspot_h);
       //const reco::BeamSpot &beamSpot = *(beamspot_h.product());

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

	 bool passConvVeto = !ConversionTools::hasMatchedConversion(*ele,conversions_h,theBeamSpot->position());

	 bool isPassEleId = (*ele_id_decisions)[ele];

	 vid::CutFlowResult fullCutFlowData = (*ele_id_cutflow_data)[ele];
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

       }
      

       if(ttrkC.size()>1){
	 AdaptiveVertexFitter fitter;
	 TransientVertex myVertex = fitter.vertex(ttrkC);
	 if(myVertex.isValid()){
	   *fvertex_x_=myVertex.position().x();
	   *fvertex_y_=myVertex.position().y();
	   *fvertex_z_=myVertex.position().z();
	   *fvertex_chi2ndof_=myVertex.normalisedChiSquared();
	   //cout<<"Position: "<<myVertex.position().x()<<", "<<myVertex.position().y()<<", "<<myVertex.position().z()<<endl;
	   //cout<<"Ndof: "<<myVertex.degreesOfFreedom()<<endl;
	   //cout<<"Normalized ChiSquared: "<<myVertex.normalisedChiSquared()<<endl;
	   //cout<<"ChiSquared: "<<myVertex.totalChiSquared()<<endl;
	   
	   uint num_close_tracks=-1;
	   //for (ttrk_It=t_tks->begin();ttrk_It != t_tks->end(); ++ttrk_It){
	   for (uint i=0; i < t_tks.size();i++){
	     //cout<<"Track pt: "<<t_tks[i].track().pt()<<endl;
	     //cout<<"Track eta: "<<t_tks[i].track().eta()<<endl;
	     TrajectoryStateClosestToPoint tS=t_tks[i].trajectoryStateClosestToPoint(myVertex.position());
	     //cout<<"Closest position on track: "<<tS.position().x()<<", "<<tS.position().y()<<", "<<tS.position().z()<<endl;
	     //believe this is all in cm
	     if(tS.isValid()){
	       float closest_pos = sqrt( pow(myVertex.position().x()-tS.position().x(),2)+pow(myVertex.position().y()-tS.position().y(),2)+pow(myVertex.position().z()-tS.position().z(),2));
	       //cout<<"Closest position: "<<closest_pos<<endl;
	       if(closest_pos<1){
		 (*fvertex_tkdist_).push_back(closest_pos);
		 (*fvertex_tkpt_).push_back(t_tks[i].track().pt());
		 (*fvertex_tketa_).push_back(t_tks[i].track().eta());
	       }//fill ntuple with tracks within 1 m
	       if(closest_pos<0.1){
		 num_close_tracks++;
	       }//end of counting tracks within 1 mm
	     }//end of making sure Trajectory state is valid
	     else{cout<<"TrajectoryStateClosestToPoint is not valid"<<endl;}
	   }//end of looping over tracks
	   *fvertex_ntracks_=num_close_tracks+1;
	 }//end of requiring valid vertex
	 else{cout<<"Fitted vertex is not valid"<<endl;
	   *fvertex_x_=-999.;
	   *fvertex_y_=-999.;
	   *fvertex_z_=-999.;
	   *fvertex_chi2ndof_=-999.;
	   cout<<"Number tracks at dimuon vertex: "<<*vertex_ntracks_<<endl;
	 }
       }//end of requirement of two tight muon tracks
     
     //}//end of looking at one specific event

   
       *run_ = iEvent.id().run();
       *ev_ = iEvent.id().event();
       *lumiblock_ = iEvent.luminosityBlock();
       
       tree_->Fill();

     }//end of looking at passing trigger

     (*muon_pt_).clear();
     (*muon_eta_).clear();
     (*muon_px_).clear();
     (*muon_py_).clear();
     (*muon_pz_).clear();
     (*muon_e_).clear();
     (*muon_charge_).clear();
     (*fvertex_tkdist_).clear();
     (*fvertex_tkpt_).clear();
     (*fvertex_tketa_).clear();
     (*rp_tracks_xraw_).clear();
     (*rp_tracks_y_).clear();
     (*rp_tracks_x_).clear();
     (*rp_tracks_xi_).clear();
     (*rp_tracks_detId_).clear();



}


// ------------ method called once each job just before starting event loop  ------------
void 
Ntupler::beginJob()
{
}

// ------------ method called once each jo bjust after ending the event loop  ------------
void 
Ntupler::endJob() 
{
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
