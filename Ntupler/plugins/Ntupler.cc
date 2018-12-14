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
  pps_token_(consumes<std::vector<CTPPSLocalTrackLite>>(edm::InputTag(("ctppsLocalTrackLiteProducer")))),
  gen_part_token_(consumes<reco::GenParticleCollection>(edm::InputTag("genParticles")))
{
  
  usesResource("TFileService");

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


   // Run and vertex multiplicity info
   *run_ = iEvent.id().run();
   *ev_ = iEvent.id().event();
   *lumiblock_ = iEvent.luminosityBlock();

      

   // Proton lite tracks
   edm::Handle<std::vector<CTPPSLocalTrackLite> > ppsTracks;
   iEvent.getByToken( pps_token_, ppsTracks );

   for ( const auto& trk : *ppsTracks ) 
     {
       const CTPPSDetId detid( trk.getRPId() );
       // transform the raw, 32-bit unsigned integer detId into the TOTEM "decimal" notation
       const unsigned short raw_id = 100*detid.arm()+10*detid.station()+detid.rp();

       (*pps_track_x_).push_back(trk.getX());
       (*pps_track_y_).push_back(trk.getY());
       (*pps_track_rpid_).push_back(raw_id);
     }

 
   

   // Fill GEN infor if running on MC
   
   edm::Handle<reco::GenParticleCollection> genP;
   //iEvent.getByLabel("prunedGenParticles",genP);
   iEvent.getByLabel("genParticles",genP);
   
   for (reco::GenParticleCollection::const_iterator mcIter=genP->begin(); mcIter != genP->end(); mcIter++ ) {
     
     if((mcIter->pdgId() == 2212) && (fabs(mcIter->pz()) > 3000) && (mcIter->status() == 1))
       {
	 double thexi = 1 - ((mcIter->energy())/(13000.0/2.0));
	 double thet = -(std::pow(mcIter->pt(), 2));
	 double thepz = mcIter->pz();
	 double thepx = mcIter->px();
	 double thepy = mcIter->py();
	 
	 (*gen_proton_xi_).push_back(thexi);
	     (*gen_proton_t_).push_back(thet);
	     (*gen_proton_pz_).push_back(thepz);
	     (*gen_proton_py_).push_back(thepy);
	     (*gen_proton_px_).push_back(thepx);
       }
   }





   tree_->Fill();


     



   (*pps_track_x_).clear();
   (*pps_track_y_).clear();
   (*pps_track_rpid_).clear();

   (*gen_proton_px_).clear();
   (*gen_proton_py_).clear();
   (*gen_proton_pz_).clear();
   (*gen_proton_xi_).clear();
   (*gen_proton_t_).clear();


}



// ------------ method called once each job just before starting event loop  ------------
void 
Ntupler::beginJob()
{

  edm::Service<TFileService> fs; 
  tree_=fs->make<TTree>("SlimmedNtuple","SlimmedNtuple");


  pps_track_x_ = new std::vector<float>;
  pps_track_y_ = new std::vector<float>;
  pps_track_rpid_ = new std::vector<int>;

  gen_proton_px_ = new std::vector<float>;
  gen_proton_py_ = new std::vector<float>;
  gen_proton_pz_ = new std::vector<float>;
  gen_proton_xi_ = new std::vector<float>;
  gen_proton_t_ = new std::vector<float>;

  ev_ = new long int;
  run_ = new int;
  lumiblock_ = new int;

  tree_->Branch("pps_track_x",&pps_track_x_);
  tree_->Branch("pps_track_y",&pps_track_y_);
  tree_->Branch("pps_track_rpid",&pps_track_rpid_);
  tree_->Branch("gen_proton_px",&gen_proton_px_);
  tree_->Branch("gen_proton_py",&gen_proton_py_);
  tree_->Branch("gen_proton_pz",&gen_proton_pz_);
  tree_->Branch("gen_proton_xi",&gen_proton_xi_);
  tree_->Branch("gen_proton_t",&gen_proton_t_);
  tree_->Branch("run",run_,"run/I");
  tree_->Branch("event",ev_,"event/L");
  tree_->Branch("lumiblock",lumiblock_,"lumiblock/I");


}

// ------------ method called once each jo bjust after ending the event loop  ------------
void 
Ntupler::endJob() 
{

  delete pps_track_x_;
  delete pps_track_y_;
  delete pps_track_rpid_;
  delete gen_proton_px_;
  delete gen_proton_py_;
  delete gen_proton_pz_;
  delete gen_proton_xi_;
  delete gen_proton_t_;

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
