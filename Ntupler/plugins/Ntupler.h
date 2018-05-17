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
#include "TH2F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"


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
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"


//#include "shared_track.h"
//#include "shared_alignment.h"
//#include "shared_reconstruction.h"
//#include "shared_fill_info.h"

#include "track_lite.h"
#include "alignment.h"
#include "fill_info.h"
#include "proton_reconstruction.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include <FWCore/Framework/interface/ESHandle.h>

//
// class declaration
//

class Ntupler : public edm::EDAnalyzer {
   public:
      explicit Ntupler(const edm::ParameterSet&);
      ~Ntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void endRun(edm::Run const &, edm::EventSetup const&) override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  bool GetTrigger(const edm::Event&,const edm::EventSetup&);
  void GetProtons(const edm::Event&);
  void GetMC(const edm::Event&);
  void GetMuons(const edm::Event&,reco::VertexRef,edm::ESHandle<TransientTrackBuilder>&,std::vector<reco::TransientTrack>&,std::vector<uint>&,std::vector<reco::TransientTrack>,int&);
  void GetElectrons(const edm::Event&,reco::VertexRef,edm::ESHandle<TransientTrackBuilder>&,std::vector<reco::TransientTrack>&,std::vector<reco::TransientTrack>&,std::vector<uint>&,std::vector<reco::TransientTrack>,int&);
  void GetTracksPrimaryVertex(reco::VertexRef,std::vector<reco::TransientTrack>,std::vector<reco::TransientTrack>);
  void GetMuonDistance(TransientVertex,std::vector<reco::TransientTrack>);
  void GetElectronDistance(TransientVertex,std::vector<reco::TransientTrack>);
  void GetTrackDistance(TransientVertex,std::vector<reco::TransientTrack>,std::vector<uint>,std::vector<uint>,uint&);
  bool FitLeptonVertex(TransientVertex&,std::vector<reco::TransientTrack>,std::vector<reco::TransientTrack>,std::vector<reco::TransientTrack>,string);

  // ----------member data ---------------------------
  
  string fp0;
  string fp1;
  AlignmentResultsCollection alignmentCollection;
  AlignmentResults *alignments;
  unsigned int prev_run;
  bool prev_pps;
  TTree * tree_;
  bool isMC;
  bool isPPS;
  string channel;
  edm::LumiReWeighting *LumiWeights;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  // ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapToken_;
  // One example of full information about the cut flow
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleIdFullInfoMapToken_;
  edm::EDGetTokenT<edm::HepMCProduct> protonsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  HLTConfigProvider hltConfig_;
  HLTPrescaleProvider hltPrescaleProvider_;
  //edm::EDGetTokenT<edm::View<CTPPSLocalTrackLite> > spTracksToken_, recoTracksToken_;

  TH1F * h_trueNumInteractions;
  TH1F * h_mu_closest;
  TH1F * h_mu_closest_chi2_10;
  TH2F * h_mu_chi2_vs_closest;
  TH1F * h_e_closest;
  TH1F * h_e_closest_chi2_10;
  TH2F * h_e_chi2_vs_closest;
  TH1F * h_lepton_pt;

  std::vector<float> * hepmc_px_;
  std::vector<float> * hepmc_py_;
  std::vector<float> * hepmc_pz_;
  std::vector<float> * hepmc_energy_;
  std::vector<float> * hepmc_vx_;
  std::vector<float> * hepmc_vy_;
  std::vector<float> * hepmc_vz_;


  std::vector<float> * muon_pt_;
  std::vector<float> * muon_eta_;
  std::vector<float> * muon_px_;
  std::vector<float> * muon_py_;
  std::vector<float> * muon_pz_;
  std::vector<float> * muon_e_;
  std::vector<float> * muon_charge_;
  std::vector<float> * muon_iso_;

  std::vector<float> * electron_pt_;
  std::vector<float> * electron_eta_;
  std::vector<float> * electron_px_;
  std::vector<float> * electron_py_;
  std::vector<float> * electron_pz_;
  std::vector<float> * electron_e_;
  std::vector<float> * electron_charge_;
  std::vector<bool> *electron_passip_;

  std::vector<float> * allvertices_z_;

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
  std::vector<float> * muon_tkdist_;
  std::vector<float> * muon_tkpt_;
  std::vector<float> * electron_tkdist_;
  std::vector<float> * electron_tkpt_;

  //std::vector<float> * rp_tracks_xraw_;
  std::vector<float> * rp_tracks_y_;
  std::vector<float> * rp_tracks_x_;
  std::vector<float> * rp_tracks_tx_;
  std::vector<float> * rp_tracks_ty_;
  std::vector<float> * rp_tracks_xi_;
  std::vector<float> * rp_tracks_xi_unc_;
  std::vector<float> * rp_tracks_detId_;
  //float * mumu_mass_;
  //float * mumu_rapidity_;

  uint * run_;
  uint * ev_;
  uint * lumiblock_;
  bool * ispps_;

  //float * Tnpv_;
  float * pileupWeight_;

};
