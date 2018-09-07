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

#include "TLorentzVector.h"

#include <iostream>       // std::cout
#include <string>         // std::string


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "CalculatePzNu.hpp"

//
// class declaration
//

class Ntupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  //class Ntupler : public edm::EDAnalyzer {
   public:
      explicit Ntupler(const edm::ParameterSet&);
      ~Ntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  //  virtual void endRun(edm::Run const &, edm::EventSetup const&) override;
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void GetProtons(const edm::Event&);
  //void GetMC(const edm::Event&);
  //void GetMuons(const edm::Event&,reco::VertexRef,edm::ESHandle<TransientTrackBuilder>&,std::vector<reco::TransientTrack>&,std::vector<uint>&,std::vector<reco::TransientTrack>,int&);
  //void GetElectrons(const edm::Event&,reco::VertexRef,edm::ESHandle<TransientTrackBuilder>&,std::vector<reco::TransientTrack>&,std::vector<reco::TransientTrack>&,std::vector<uint>&,std::vector<reco::TransientTrack>,int&);
  //void GetJets(const edm::Event&);
  bool isJetLepton(double,double);

  // ----------member data ---------------------------
  boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jet_token_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jetAK8_token_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muon_token_;
  edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite> > pps_token_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertex_token_;
  edm::EDGetTokenT<double> rho_token_;
  edm::EDGetTokenT<edm::TriggerResults> hlt_token_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > > pu_token_;
  edm::EDGetTokenT<reco::GenParticleCollection> gen_part_token_;
  edm::EDGetTokenT<reco::GenJetCollection> gen_jet_token_;
  edm::EDGetTokenT<edm::View<pat::MET>> met_token_;
  edm::EDGetTokenT<edm::View<pat::Electron>> electron_token_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfcand_token_;

  TTree * tree_;

  std::vector<float> * muon_pt_;
  std::vector<float> * muon_eta_;
  std::vector<float> * muon_phi_;
  std::vector<float> * muon_px_;
  std::vector<float> * muon_py_;
  std::vector<float> * muon_pz_;
  std::vector<float> * muon_e_;
  std::vector<float> * muon_charge_;
  std::vector<float> * muon_iso_;
  std::vector<float> * muon_dxy_;
  std::vector<float> * muon_dz_;

  float * met_;
  float * met_x_;
  float * met_y_;

  std::vector<float> * jet_pt_;
  std::vector<float> * jet_energy_;
  std::vector<float> * jet_phi_;
  std::vector<float> * jet_eta_;
  std::vector<float> * jet_px_;
  std::vector<float> * jet_py_;
  std::vector<float> * jet_pz_;
  std::vector<float> * jet_mass_;
  std::vector<float> * jet_tau1_;
  std::vector<float> * jet_tau2_;
  std::vector<float> * jet_corrmass_;
  std::vector<float> * jet_vertexz_;

  std::vector<float> * dijet_mass_;
  std::vector<float> * dijet_pt_;
  std::vector<float> * dijet_y_;
  std::vector<float> * dijet_phi_;
  std::vector<float> * dijet_dphi_;

  std::vector<float> * pps_track_x_;
  std::vector<float> * pps_track_y_;
  std::vector<int> * pps_track_rpid_;

  std::vector<float> * gen_proton_px_;
  std::vector<float> * gen_proton_py_;
  std::vector<float> * gen_proton_pz_;
  std::vector<float> * gen_proton_xi_;
  std::vector<float> * gen_proton_t_;

  std::vector<string> * hlt_;

  int * run_;
  long int * ev_;
  int * lumiblock_;
  int * nVertices_;
  float * pileupWeight_;
  int * pfcand_nextracks_;

  float * recoMWhad_;
  float * recoMWlep_;
  float * dphiWW_;
  float * recoMWW_;
  float * WLeptonicPt_;

  HLTConfigProvider hltConfig_;
  HLTPrescaleProvider hltPrescaleProvider_;
  edm::LumiReWeighting *LumiWeights;

  bool isMC;
  int year;
  std::string era;

  

};
