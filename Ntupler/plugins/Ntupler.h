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
#include "TRandom3.h"
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
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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
  bool isJetLepton(double,double);
  bool isJetLeptonAK8(double,double);
  double calculateJERFactor(double,double,double,double,double,bool);

  // ----------member data ---------------------------
  boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  boost::shared_ptr<FactorizedJetCorrector> jecAK8_withL1_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jet_token_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jetAK8_token_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muon_token_;
  edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite> > pps_token_;
  edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsSingleRPToken_;
  edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsMultiRPToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertex_token_;
  edm::EDGetTokenT<double> rho_token_;
  edm::EDGetTokenT<edm::TriggerResults> hlt_token_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > > pu_token_;
  edm::EDGetTokenT<reco::GenParticleCollection> gen_part_token_;
  edm::EDGetTokenT<reco::GenJetCollection> gen_jet_token_;
  edm::EDGetTokenT<pat::METCollection> met_token_;
  edm::EDGetTokenT< bool >ecalBadCalibFilterUpdate_token;
  edm::EDGetTokenT<edm::View<pat::Electron>> electron_token_;
  //edm::EDGetTokenT<edm::View<reco::GsfElectron>> electron_token_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfcand_token_;
  edm::EDGetTokenT<GenEventInfoProduct> mcweight_token_;

  //edm::EDGetToken electronsMiniAODToken_;
  // ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapToken_veto_;

  std::string jerAK8chsName_res_ ;
  std::string jerAK8chsName_sf_ ;
  std::string jerAK4chsName_res_ ;
  std::string jerAK4chsName_sf_ ;


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

  std::vector<float> * electron_pt_;
  std::vector<float> * electron_eta_;
  std::vector<float> * electron_phi_;
  std::vector<float> * electron_dxy_;
  std::vector<float> * electron_dz_;
  std::vector<float> * electron_px_;
  std::vector<float> * electron_py_;
  std::vector<float> * electron_pz_;
  std::vector<float> * electron_e_;
  std::vector<float> * electron_charge_;

  float * met_;
  float * met_x_;
  float * met_y_;
  float * met_phi_;

  int * num_bjets_ak8_;
  int * num_bjets_ak4_;
  int * num_jets_ak4_;

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
  std::vector<float> * jet_jer_res_;
  std::vector<float> * jet_jer_sf_;
  std::vector<float> * jet_jer_sfup_;
  std::vector<float> * jet_jer_sfdown_;

  std::vector<float> * pps_track_x_;
  std::vector<float> * pps_track_y_;
  std::vector<int> * pps_track_rpid_;

  std::vector<float> * gen_W_pt_;
  std::vector<float> * gen_W_charge_;

  std::vector<float> * gen_proton_px_;
  std::vector<float> * gen_proton_py_;
  std::vector<float> * gen_proton_pz_;
  std::vector<float> * gen_proton_energy_;
  std::vector<float> * gen_proton_xi_;
  std::vector<float> * gen_proton_t_;

  std::vector<float> * proton_xi_;
  std::vector<float> * proton_thy_;
  std::vector<float> * proton_thx_;
  std::vector<float> * proton_t_;
  std::vector<int> * proton_ismultirp_;
  std::vector<int> * proton_rpid_;
  std::vector<int> * proton_arm_;

  std::vector<float> * proton_time_;
  std::vector<float> * proton_trackx1_;
  std::vector<float> * proton_tracky1_;
  std::vector<float> * proton_trackx2_;
  std::vector<float> * proton_tracky2_;
  std::vector<int> * proton_trackpixshift1_;
  std::vector<int> * proton_trackpixshift2_;
  std::vector<int> * proton_rpid1_;
  std::vector<int> * proton_rpid2_;

  std::vector<float> * gen_jet_pt_;
  std::vector<float> * gen_jet_eta_;
  std::vector<float> * gen_jet_phi_;
  std::vector<float> * gen_jet_energy_;

  std::vector<string> * hlt_;

  int * run_;
  long int * ev_;
  int * lumiblock_;
  int * nVertices_;
  float * pileupWeight_;
  float * mc_pu_trueinteractions_;
  float * mcWeight_;
  int * pfcand_nextracks_;
  int * pfcand_nextracks_noDRl_;

  float * recoMWhad_;
  float * recoMWlep_;
  float * dphiWW_;
  float * recoMWW_;
  float * recoRapidityWW_;
  float * WLeptonicPt_;
  float * WLeptonicPhi_;

  //bool * ecalBadCalFilter_;

  HLTConfigProvider hltConfig_;
  HLTPrescaleProvider hltPrescaleProvider_;
  edm::LumiReWeighting *LumiWeights;

  bool isMC;
  bool isSignalMC;
  bool isInteractive;
  int year;
  std::string era;
  std::string mcName;

  

};
