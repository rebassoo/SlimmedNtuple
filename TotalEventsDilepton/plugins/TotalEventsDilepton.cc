// -*- C++ -*-
//
// Package:    SlimmedNtuple/TotalEventsDilepton
// Class:      TotalEventsDilepton
// 
/**\class TotalEventsDilepton TotalEventsDilepton.cc SlimmedNtuple/TotalEventsDilepton/plugins/TotalEventsDilepton.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Fri, 21 Sep 2018 20:37:34 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include "TH2F.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TotalEventsDilepton : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TotalEventsDilepton(const edm::ParameterSet&);
      ~TotalEventsDilepton();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  

      // ----------member data ---------------------------
      TH1F * h_total_events;
  
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
TotalEventsDilepton::TotalEventsDilepton(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


TotalEventsDilepton::~TotalEventsDilepton()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TotalEventsDilepton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   h_total_events->Fill(1.);
   //cout<<"Get to this event"<<endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
TotalEventsDilepton::beginJob()
{
  edm::Service<TFileService> fs;
  h_total_events = fs->make<TH1F>("h_total_events","",2 ,-0.5 ,1.5);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TotalEventsDilepton::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TotalEventsDilepton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TotalEventsDilepton);
