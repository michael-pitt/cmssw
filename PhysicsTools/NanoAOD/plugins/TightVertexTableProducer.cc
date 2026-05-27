// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      TightVertexTableProducer
//
/**\class TightVertexTableProducer TightVertexTableProducer.cc PhysicsTools/NanoAOD/plugins/TightVertexTableProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Pitt
//         Created:  Sun, 24 May 2026 13:12:56 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
//
// class declaration
//

class TightVertexTableProducer : public edm::stream::EDProducer<> {
public:
  explicit TightVertexTableProducer(const edm::ParameterSet&);
  ~TightVertexTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  const StringCutObjectSelector<reco::Vertex> goodPvCut_;
};

//
// constructors and destructor
//
TightVertexTableProducer::TightVertexTableProducer(const edm::ParameterSet& iConfig):
    // Initialize the tokens to read from the python configuration
    vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    goodPvCut_(iConfig.getParameter<std::string>("goodPvCut"), true)
  {
  //register your products
  produces<nanoaod::FlatTable>("tightPV");

  //now do what ever other initialization is needed
}

TightVertexTableProducer::~TightVertexTableProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void TightVertexTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Get the primary vertices from the event
  Handle<std::vector<reco::Vertex>> vtxs;
  iEvent.getByToken(vtxToken_, vtxs);

  int tight = 0;

  // Loop over vertices and count
  if (vtxs.isValid()) {
    for (const auto& v : *vtxs) {
      // Check if the vertex passes the standard CMS "Good" criteria first
      if (goodPvCut_(v)) {
        // Apply diffractive-specific tracking cuts
        if (v.ndof() > 10.0) tight++;
      }
    }
  }

  // Fill Tables
  auto pvTableExtension = std::make_unique<nanoaod::FlatTable>(1, "PV", true, true);
  pvTableExtension->addColumnValue<uint8_t>("npvsTight", tight, "number of good PVs with ndof > 10");

  // Put it in the event
  iEvent.put(std::move(pvTableExtension), "tightPV");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TightVertexTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlineSlimmedPrimaryVertices"))->setComment("primary vertex input collection");
  desc.add<std::string>("goodPvCut", "!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2")->setComment("selection on the primary vertex");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TightVertexTableProducer);
