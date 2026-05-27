// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      ObjectTrackCountProducer
//
/**\class ObjectTrackCountProducer ObjectTrackCountProducer.cc PhysicsTools/NanoAOD/plugins/ObjectTrackCountProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Pitt
//         Created:  Mon, 30 Mar 2026 11:11:03 GMT
//
//

// system include files
#include <memory>
#include <vector>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class ObjectTrackCountProducer : public edm::stream::EDProducer<> {
public:
  explicit ObjectTrackCountProducer(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  const edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  const edm::EDGetTokenT<pat::JetCollection> jetToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<pat::ElectronCollection> elecToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfcToken_;

  const double dzOutCut_;
  const double dzCut_;
  const double etaCut_;
  const double drJet_;
  const double drEl_;  
  const double drMu_;  

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
ObjectTrackCountProducer::ObjectTrackCountProducer(const edm::ParameterSet& iConfig) 
    : vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxs"))),
      jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
      muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      elecToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
      pfcToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcands"))),
      dzOutCut_(iConfig.getParameter<double>("dzOutCut")),
      dzCut_(iConfig.getParameter<double>("dzCut")),
      etaCut_(iConfig.getParameter<double>("etaCut")),
      drJet_(iConfig.getParameter<double>("drJet")),
      drEl_(iConfig.getParameter<double>("drEl")),
      drMu_(iConfig.getParameter<double>("drMu")) {

  produces<nanoaod::FlatTable>("pvTrk");
  produces<nanoaod::FlatTable>("otherPVTrk");
  produces<nanoaod::FlatTable>("jetTrk");
  produces<nanoaod::FlatTable>("muonTrk");
  produces<nanoaod::FlatTable>("elecTrk");
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void ObjectTrackCountProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  
  auto const& vtxs = iEvent.get(vtxToken_);
  auto const& jets = iEvent.get(jetToken_);
  auto const& muons = iEvent.get(muonToken_);
  auto const& elecs = iEvent.get(elecToken_);
  auto const& pfcands = iEvent.get(pfcToken_);
  
  if (vtxs.empty()) {
      iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "PV", true), "pvTrk");
      iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "OtherPV", true), "otherPVTrk");
      iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "Jet", true), "jetTrk");
      iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "Muon", true), "muonTrk");
      iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "Electron", true), "elecTrk");
      return;
  }
  const auto& pv0 = vtxs.front();

  // Count tracks for vertices:
  std::vector<int> v05, v09;
  std::vector<int> vOut05, vOut09;
  size_t nVtxToProcess = std::min<size_t>(vtxs.size(), 4);
  
  for (size_t i = 0; i < nVtxToProcess; ++i) {
	  const auto& v = vtxs[i];
	  int n05 = 0, n09 = 0;
	  int nOut05 = 0, nOut09 = 0;
	  for (const auto& pfc : pfcands) {
        if (pfc.charge() == 0 || pfc.pt() < 0.5 || std::abs(pfc.eta()) > etaCut_) continue;
		
		// Count tracks that are strictly OUTSIDE the main interaction window
		if (std::abs(pfc.dz(v.position())) > dzOutCut_) {
			nOut05++;
			if (pfc.pt() > 0.9) nOut09++;
		}
		
		// track - vertex association (NoPV = 0, PVLoose = 1, PVTight = 2, PVUsedInFit = 3):
		if (pfc.fromPV(i)<pat::PackedCandidate::PVTight) continue;
      
        // Master dZ cut relative to THIS specific vertex (v)
        if (std::abs(pfc.dz(v.position())) < dzCut_) {
          n05++;
          if (pfc.pt() > 0.9) n09++;
        }
      }
	  v05.push_back(n05);
	  v09.push_back(n09);
	  vOut05.push_back(nOut05);
	  vOut09.push_back(nOut09);
  }
  
  
  // Helper lambda for footprint counting
  auto countFootprint = [&](const auto& obj, double dRMax) {
    int n05 = 0, n09 = 0;
    for (const auto& pfc : pfcands) {
      if (pfc.charge() == 0 || pfc.pt() < 0.5) continue;
	  
      // Require association to the primary vertex
	  if (pfc.fromPV(0)<pat::PackedCandidate::PVTight) continue;
	  
	  // The Master dz Cut
	  if (std::abs(pfc.dz(pv0.position())) >= dzCut_) continue;
	  
	  // Reject tracks outside a dR cone from the object
      if (reco::deltaR(obj, pfc) >= dRMax) continue;
	  
	  // Counter
	  n05++;
	  if (pfc.pt() > 0.9) n09++;
    }
    return std::make_pair(n05, n09);
  };

  // Process Collections
  std::vector<int> j05, j09, m05, m09, e05, e09;
  for (const auto& j : jets) { auto res = countFootprint(j, drJet_); j05.push_back(res.first); j09.push_back(res.second); }
  for (const auto& m : muons) { auto res = countFootprint(m, drMu_); m05.push_back(res.first); m09.push_back(res.second); }
  for (const auto& e : elecs) { auto res = countFootprint(e, drEl_); e05.push_back(res.first); e09.push_back(res.second); }
  
  // Fill Tables
  auto pvTab = std::make_unique<nanoaod::FlatTable>(1, "PV", true, true);
  pvTab->addColumnValue<int>("ntrk0p5", v05.at(0), "Total PV tracks pt>0.5");
  pvTab->addColumnValue<int>("ntrk0p9", v09.at(0), "Total PV tracks pt>0.9");
  pvTab->addColumnValue<int>("nPUtrk0p5", vOut05.at(0), "Total PV tracks pt>0.5 outside |dz| > 2cm from PV0");
  pvTab->addColumnValue<int>("nPUtrk0p9", vOut09.at(0), "Total PV tracks pt>0.9 outside |dz| > 2cm from PV0");
  
  size_t nOther = (v05.size() > 1) ? v05.size() - 1 : 0;
  auto otherTab = std::make_unique<nanoaod::FlatTable>(nOther, "OtherPV", false, true);
  std::vector<int> ov05(v05.begin() + 1, v05.end());
  std::vector<int> ov09(v09.begin() + 1, v09.end());
  std::vector<int> ovOut05(vOut05.begin() + 1, vOut05.end());
  std::vector<int> ovOut09(vOut09.begin() + 1, vOut09.end());
  if (nOther > 0) {
	  ov05.assign(v05.begin() + 1, v05.end());
	  ov09.assign(v09.begin() + 1, v09.end());
	  ovOut05.assign(vOut05.begin() + 1, vOut05.end());
	  ovOut09.assign(vOut09.begin() + 1, vOut09.end());
  }
  otherTab->addColumn<int>("ntrk0p5", ov05, "Other PV tracks pt>0.5");
  otherTab->addColumn<int>("ntrk0p9", ov09, "Other PV tracks pt>0.9");
  otherTab->addColumn<int>("nPUtrk0p5", ovOut05, "Other PV tracks pt>0.5 outside |dz| > 2cm from PV");
  otherTab->addColumn<int>("nPUtrk0p9", ovOut09, "Other PV tracks pt>0.9 outside |dz| > 2cm from PV");

  auto jetTab = std::make_unique<nanoaod::FlatTable>(jets.size(), "Jet", false, true);
  jetTab->addColumn<int>("ntrk0p5", j05, "Jet track footprint pt>0.5");
  jetTab->addColumn<int>("ntrk0p9", j09, "Jet track footprint pt>0.9");

  auto muonTab = std::make_unique<nanoaod::FlatTable>(muons.size(), "Muon", false, true);
  muonTab->addColumn<int>("ntrk0p5", m05, "Muon track footprint pt>0.5");
  muonTab->addColumn<int>("ntrk0p9", m09, "Muon track footprint pt>0.9");

  auto elecTab = std::make_unique<nanoaod::FlatTable>(elecs.size(), "Electron", false, true);
  elecTab->addColumn<int>("ntrk0p5", e05, "Elec track footprint pt>0.5");
  elecTab->addColumn<int>("ntrk0p9", e09, "Elec track footprint pt>0.9");

  iEvent.put(std::move(pvTab),    "pvTrk");
  iEvent.put(std::move(otherTab), "otherPVTrk");
  iEvent.put(std::move(jetTab),   "jetTrk");
  iEvent.put(std::move(muonTab),  "muonTrk");
  iEvent.put(std::move(elecTab),  "elecTrk");
}


// ------------ method called when starting to processes a run  ------------
/*
void
ObjectTrackCountProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
ObjectTrackCountProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ObjectTrackCountProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ObjectTrackCountProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ObjectTrackCountProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("vtxs", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("pfcands", edm::InputTag("packedPFCandidates"));
  
  // Define default values here
  desc.add<double>("dzOutCut", 0.2);
  desc.add<double>("dzCut", 0.1);
  desc.add<double>("etaCut", 2.1);
  desc.add<double>("drJet", 0.4);
  desc.add<double>("drEl", 0.2);
  desc.add<double>("drMu", 0.1);
  
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ObjectTrackCountProducer);
