// -*- C++ -*-
//
// Package:    L1Trigger/L1TCommon
// Class:      MyNTuplizer
// 
/**\class MyNTuplizer MyNTuplizer.cc L1Trigger/L1TCommon/plugins/MyNTuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Pitt
//         Created:  Thu, 20 Feb 2020 11:30:38 GMT
//
//

#include <iostream>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/L1TCorrelator/interface/TkElectron.h"
#include "DataFormats/L1TCorrelator/interface/TkElectronFwd.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"

// ROOT includes
#include <TTree.h>

#define ETCUT 0 // 1.47 
#define MAXETA 3.5  
#define ZERO_PT 0.05

using namespace std;
using namespace edm;
using namespace l1t;
using namespace reco;

namespace {
	enum ParticleType{ Electron = 11, Photon = 22 };
	enum HWQUALITY{ ENDCUPNOBREM = 4, ENDCUPBREM = 5 };  
}

class MyNTuplizer : public EDAnalyzer {
public:
  explicit MyNTuplizer(const ParameterSet&);
  ~MyNTuplizer() override;

  static void fillDescriptions(ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(Event const&, EventSetup const&) override;
  void endJob() override;

  void beginRun(Run const&, EventSetup const&) override;
  void endRun(Run const&, EventSetup const&) override;
  void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
  void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;

  // EDM tokens:
  EDGetTokenT<std::vector<GenParticle>> genparticles_;
  EDGetTokenT<EGammaBxCollection> egToken_;
  EDGetTokenT<CandidateView> epfToken_;
  EDGetTokenT<CandidateView> eegpfToken_;
  EDGetTokenT<TkElectronCollection> tkElToken_;
  
  TTree *tree_; 
  uint32_t run_, lumi_; uint64_t event_;
  
  Service<TFileService> fs;
  
	  struct McVars {
         int n = 0;
         std::vector<float> pt;
		 std::vector<float> e;
		 std::vector<float> caloeta;
		 std::vector<float> eta;
		 std::vector<float> phi;
		 std::vector<float> calophi;
		 std::vector<int> id;
		 std::vector<int> q;
         void makeBranches(TTree *tree) {
		    tree->Branch("mc_n",     &n, "mc_n/I");
		    tree->Branch("mc_pt",     &pt);
		    tree->Branch("mc_eta",    &eta);
		    tree->Branch("mc_caloeta",    &caloeta);
		    tree->Branch("mc_phi",    &phi);
		    tree->Branch("mc_calophi",    &calophi);
		    tree->Branch("mc_q", &q);
		    tree->Branch("mc_id", &id);
         }
		 bool isBarrel(float eta){ return (fabs(eta)<ETCUT);}
		 void clear(){n=0; pt.clear(); caloeta.clear(); eta.clear(); phi.clear(); calophi.clear(); id.clear(); q.clear();}		 
         void fill(const reco::GenParticleCollection &genparticles) {
			 clear();
			 for (const reco::GenParticle &gen : genparticles) {
				 if(gen.status()!=1) continue; // keep only stable particles
				 float _eta = gen.eta();
				 if(isBarrel(_eta)) continue;        // only EndCup objects
				 if(fabs(_eta)>MAXETA) continue;        // remove forward particles
				 float _pt = gen.pt();
				 if(_pt<ZERO_PT) continue;				 // skip 
				 n++;
                 pt.emplace_back(_pt);
			     eta.emplace_back(_eta);
			     phi.emplace_back(gen.phi());
			     id.emplace_back(gen.pdgId());
				 int _q = gen.charge();
				 float _caloeta = _eta;
				 float _calophi = gen.phi();
			     q.emplace_back(_q);
				 if(_q){
					 math::XYZTLorentzVector vertex(gen.vx(),gen.vy(),gen.vz(),0.);
					 auto caloetaphi = l1tpf::propagateToCalo(gen.p4(),vertex, _q, 3.8112);
					 _caloeta = caloetaphi.first; _calophi = caloetaphi.second;
				 }
			     caloeta.emplace_back(_caloeta);
			     calophi.emplace_back(_calophi);				 				 
			 }
         }
      } mc_;	

	  struct PFlowVars {
		 int n = 0;
         std::vector<float> pt;
		 std::vector<float> e;
		 std::vector<float> eta;
		 std::vector<float> phi;
		 std::vector<int> q;
		 std::vector<int> id;
         void makeBranches(TTree *tree) {
		    tree->Branch("PFEgamma_n",     &n, "PFEgamma_n/I");
		    tree->Branch("PFEgamma_pt",     &pt);
		    tree->Branch("PFEgamma_energy", &e);
		    tree->Branch("PFEgamma_eta",    &eta);
		    tree->Branch("PFEgamma_phi",    &phi);
		    tree->Branch("PFEgamma_q", &q);
		    tree->Branch("PFEgamma_type", &id);
         }
		 bool isBarrel(float eta){ return (fabs(eta)<ETCUT);}
		 void clear(){n=0; pt.clear(); e.clear(); eta.clear(); phi.clear(); q.clear(); id.clear();}
         void fill(const reco::CandidateView &ecals) {
			 clear();
			 for (const reco::Candidate & c : ecals) {
				 //if(!c.isElectron()) continue; // not working!! :/
				 float _eta = c.eta();
				 if(fabs(c.pdgId())!=11 && fabs(c.pdgId())!=22) continue; // keep only Egamma
				 //if(isBarrel(_eta)) continue;        // only EndCup objects
				 float _pt = c.pt();
				 if(_pt<ZERO_PT) continue;				 // skip zero energy clusters
				 n++;
                 pt.emplace_back(_pt);
			     e.emplace_back(c.energy());
			     eta.emplace_back(_eta);
			     phi.emplace_back(c.phi());
			     q.emplace_back(c.charge());
			     id.emplace_back(fabs(c.pdgId()));
			 }
         }
      } pf_;	
	  
      struct EGPFlowVars {
		 int n = 0;
         std::vector<float> pt;
		 std::vector<float> e;
		 std::vector<float> eta;
		 std::vector<float> phi;
		 std::vector<int> q;
		 std::vector<int> id;
         void makeBranches(TTree *tree) {
		    tree->Branch("EGPFEgamma_n",     &n, "EGPFEgamma_n/I");
		    tree->Branch("EGPFEgamma_pt",     &pt);
		    tree->Branch("EGPFEgamma_energy", &e);
		    tree->Branch("EGPFEgamma_eta",    &eta);
		    tree->Branch("EGPFEgamma_phi",    &phi);
		    tree->Branch("EGPFEgamma_q", &q);
		    tree->Branch("EGPFEgamma_type", &id);
         }
		 bool isBarrel(float eta){ return (fabs(eta)<ETCUT);}
		 void clear(){n=0; pt.clear(); e.clear(); eta.clear(); phi.clear(); q.clear(); id.clear();}
         void fill(const reco::CandidateView &ecals) {
			 clear();
			 for (const reco::Candidate & c : ecals) {
				 float _eta = c.eta();
				 float _pt = c.pt();
				 if(_pt<ZERO_PT) continue;				 // skip zero energy clusters
				 n++;
                 pt.emplace_back(_pt);
			     e.emplace_back(c.energy());
			     eta.emplace_back(_eta);
			     phi.emplace_back(c.phi());
			     q.emplace_back(c.charge());
			     id.emplace_back(fabs(c.pdgId()));
			 }
         }
      } egpf_;		  
	  
	  struct EgammaEEVars {
		 int n = 0;
         vector<float> pt;
		 vector<float> e;
		 vector<float> eta;
		 vector<float> phi;
		 vector<float> hwQ;
         void makeBranches(TTree *tree) {
		    tree->Branch("egammaEE_n",     &n, "egammaEE_n/I");
		    tree->Branch("egammaEE_pt",     &pt);
		    tree->Branch("egammaEE_energy", &e);
		    tree->Branch("egammaEE_eta",    &eta);
		    tree->Branch("egammaEE_phi",    &phi);
		    tree->Branch("egammaEE_hwQual", &hwQ);
         }
		 void clear() {n=0; pt.clear(); e.clear(); eta.clear(); phi.clear();   hwQ.clear();}
         void fill(const EGammaBxCollection &egamma_ee_collection) {
			 clear();
			 for(auto egee_itr=egamma_ee_collection.begin(0); egee_itr!=egamma_ee_collection.end(0); egee_itr++) {
			    // for endcup using hwQual=5 or 4 (bremm recovery)
				if(egee_itr->hwQual()!=4 && egee_itr->hwQual()!=5) continue;
				float _pt = egee_itr->pt();
				if(_pt<ZERO_PT) continue;				 // skip zero energy clusters
				n++;
                pt.emplace_back(_pt);
			    e.emplace_back(egee_itr->energy());
			    eta.emplace_back(egee_itr->eta());
			    phi.emplace_back(egee_itr->phi());
			    hwQ.emplace_back(egee_itr->hwQual());
			 }
         }
      } egamma_;  

	  struct l1tkElectron {
		 int n = 0;
         std::vector<float> pt;
		 std::vector<float> eta;
		 std::vector<float> phi;
		 std::vector<int> 	q;
		 std::vector<float> hwQ;
         void makeBranches(TTree *tree) {
		    tree->Branch("l1tkElectron_n",     &n, "l1tkElectron_n/I");
		    tree->Branch("l1tkElectron_pt",     &pt);
		    tree->Branch("l1tkElectron_eta",    &eta);
		    tree->Branch("l1tkElectron_phi",    &phi);
		    tree->Branch("l1tkElectron_q	",    &q);
		    tree->Branch("l1tkElectron_hwQual	",    &hwQ);
         }
		 void clear() {n=0; pt.clear(); eta.clear(); phi.clear(); q.clear();  hwQ.clear();}
         void fill(const TkElectronCollection &egamma_collection) {
			 // for endcup using hwQual=5
			 clear();
			 for(auto egTrkIter=egamma_collection.begin(); egTrkIter!=egamma_collection.end(); egTrkIter++) {
				 //if(egTrkIter->getEGRef()->hwQual()!=hwQual) continue;
				 int hwQual = egTrkIter->EGRef()->hwQual();
				 if(5!=hwQual && 4!=hwQual) continue;
				 float _pt = egTrkIter->pt();
				 if(_pt<ZERO_PT) continue;				 // skip zero energy tracks
				 n++;
                 pt.emplace_back(_pt);
			     eta.emplace_back(egTrkIter->eta());
			     phi.emplace_back(egTrkIter->phi());
			     q.emplace_back(egTrkIter->charge());
			     hwQ.emplace_back(float(hwQual));
			 }
         }
      } l1tkElectron_;
};

MyNTuplizer::MyNTuplizer(const ParameterSet& iConfig) : 
genparticles_(consumes<vector<GenParticle>>(iConfig.getParameter<InputTag>("genParticles"))),
egToken_(consumes<EGammaBxCollection>(iConfig.getParameter<InputTag>("EgammaEEInputTag"))),
epfToken_(consumes<CandidateView>(iConfig.getParameter<InputTag>("PFEgammaInputTag"))),
eegpfToken_(consumes<CandidateView>(iConfig.getParameter<InputTag>("EGPFEgammaInputTag"))),
tkElToken_(consumes< TkElectronCollection > (iConfig.getParameter<InputTag>("L1TkElectronInputTag")))
{

  tree_ = fs->make<TTree>("tree","tree");
  tree_->Branch("run",  &run_,  "run/i");
  tree_->Branch("lumi", &lumi_, "lumi/i");
  tree_->Branch("event", &event_, "event/l");
}

MyNTuplizer::~MyNTuplizer() {}

void MyNTuplizer::analyze(Event const& iEvent, EventSetup const& iSetup) {

    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();

	Handle<EGammaBxCollection> egamma_ee_h;
	Handle<GenParticleCollection> genparticles;
	Handle<CandidateView> pf_egamma_;
	Handle<CandidateView> egpf_egamma_;
	Handle<TkElectronCollection> l1tkElectronHandle;

	iEvent.getByToken(egToken_, egamma_ee_h);
	iEvent.getByToken(epfToken_, pf_egamma_);
	iEvent.getByToken(eegpfToken_, egpf_egamma_);
	iEvent.getByToken(genparticles_, genparticles);
	iEvent.getByToken(tkElToken_, l1tkElectronHandle);
	
	mc_.fill(*genparticles);
	egamma_.fill(*egamma_ee_h);
	pf_.fill(*pf_egamma_);
	egpf_.fill(*egpf_egamma_);
	l1tkElectron_.fill(*l1tkElectronHandle.product());
	
    //if(mc_.n || egamma_.n || pf_.n) 
    tree_->Fill(); // so that we don't write events with zero leptons

}

void MyNTuplizer::beginJob() { 
	cout << "INFO:  MyNTuplizer module beginJob called.\n";
	egamma_.makeBranches(tree_);
	l1tkElectron_.makeBranches(tree_);
	mc_.makeBranches(tree_);
	pf_.makeBranches(tree_);
	egpf_.makeBranches(tree_);
}


void MyNTuplizer::endJob() {
}

void MyNTuplizer::beginRun(Run const& run, EventSetup const& iSetup) {
}

void MyNTuplizer::endRun(Run const&, EventSetup const&) {}

void MyNTuplizer::beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) {}

void MyNTuplizer::endLuminosityBlock(LuminosityBlock const&, EventSetup const&) {}

void MyNTuplizer::fillDescriptions(ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MyNTuplizer);
