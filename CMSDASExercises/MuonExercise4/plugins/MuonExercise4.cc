// -*- C++ -*-
//
// Package:    CMSDASExercises/MuonExercise4
// Class:      MuonExercise2015
// 
/**\class MuonExercise4 MuonExercise4.cc CMSDASExercises/MuonExercise4/plugins/MuonExercise4.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author: 
//         Created:  Thu, 10 Dec 2015 18:16:22 GMT
//
//


// system include files
#include <memory>

#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;

//
// class declaration
//

class MuonExercise4 : public edm::EDAnalyzer {

  public:

    explicit MuonExercise4(const edm::ParameterSet&);
    ~MuonExercise4();
  
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  private:

    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
  
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_; 
  // function declaration
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  bool MatchedToTriggerObject(const pat::Muon& mu, const edm::Event& ev, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<edm::TriggerResults> triggerBits, int & index );
  bool MatchedToGenObject(const pat::Muon& mu, Handle<edm::View<pat::PackedGenParticle> > packed,int & genindex );
  bool SoftMuon(const pat::Muon&  , const reco::Vertex& ) ;
  const reco::GenParticle getMother(const reco::GenParticle& gp);
  // histogram declaration
  edm::Service<TFileService> fs;
  const reco::GenParticle* mother;
  
  TH1F* hFillProbePt;
  TH1F* hFillProbeEta;
  TH1F* hFillProbePhi;
  
  TH1F* hFillProbePtLoose;
  TH1F* hFillProbeEtaLoose;
  TH1F* hFillProbePhiLoose;
  
  TH1F* hFillProbePtSoft;
  TH1F* hFillProbeEtaSoft;
  TH1F* hFillProbePhiSoft;

  TH1F* hFillProbePtTight;
  TH1F* hFillProbeEtaTight;
  TH1F* hFillProbePhiTight;
  
  TH1F* hFillProbePtGlob;
  TH1F* hFillProbeEtaGlob;
  TH1F* hFillProbePhiGlob;
  
  TH1F* hFillProbePtPF;
  TH1F* hFillProbeEtaPF;
  TH1F* hFillProbePhiPF;
  
  
  TH1F* hFillProbePtIso;
  TH1F* hFillProbeEtaIso;
  TH1F* hFillProbePhiIso;

 TH1F* hFillProbePtTrkIso;
 TH1F* hFillProbeEtaTrkIso;
 TH1F* hFillProbePhiTrkIso;
  
  TH1F* hFillProbePtTightIso;
  TH1F* hFillProbeEtaTightIso;
  TH1F* hFillProbePhiTightIso;

  TH1F* hFillProbePtTightTrkIso;
  TH1F* hFillProbeEtaTightTrkIso;
  TH1F* hFillProbePhiTightTrkIso;
 
  TH2F *hFillIsoVsMuPt;
  TH2F *hFillIsoVsMuPtTag;
  TH1F *hFillZMass;
  TH1F *hFillZMassTight;

  TH1F *hFillNVtx;
  TH1F *hFillNVtxPF;
  TH1F *hFillNVtxGlob;
  TH1F *hFillNVtxLoose;
  TH1F *hFillNVtxSoft;
  TH1F *hFillNVtxTight;
  TH1F *hFillNVtxTrkIso;
  TH1F *hFillNVtxIso;
  TH1F *hFillNVtxTightTrkIso;
  TH1F *hFillNVtxTightIso;

  double coriso= 99; 
  double coriso_p = 99; 
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

MuonExercise4::MuonExercise4(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
{
}


MuonExercise4::~MuonExercise4()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonExercise4::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  
  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);
  
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);
  
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  if (vertices->empty()) return;
  
  VertexCollection::const_iterator PV = vertices->end();
  
  for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
    if ( !(vtx->isFake())
	 && vtx->ndof()>=4. && vtx->position().Rho() < 2.0
	 && fabs(vtx->position().Z()) < 24.0) {
      PV = vtx;
      break;
    }
  }

  if ( PV==vertices->end() ) return; 

  // to calculate number of good reconstructed primary vertices.
  
  int ngood = 0;

  for (VertexCollection::const_iterator vtx = vertices->begin();vtx != vertices->end(); ++vtx) {
    if ( !(vtx->isFake())
         && vtx->ndof()>=4. && vtx->position().Rho()<=2.0
         && fabs(vtx->position().Z())<=24.0) {
      ngood++;
    }
  }
  
  
  // for muon id efficiency
  
  for (pat::MuonCollection::const_iterator muontag = muons->begin() ; muontag != muons->end() ; ++muontag) { 
    
    // for nutuples
    //
    // first lets find a tag ////////////////////
    
    if (muontag->pt() < 20.0) continue; 
    if (!(fabs(muontag->eta() < 2.4))) continue; 
    if (!(muontag->isTightMuon(*PV))) continue;
    
    // tag is isolated
    if (muontag->isIsolationValid())
      {
	reco::MuonPFIsolation pfR04 = muontag->pfIsolationR04();
	coriso = pfR04.sumChargedHadronPt + std::max(0., pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt);
	hFillIsoVsMuPtTag->Fill(muontag->pt(),coriso/muontag->pt());
      }
    
    if (!(coriso/muontag->pt() <  0.12 )) continue;
    
    // lets find all probes having Z mass window 70 < Mass < 110
    for (pat::MuonCollection::const_iterator muonprobe = muontag+1 ; muonprobe != muons->end() ; ++muonprobe) {

      
      if((muontag->charge() * muonprobe->charge()) != -1) continue; 
     
      double ZMass = ( muontag->p4() + muonprobe->p4()).M();
      
      if ( ZMass > 70 && ZMass < 110 && (muonprobe->innerTrack().isNonnull())) {

	// gen matching of probe

	if ((muonprobe->pt() > 20.) && (fabs(muonprobe->eta()) < 2.4)) {
	  
	  hFillProbePt->Fill(muonprobe->pt());
	  hFillProbeEta->Fill(muonprobe->eta());
	  hFillProbePhi->Fill(muonprobe->phi());
	  hFillZMass->Fill(ZMass);
	  hFillNVtx->Fill(ngood);
	  
	  if (muonprobe->isGlobalMuon()){  
	    hFillProbePtGlob->Fill(muonprobe->pt());
	    hFillProbeEtaGlob->Fill(muonprobe->eta());
	    hFillProbePhiGlob->Fill(muonprobe->phi());
	    hFillNVtxGlob->Fill(ngood); 
	  }// global muon
	  
	  if (muonprobe->isPFMuon()) {
	    hFillProbePtPF->Fill(muonprobe->pt());
	    hFillProbeEtaPF->Fill(muonprobe->eta());
	    hFillProbePhiPF->Fill(muonprobe->phi());
	    hFillNVtxPF->Fill(ngood);
	  } // pf muon
	  
	  if (muonprobe->isLooseMuon()) {
	     hFillProbePtLoose->Fill(muonprobe->pt());
	     hFillProbeEtaLoose->Fill(muonprobe->eta());
	     hFillProbePhiLoose->Fill(muonprobe->phi());  
	     hFillNVtxLoose->Fill(ngood);
	  } //loose
	   
	  if (SoftMuon(*muonprobe, *PV)){
	     hFillProbePtSoft->Fill(muonprobe->pt());
	     hFillProbeEtaSoft->Fill(muonprobe->eta());
	     hFillProbePhiSoft->Fill(muonprobe->phi());
	     hFillNVtxSoft->Fill(ngood);
	   } // soft
	   
	   if (muonprobe->isTightMuon(*PV)){
	     hFillProbePtTight->Fill(muonprobe->pt());
	     hFillProbeEtaTight->Fill(muonprobe->eta());
	     hFillProbePhiTight->Fill(muonprobe->phi());
	     hFillZMassTight->Fill(ZMass);
	     hFillNVtxTight->Fill(ngood);
	   }  
	   
	   // isolation efficiency is calculated for muons passing tight muon efficiency
	   
	   if (muonprobe->isIsolationValid()) {
	      
	       reco::MuonPFIsolation pfR04_p = muonprobe->pfIsolationR04();
	       coriso_p = pfR04_p.sumChargedHadronPt + std::max(0., pfR04_p.sumNeutralHadronEt+pfR04_p.sumPhotonEt-0.5*pfR04_p.sumPUPt);
	       hFillIsoVsMuPt->Fill(muonprobe->pt(),coriso_p/muonprobe->pt());
	       
	       if (coriso_p/muonprobe->pt() <  0.12 ) {
		 hFillProbePtIso->Fill(muonprobe->pt());
		 hFillProbeEtaIso->Fill(muonprobe->eta());
		 hFillProbePhiIso->Fill(muonprobe->phi());
		 hFillNVtxIso->Fill(ngood);
	       }
	       
	       if ( (muonprobe->isTightMuon(*PV)) && (coriso_p/muonprobe->pt() <  0.12) ) {
		 hFillProbePtTightIso->Fill(muonprobe->pt());
		 hFillProbeEtaTightIso->Fill(muonprobe->eta());
		 hFillProbePhiTightIso->Fill(muonprobe->phi());
		 hFillNVtxTightIso->Fill(ngood);
	       }
	       
	     } // closing isolation
	   
	   // track isolation
	   if (muonprobe->trackIso()/muonprobe->pt() <  0.05 ) {
	     hFillProbePtTrkIso->Fill(muonprobe->pt());
	     hFillProbeEtaTrkIso->Fill(muonprobe->eta());
	     hFillProbePhiTrkIso->Fill(muonprobe->phi());
	     hFillNVtxTrkIso->Fill(ngood);
	   }
	   
	   if ( (muonprobe->isTightMuon(*PV)) &&  (muonprobe->trackIso()/muonprobe->pt() <  0.05 )) {
	     hFillProbePtTightTrkIso->Fill(muonprobe->pt());
	     hFillProbeEtaTightTrkIso->Fill(muonprobe->eta());
	     hFillProbePhiTightTrkIso->Fill(muonprobe->phi());
	     hFillNVtxTightTrkIso->Fill(ngood);
	   }
	   
	} // probe pt and eta 
	
      }
    } // muon probe
  }// muon tag loop
  
}



// ------------ method called once each job just before starting event loop  ------------
void 
MuonExercise4::beginJob()
{

  // histograms here
  //
  TH1::SetDefaultSumw2();

  hFillIsoVsMuPtTag = fs->make<TH2F>("hFillIsoVsMuPtTag","hFillIsoVsMuPtTag",200,0,200,1000,0,2); 
  hFillIsoVsMuPt = fs->make<TH2F>("hFillIsoVsMuPt","hFillIsoVsMuPt",200,0,200,1000,0,2);

  hFillZMass = fs->make<TH1F>("hFillZMass","hFillZMass",500,60,120); 
  hFillZMassTight = fs->make<TH1F>("hFillZMassTight","hFillZMassTight",500,60,120);

  hFillProbePt  = fs->make<TH1F>("hFillProbePt","hFillProbePt",200,0,200); 
  hFillProbeEta = fs->make<TH1F>("hFillProbeEta","hFillProbeEta",200,-3,3);
  hFillProbePhi = fs->make<TH1F>("hFillProbePhi","hFillProbePhi",200,-3,3);
  
  hFillProbePtPF  = fs->make<TH1F>("hFillProbePtPF","hFillProbePtPF",200,0,200);
  hFillProbeEtaPF = fs->make<TH1F>("hFillProbeEtaPF","hFillProbeEtaPF",200,-3,3);
  hFillProbePhiPF = fs->make<TH1F>("hFillProbePhiPF","hFillProbePhiPF",200,-3,3);
  
  hFillProbePtGlob  = fs->make<TH1F>("hFillProbePtGlob","hFillProbePtGlob",200,0,200);
  hFillProbeEtaGlob = fs->make<TH1F>("hFillProbeEtaGlob","hFillProbeEtaGlob",200,-3,3);
  hFillProbePhiGlob = fs->make<TH1F>("hFillProbePhiGlob","hFillProbePhiGlob",200,-3,3);
  
  hFillProbePtLoose  = fs->make<TH1F>("hFillProbePtLoose","hFillProbePtLoose",200,0,200);
  hFillProbeEtaLoose = fs->make<TH1F>("hFillProbeEtaLoose","hFillProbeEtaLoose",200,-3,3);
  hFillProbePhiLoose = fs->make<TH1F>("hFillProbePhiLoose","hFillProbePhiLoose",200,-3,3);

  hFillProbePtSoft  = fs->make<TH1F>("hFillProbePtSoft","hFillProbePtSoft",200,0,200);
  hFillProbeEtaSoft = fs->make<TH1F>("hFillProbeEtaSoft","hFillProbeEtaSoft",200,-3,3);
  hFillProbePhiSoft = fs->make<TH1F>("hFillProbePhiSoft","hFillProbePhiSoft",200,-3,3);
  
  hFillProbePtTight  = fs->make<TH1F>("hFillProbePtTight","hFillProbePtTight",200,0,200);
  hFillProbeEtaTight = fs->make<TH1F>("hFillProbeEtaTight","hFillProbeEtaTight",200,-3,3);
  hFillProbePhiTight = fs->make<TH1F>("hFillProbePhiTight","hFillProbePhiTight",200,-3,3);
  
  
  hFillProbePtIso  = fs->make<TH1F>("hFillProbePtIso","hFillProbePtIso",200,0,200);
  hFillProbeEtaIso = fs->make<TH1F>("hFillProbeEtaIso","hFillProbeEtaIso",200,-3,3);
  hFillProbePhiIso = fs->make<TH1F>("hFillProbePhiIso","hFillProbePhiIso",200,-3,3);
  
  hFillProbePtTightIso  = fs->make<TH1F>("hFillProbePtTightIso","hFillProbePtTightIso",200,0,200);
  hFillProbeEtaTightIso = fs->make<TH1F>("hFillProbeEtaTightIso","hFillProbeEtaTightIso",200,-3,3);
  hFillProbePhiTightIso = fs->make<TH1F>("hFillProbePhiTightIso","hFillProbePhiTightIso",200,-3,3);
  
  hFillProbePtTrkIso  = fs->make<TH1F>("hFillProbePtTrkIso","hFillProbePtTrkIso",200,0,200);
  hFillProbeEtaTrkIso = fs->make<TH1F>("hFillProbeEtaTrkIso","hFillProbeEtaTrkIso",200,-3,3);
  hFillProbePhiTrkIso = fs->make<TH1F>("hFillProbePhiTrkIso","hFillProbePhiTrkIso",200,-3,3);
  
  hFillProbePtTightTrkIso  = fs->make<TH1F>("hFillProbePtTightTrkIso","hFillProbePtTightTrkIso",200,0,200);
  hFillProbeEtaTightTrkIso = fs->make<TH1F>("hFillProbeEtaTightTrkIso","hFillProbeEtaTightTrkIso",200,-3,3);
  hFillProbePhiTightTrkIso = fs->make<TH1F>("hFillProbePhiTightTrkIso","hFillProbePhiTightTrkIso",200,-3,3);

  hFillNVtx = fs->make<TH1F>("hFillNVtx","hFillNVtx",100,0,100);
  hFillNVtxPF = fs->make<TH1F>("hFillNVtxPF","hFillNVtxPF",100,0,100);
  hFillNVtxGlob = fs->make<TH1F>("hFillNVtxGlob","hFillNVtxGlob",100,0,100);
  hFillNVtxLoose = fs->make<TH1F>("hFillNVtxLoose","hFillNVtxLoose",100,0,100);
  hFillNVtxSoft = fs->make<TH1F>("hFillNVtxSoft","hFillNVtxSoft",100,0,100);
  hFillNVtxTight = fs->make<TH1F>("hFillNVtxTight","hFillNVtxTight",100,0,100);
  hFillNVtxIso = fs->make<TH1F>("hFillNVtxIso","hFillNVtxIso",100,0,100);
  hFillNVtxTrkIso = fs->make<TH1F>("hFillNVtxTrkIso","hFillNVtxTrkIso",100,0,100);
  hFillNVtxTightIso = fs->make<TH1F>("hFillNVtxTightIso","hFillNVtxTightIso",100,0,100);
  hFillNVtxTightTrkIso = fs->make<TH1F>("hFillNVtxTightTrkIso","hFillNVtxTightTrkIso",100,0,100);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise4::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonExercise4::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

// soft muon new
bool MuonExercise4::SoftMuon(const pat::Muon& mu , const reco::Vertex& vertex) {

  if (!(muon::isGoodMuon(mu, muon::TMOneStationTight))) return false;
  if (!(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)) return false;
  if (!(mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0)) return false;
  if (!(mu.innerTrack()->quality(reco::TrackBase::highPurity))) return false;
  if (!((fabs(mu.innerTrack()->dxy(vertex.position())) < 0.3) && (fabs(mu.innerTrack()->dz(vertex.position())) < 20.))) return false;
  return true;

}


bool MuonExercise4::MatchedToTriggerObject(const pat::Muon& mu, const edm::Event& ev, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<edm::TriggerResults> triggerBits, int & index ) {

  bool trigflag = false;
  index = -1;
  double dR = 0.1;
  const edm::TriggerNames &triggerNames = ev.triggerNames(*triggerBits);
  for(pat::TriggerObjectStandAlone trig : *triggerObjects){
    trig.unpackPathNames(triggerNames);
    if ( !(trig.hasTriggerObjectType(trigger::TriggerMuon))) continue;
    double dRtmp = deltaR(trig, *(mu.innerTrack())); 
    if(dRtmp < dR){
      //  dR = dRtmp;
      //   index = std::distance(triggerObjects->begin(), &trig);
      trigflag = true;
    }
  }

  return trigflag;

}


bool MuonExercise4::MatchedToGenObject(const pat::Muon& mu, Handle<edm::View<pat::PackedGenParticle> > packed, int & genindex ) {
  
  bool genflag = false;
  genindex = -1;
  double dRgen = 0.1;
  for(edm::View<pat::PackedGenParticle>::const_iterator genmuon = packed->begin(); genmuon != packed->end(); ++genmuon){
    if(!(fabs(genmuon->pdgId()) == 13)) continue;
    
    if(genmuon->numberOfMothers() > 0) {
      const reco::GenParticle* firstMother = (const reco::GenParticle*)genmuon->mother(0);
      
      if (firstMother != 0) {
	if (firstMother->pdgId() != genmuon->pdgId()){ mother = firstMother; }
	else{
	  const reco::GenParticle newMother = getMother(*firstMother);
	  mother = &newMother;
	}
      }
    }
    
    double dRtmp1 = deltaR(*genmuon, *(mu.innerTrack()));
    if(mother->pdgId() == 23 && dRtmp1 < dRgen ){
      dRgen = dRtmp1;
      genindex = std::distance(packed->begin(), genmuon);
      genflag = true;
    } 
  }
  
  return genflag;

}	  


const reco::GenParticle MuonExercise4::getMother(const reco::GenParticle& gp)
{

  const reco::GenParticle* mom = &gp;
  while (mom->numberOfMothers() > 0)
    {
      for (unsigned int idx = 0; idx < mom->numberOfMothers(); idx++)
	{
	  mom = dynamic_cast<const reco::GenParticle*>(mom->mother(idx));
	  if (mom->pdgId() != gp.pdgId())
	    return (*mom);
	}
    }

  return (*mom);

}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise4);
