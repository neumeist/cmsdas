// -*- C++ -*-
//
// Package:    CMSDASExercises/MuonExercise2
// Class:      MuonExercise2
// 
/**\class MuonExercise3 MuonExercise2.cc CMSDASExercises/MuonExercise3/plugins/MuonExercise2.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Norbert Neumeister
//         Created:  Thu, 10 Dec 2016 21:10:01 GMT
//
//

// system include files
#include <memory>
#include <iomanip>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TProfile.h>
#include <TTree.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class MuonExercise2 : public edm::EDAnalyzer {

   public:

      explicit MuonExercise2(const edm::ParameterSet&);
      ~MuonExercise2();

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
      
      edm::EDGetTokenT<pat::MuonCollection> muonCollToken;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> genCollToken;
  
      TH1F* h_RecDiMuonM;
      TH1F* h_GenDiMuonM;
      TH1F* h_MassRes;
      TH1F* h_MupTRes;
  
      // ProfileHistograms declaration for scale estimation
      TProfile* prof_MuPlusPhivsDiMuonM;
      TProfile* prof_MuMinusPhivsDiMuonM;
      TProfile* prof_MuEtavsDiMuonM;
  
      // ProfileHistograms declaration for resolution study
      TProfile* prof_MuEtavspTRes;
      TProfile* prof_MupTvspTRes;
  
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
MuonExercise2::MuonExercise2(const edm::ParameterSet& iConfig)
{

  edm::InputTag theMuonLabel("slimmedMuons");
  edm::InputTag theGenMuonLabel("packedGenParticles");
  
  muonCollToken = consumes<pat::MuonCollection>(theMuonLabel);
  genCollToken = consumes<pat::PackedGenParticleCollection>(theGenMuonLabel);

  edm::Service<TFileService> fs;
  
  h_RecDiMuonM = fs->make<TH1F>("h_RecDiMuonM",";m_{#mu^{+}#mu^{-}};",80,70,110);
  h_GenDiMuonM = fs->make<TH1F>("h_GenDiMuonM",";m_{#mu^{+}#mu^{-}};",80,70,110);
  h_MassRes = fs->make<TH1F>("h_MassRes","Mass Resolution",80,-0.15,0.15);
  h_MupTRes = fs->make<TH1F>("h_MupTRes","Muon p_{T} resolution;#Delta p_{T}/p_{T};",80,-0.2,0.2);
  
  // ProfileHistograms declaration for scale estimation
  prof_MuPlusPhivsDiMuonM = fs->make<TProfile>("prof_MuPlusPhivsDiMuonM","#mu^{+} #phi vs m_{#mu^{+}#mu^{-}};Reco muon(+) #phi[rad]; Z peak position [GeV/c^2]",16,-3.14,3.14,88,93);
  prof_MuMinusPhivsDiMuonM = fs->make<TProfile>("prof_MuMinusPhivsDiMuonM","#mu^{-} #phi vs m_{#mu^{+}#mu^{-}};Reco muon(-) #phi[rad];Z peak position [GeV/c^2]",16,-3.14,3.14,88,93);
  prof_MuEtavsDiMuonM = fs->make<TProfile>("prof_MuEtavsDiMuonM","Muon #eta vs m_{#mu^{+}#mu^{-}};Reco Muon #eta; Z peak position [GeV/c^2]",50,-2.4,2.4,88,93);
  
  // ProfileHistograms declaration for resolution study
  prof_MuEtavspTRes = fs->make<TProfile>("prof_MuEtavspTRes",";Gen Muon #eta;#Delta p_{T}/p_{T}",50,-2.4,2.4,0,1);
  prof_MupTvspTRes = fs->make<TProfile>("prof_MupTvspTRes",";Gen Muon p_{T} [GeV];#Delta p_{T}/p_{T}",25,20,100,0,1);

}


MuonExercise2::~MuonExercise2()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonExercise2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
   
  edm::Handle<vector<pat::Muon>> muons;
  iEvent.getByToken(muonCollToken, muons);
   
  edm::Handle <pat::PackedGenParticleCollection> genColl;
  iEvent.getByToken(genCollToken, genColl);
   
  // find mu+
  for (auto mup = muons->cbegin(); mup != muons->cend(); ++mup) {
    if ( not mup->isGlobalMuon() ) continue;
    if ( not mup->charge() > 0 ) continue;
    if ( not mup->pt() > 20.0 ) continue; 
    if ( fabs(mup->eta()) > 2.4 ) continue; 
    if ( not mup->chargedHadronIso() < 0.15 ) continue;
    if ( fabs(mup->eta()) > 2.4 ) continue; 
      
    // find mu-
    for (auto mum = mup+1; mum != muons->cend(); ++mum) {
      if ( not mum->isGlobalMuon() ) continue;
      if ( not mum->charge() < 0 ) continue;
      if ( not mum->pt() > 20.0 ) continue; 
      if ( fabs(mum->eta()) > 2.4 ) continue; 
      if ( not mum->chargedHadronIso() < 0.15 ) continue;

      // dimuon invariant mass 
      double mRec = (mup->p4() + mum->p4()).M();
      if ( mRec < 70 || mRec > 110) continue;
      h_RecDiMuonM->Fill(mRec);
       
      // fill TProfile Histograms
      prof_MuPlusPhivsDiMuonM->Fill(mup->phi(),mRec,1);
      prof_MuMinusPhivsDiMuonM->Fill(mum->phi(),mRec,1);
      prof_MuEtavsDiMuonM->Fill(mup->eta(),mRec,1);
       
      // Gen Matching       
      int idxmup_G = -1;
      int idxmum_G = -1;   

      for (auto mu_G = genColl->cbegin(); mu_G != genColl->cend(); ++mu_G) {
        const pat::PackedGenParticle& mcMuon = (*mu_G);
        if ( not mcMuon.pdgId() == -13 ) continue;
        if ( fabs(mcMuon.eta()) > 2.4 ) continue;
        if ( not mcMuon.pt() > 1.5 ) continue;	 
        if ( deltaR(*mu_G, *(mup->innerTrack())) < 0.1 && mcMuon.charge() > 0 ) idxmup_G = std::distance(genColl->cbegin(), mu_G);
        if ( deltaR(*mu_G, *(mum->innerTrack())) < 0.1 && mcMuon.charge() < 0 ) idxmum_G = std::distance(genColl->cbegin(), mu_G);
      }

      if ( idxmup_G > -1 && idxmum_G > -1) {
        double mGen = (genColl->at(idxmup_G).p4() + genColl->at(idxmum_G).p4()).M();
        h_GenDiMuonM->Fill(mGen);
        h_MassRes->Fill((mRec - mGen)/mGen);
        h_MupTRes->Fill((1/mup->pt()-1/genColl->at(idxmup_G).pt())/(1/genColl->at(idxmup_G).pt()));
        h_MupTRes->Fill((1/mum->pt()-1/genColl->at(idxmum_G).pt())/(1/genColl->at(idxmum_G).pt()));
       
        // ProfileHistograms declaration for resolution study
        prof_MuEtavspTRes->Fill(genColl->at(idxmup_G).eta(),(1/mup->pt()-1/genColl->at(idxmup_G).pt())/(1/genColl->at(idxmup_G).pt()));
        prof_MuEtavspTRes->Fill(genColl->at(idxmum_G).eta(),(1/mum->pt()-1/genColl->at(idxmum_G).pt())/(1/genColl->at(idxmum_G).pt()));
       
        prof_MupTvspTRes->Fill(genColl->at(idxmup_G).pt(),(1/mup->pt()-1/genColl->at(idxmup_G).pt())/(1/genColl->at(idxmup_G).pt()));
        prof_MupTvspTRes->Fill(genColl->at(idxmum_G).pt(),(1/mum->pt()-1/genColl->at(idxmum_G).pt())/(1/genColl->at(idxmum_G).pt()));
      }
   }
 }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonExercise2::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise2::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MuonExercise2::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MuonExercise2::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MuonExercise2::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MuonExercise2::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonExercise2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise2);
