#define tth_analysis_cxx
#include "tth_analysis.h"
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TSystem.h>

#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCRelation.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "IOIMPL/LCCollectionIOVec.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;

std::vector<EVENT::MCParticle*> getMCTaus(EVENT::LCEvent* event){
  std::vector<EVENT::MCParticle*> MC_tau_vector;
  EVENT::LCCollection* mcp = event->getCollection("MCParticle");

  for (int i = 0; i < mcp->getNumberOfElements(); i++){
    EVENT::MCParticle* mcpi = (EVENT::MCParticle*) mcp->getElementAt(i);
    if (abs(mcpi->getPDG()) == 15){

      //cout << "->got an MC tau " << mcpi->getPDG() << " energy " << mcpi->getEnergy() << " with daughters: " << mcpi->getDaughters().size() << " generator status " << mcpi->getGeneratorStatus() << endl;
      if (mcpi->getParents().size() > 0){
	//cout << "tau from mother: " << mcpi->getParents()[0]->getPDG() << endl;
	if (mcpi->getParents()[0]->getParents().size() > 0){
	  //cout << "tau mother from mother: " << mcpi->getParents()[0]->getParents()[0]->getPDG() << endl;

	  //if (abs(mcpi->getParents()[0]->getPDG()) == 24 && abs(mcpi->getParents()[0]->getParents()[0]->getPDG()) == 94){
	  if (abs(mcpi->getParents()[0]->getPDG()) == 24 && (abs(mcpi->getParents()[0]->getParents()[0]->getPDG()) == 94 || abs(mcpi->getParents()[0]->getParents()[0]->getPDG()) == 15 || abs(mcpi->getParents()[0]->getParents()[0]->getPDG()) == 16)){
	    if (mcpi->getDaughters().size() > 0){
	      MC_tau_vector.push_back(mcpi);
	    }
	  }
	}
      }
    }
  }
  return MC_tau_vector;
}

std::vector<EVENT::MCParticle*> getMCTauDaughters(EVENT::MCParticle* mcpi){
  std::vector<EVENT::MCParticle*> MC_tau_daughter_vector;

  for (unsigned int k = 0; k < mcpi->getDaughters().size(); k++){
    EVENT::MCParticle* daughterk = (EVENT::MCParticle*) mcpi->getDaughters()[k];
    if (daughterk->getGeneratorStatus() == 1){
      MC_tau_daughter_vector.push_back(daughterk);
    }
    if (daughterk->getGeneratorStatus() == 2){
      for(unsigned int m = 0; m < daughterk->getDaughters().size(); m++){
	EVENT::MCParticle* daughterm = (EVENT::MCParticle*) daughterk->getDaughters()[m];
	if (daughterm->getGeneratorStatus() == 1){
	  MC_tau_daughter_vector.push_back(daughterm);
	}
	if (daughterm->getGeneratorStatus() == 2){
	  for(unsigned int n = 0; n < daughterm->getDaughters().size(); n++){
	    EVENT::MCParticle* daughtern = (EVENT::MCParticle*) daughterm->getDaughters()[n];
	    if (daughtern->getGeneratorStatus() == 1){
	      MC_tau_daughter_vector.push_back(daughtern);
	    }
	  }
	}
      }
    }
  }

  return MC_tau_daughter_vector;
}

std::vector<EVENT::MCParticle*> getAllMCTauDaughters(EVENT::MCParticle* mcpi, std::vector<EVENT::MCParticle*> daughters){

  for (unsigned int k = 0; k < mcpi->getDaughters().size(); k++){
    EVENT::MCParticle* daughterk = (EVENT::MCParticle*) mcpi->getDaughters()[k];
    if (daughterk->getGeneratorStatus() == 1){
      daughters.push_back(daughterk);
    }
    if (daughterk->getGeneratorStatus() == 2){
      for(unsigned int m = 0; m < daughterk->getDaughters().size(); m++){
	EVENT::MCParticle* daughterm = (EVENT::MCParticle*) daughterk->getDaughters()[m];
	if (daughterm->getGeneratorStatus() == 1){
	  daughters.push_back(daughterm);
	}
	if (daughterm->getGeneratorStatus() == 2){
	  for(unsigned int n = 0; n < daughterm->getDaughters().size(); n++){
	    EVENT::MCParticle* daughtern = (EVENT::MCParticle*) daughterm->getDaughters()[n];
	    if (daughtern->getGeneratorStatus() == 1){
	      daughters.push_back(daughtern);
	    }
	  }
	}
      }
    }
  }
  return daughters;
}

std::vector<EVENT::ReconstructedParticle*> getRecoTaus(EVENT::LCEvent* event){
  std::vector<EVENT::ReconstructedParticle*> Reco_tau_vector;
  EVENT::LCCollection* taus = event->getCollection("Taus");

    for(int i = 0; i < taus->getNumberOfElements(); i++){
      EVENT::ReconstructedParticle* taui = (EVENT::ReconstructedParticle*) taus->getElementAt(i);
      Reco_tau_vector.push_back(taui);
    }

  return Reco_tau_vector;
}

std::vector<EVENT::ReconstructedParticle*> getRecoTauDaughters(EVENT::ReconstructedParticle* taui){
  std::vector<EVENT::ReconstructedParticle*> Reco_tau_daughter_vector;

  for(unsigned int j = 0; j < taui->getParticles().size(); j++){
    EVENT::ReconstructedParticle* tau_dj = (EVENT::ReconstructedParticle*) taui->getParticles()[j];
    Reco_tau_daughter_vector.push_back(tau_dj);
  }

  return Reco_tau_daughter_vector;
}

std::vector<EVENT::ReconstructedParticle*> getAllRecoTauDaughters(EVENT::ReconstructedParticle* taui, std::vector<EVENT::ReconstructedParticle*> daughters){
  
  for(unsigned int j = 0; j < taui->getParticles().size(); j++){
    EVENT::ReconstructedParticle* tau_dj = (EVENT::ReconstructedParticle*) taui->getParticles()[j];
    daughters.push_back(tau_dj);
  }

  return daughters;
}

void postLeptonFinder(std::string recfileName, TString outfileName){
  // to see how well the marlin processor did
  IO::LCReader* reclcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
  reclcReader->open( recfileName ) ;

  TFile outfile( outfileName,"RECREATE");

  TTree *Event_tree = new TTree("Event_tree","");
  Int_t N_MC_Taus, N_Found_Taus, N_Found_Isolep, GenLeptonPDG, MCLeptonPDG, RecoLeptonType;
  Event_tree->Branch("N_MC_Taus",&N_MC_Taus,"N_MC_Taus/I");
  Event_tree->Branch("N_Found_Taus",&N_Found_Taus,"N_Found_Taus/I");
  Event_tree->Branch("N_Found_Isolep",&N_Found_Isolep,"N_Found_Isolep/I");
  Event_tree->Branch("GenLeptonPDG",&GenLeptonPDG,"GenLeptonPDG/I");
  Event_tree->Branch("MCLeptonPDG",&MCLeptonPDG,"MCLeptonPDG/I");
  Event_tree->Branch("RecoLeptonType",&RecoLeptonType,"RecoLeptonType/I");

  TTree *MCIsolep_tree = new TTree("MCIsolep_tree","");
  Int_t MCIsolep_Found, MCIsolep_PDG;
  MCIsolep_tree->Branch("MCIsolep_Found",&MCIsolep_Found,"MCIsolep_Found/I");
  MCIsolep_tree->Branch("MCIsolep_PDG",&MCIsolep_PDG,"MCIsolep_PDG/I");

  TTree *RecoIsolep_tree = new TTree("RecoIsolep_tree","");
  Int_t RecoIsolep_Fake, RecoIsolep_Type;
  Double_t RecoIsolep_Energy, RecoIsolep_Pt, RecoIsolep_CosTheta, RecoIsolep_D0, RecoIsolep_Z0, RecoIsolep_R0, RecoIsolep_CalE;
  RecoIsolep_tree->Branch("RecoIsolep_Fake",&RecoIsolep_Fake,"RecoIsolep_Fake/I");
  RecoIsolep_tree->Branch("RecoIsolep_Type",&RecoIsolep_Type,"RecoIsolep_Type/I");
  RecoIsolep_tree->Branch("RecoIsolep_Energy",&RecoIsolep_Energy,"RecoIsolep_Energy/D");
  RecoIsolep_tree->Branch("RecoIsolep_Pt",&RecoIsolep_Pt,"RecoIsolep_Pt/D");
  RecoIsolep_tree->Branch("RecoIsolep_CosTheta",&RecoIsolep_CosTheta,"RecoIsolep_CosTheta/D");
  RecoIsolep_tree->Branch("RecoIsolep_D0",&RecoIsolep_D0,"RecoIsolep_D0/D");
  RecoIsolep_tree->Branch("RecoIsolep_Z0",&RecoIsolep_Z0,"RecoIsolep_Z0/D");
  RecoIsolep_tree->Branch("RecoIsolep_R0",&RecoIsolep_R0,"RecoIsolep_R0/D");
  RecoIsolep_tree->Branch("RecoIsolep_CalE",&RecoIsolep_CalE,"RecoIsolep_CalE/D");

  TTree *MCTau_tree = new TTree("MCTau_tree","");
  Int_t MCTau_Found, MCTau_NNT, MCTau_NCT, MCTau_NDaughters,MCTau_NRingParticles;
  Double_t  MCTau_MaxDAngle, MCTau_MaxDPt, MCTau_Mass,MCTau_RingEnergy;
  MCTau_tree->Branch("MCTau_Found",&MCTau_Found,"MCTau_Found/I");
  MCTau_tree->Branch("MCTau_NNT",&MCTau_NNT,"MCTau_NNT/I");
  MCTau_tree->Branch("MCTau_NCT",&MCTau_NCT,"MCTau_NCT/I");
  MCTau_tree->Branch("MCTau_NDaughters",&MCTau_NDaughters,"MCTau_NDaughters/I");
  MCTau_tree->Branch("MCTau_NRingParticles",&MCTau_NRingParticles,"MCTau_NRingParticles/I");
  MCTau_tree->Branch("MCTau_MaxDAngle",&MCTau_MaxDAngle,"MCTau_MaxDAngle/D");
  MCTau_tree->Branch("MCTau_MaxDPt",&MCTau_MaxDPt,"MCTau_MaxDPt/D");
  MCTau_tree->Branch("MCTau_Mass",&MCTau_Mass,"MCTau_Mass/D");
  MCTau_tree->Branch("MCTau_RingEnergy",&MCTau_RingEnergy,"MCTau_RingEnergy/D");

  TTree *MCTauDaughter_tree = new TTree("MCTauDaughter_tree","");
  Double_t MCTauDaughter_Pt, MCTauDaughter_CosTheta, MCTauDaughter_TauAngle;
  Int_t MCTau_Found2;
  MCTauDaughter_tree->Branch("MCTauDaughter_Pt",&MCTauDaughter_Pt,"MCTauDaughter_Pt/D");
  MCTauDaughter_tree->Branch("MCTauDaughter_CosTheta",&MCTauDaughter_CosTheta,"MCTauDaughter_CosTheta/D");
  MCTauDaughter_tree->Branch("MCTauDaughter_TauAngle",&MCTauDaughter_TauAngle,"MCTauDaughter_TauAngle/D");
  MCTauDaughter_tree->Branch("MCTau_Found2",&MCTau_Found2,"MCTau_Found2/I");

  TTree *RecoTau_tree = new TTree("RecoTau_tree","");
  Int_t RecoTau_Fake, RecoTau_NRingParticles, RecoTau_NDaughters, RecoTau_NCT, RecoTau_NNT;
  Double_t RecoTau_RingEnergy, RecoTau_Mass, RecoTau_MaxDPt, RecoTau_MaxDAngle, RecoTau_MinPFOAngle;
  RecoTau_tree->Branch("RecoTau_Fake",&RecoTau_Fake,"RecoTau_Fake/I");
  RecoTau_tree->Branch("RecoTau_NRingParticles",&RecoTau_NRingParticles,"RecoTau_NRingParticles/I");
  RecoTau_tree->Branch("RecoTau_NDaughters",&RecoTau_NDaughters,"RecoTau_NDaughters/I"); 
  RecoTau_tree->Branch("RecoTau_NCT",&RecoTau_NCT,"RecoTau_NCT/I");
  RecoTau_tree->Branch("RecoTau_NNT",&RecoTau_NNT,"RecoTau_NNT/I");
  RecoTau_tree->Branch("RecoTau_RingEnergy",&RecoTau_RingEnergy,"RecoTau_RingEnergy/D");
  RecoTau_tree->Branch("RecoTau_Mass",&RecoTau_Mass,"RecoTau_Mass/D");
  RecoTau_tree->Branch("RecoTau_MaxDPt",&RecoTau_MaxDPt,"RecoTau_MaxDPt/D");
  RecoTau_tree->Branch("RecoTau_MaxDAngle",&RecoTau_MaxDAngle,"RecoTau_MaxDAngle/D");
  RecoTau_tree->Branch("RecoTau_MinPFOAngle",&RecoTau_MinPFOAngle,"RecoTau_MinPFOAngle/D");

  TTree *RecoTauDaughter_tree = new TTree("RecoTauDaughter_tree","");
  Double_t RecoTauDaughter_Pt, RecoTauDaughter_CosTheta, RecoTauDaughter_TauAngle;
  Double_t RecoTauDaughter_D0, RecoTauDaughter_Z0, RecoTauDaughter_R0, RecoTauDaughter_CalE;
  Int_t RecoTau_Fake2;
  RecoTauDaughter_tree->Branch("RecoTauDaughter_Pt",&RecoTauDaughter_Pt,"RecoTauDaughter_Pt/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_CosTheta",&RecoTauDaughter_CosTheta,"RecoTauDaughter_CosTheta/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_TauAngle",&RecoTauDaughter_TauAngle,"RecoTauDaughter_TauAngle/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_D0",&RecoTauDaughter_D0,"RecoTauDaughter_D0/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_Z0",&RecoTauDaughter_Z0,"RecoTauDaughter_Z0/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_R0",&RecoTauDaughter_R0,"RecoTauDaughter_R0/D");
  RecoTauDaughter_tree->Branch("RecoTau_Fake2",&RecoTau_Fake2,"RecoTau_Fake2/I");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_CalE",&RecoTauDaughter_CalE,"RecoTauDaughter_CalE/D");

  EVENT::LCEvent* recevent = 0;

  int total_mc_taus = 0;
  int total_true_taus = 0;
  int total_fake_taus = 0;
  int total_taus_MCEMU = 0; // total found taus which were e/mus

  int total_mc_emus = 0;
  int total_true_emus = 0;
  int total_fake_emus = 0;
  int total_emus_MCTD = 0; // total found isoleps which were tau daughters

  while( (recevent = reclcReader->readNextEvent()) != 0){

    //cout << "=========================" << endl;

    EVENT::LCCollection* links = recevent->getCollection("RecoMCTruthLink");
    EVENT::LCCollection* pfo = recevent->getCollection("SelectedPandoraPFOCollection");

    EVENT::MCParticle* gen_lepton = getGeneratedLepton(recevent);
    //cout << "generated lepton pdg: " << gen_lepton->getPDG() << endl; // 11, 13 or 15
    GenLeptonPDG = gen_lepton->getPDG();
    
    EVENT::MCParticle* mc_lepton = getMCLepton(recevent); // 11 or 13 only // not anymore!
    if (mc_lepton==0){
      //cout << "couldn't find generated lepton in mc particles" << endl;
      MCLeptonPDG = 0;
    }
    else{
      //cout << "MC lepton pdg: " << mc_lepton->getPDG() << endl;
      MCLeptonPDG = mc_lepton->getPDG();
    }

    EVENT::ReconstructedParticle* recod_lepton = getReconstructedLepton(recevent); // 11 or 13 only
    if (recod_lepton==0){
      //cout << "couldn't find generated lepton in reconstructed particles" << endl;
      RecoLeptonType = 0;
    }
    else{
      //cout << "found generated lepton in reconstructed particles" << endl;
      RecoLeptonType = recod_lepton->getType();
      total_mc_emus++;
    }
    
    // moving up
    // get the MC taus which have daughters
    std::vector<EVENT::MCParticle*> MCTaus;
    MCTaus.clear();
    MCTaus = getMCTaus(recevent);
    //cout << "found " << MCTaus.size() << " MC taus which have daughters" << endl;
    N_MC_Taus = MCTaus.size();
    total_mc_taus = total_mc_taus + N_MC_Taus;

    // get every MC tau daughter in the event
    std::vector<EVENT::MCParticle*> AllMCTauDaughters;
    AllMCTauDaughters.clear();
    for(unsigned int p = 0; p < MCTaus.size(); p++){
      EVENT::MCParticle* MCTau = MCTaus[p];
      AllMCTauDaughters = getAllMCTauDaughters(MCTau, AllMCTauDaughters);
    }
    //cout << "got a total of " << AllMCTauDaughters.size() << " MC tau daughter particles in event" << endl;

    // get the reconstructed taus
    std::vector<EVENT::ReconstructedParticle*> RecoTaus;
    RecoTaus.clear();
    RecoTaus = getRecoTaus(recevent);
    //cout << "found " << RecoTaus.size() << " reconstructed taus" << endl;
    N_Found_Taus = RecoTaus.size();

    // get every recontructed tau daughter in the event
    std::vector<EVENT::ReconstructedParticle*> AllRecoTauDaughters;
    AllRecoTauDaughters.clear();
    for(unsigned int p = 0; p < RecoTaus.size(); p++){
      EVENT::ReconstructedParticle* RecoTau = RecoTaus[p];
      AllRecoTauDaughters = getAllRecoTauDaughters(RecoTau, AllRecoTauDaughters);
    }
    //cout << "got a total of " << AllRecoTauDaughters.size() << " reco tau daughter particles in event" << endl;
    // end of moving up


    // get the found isolated leptons
    EVENT::LCCollection* found_leptons = recevent->getCollection("Isolep_Selected");
    //cout << "found " << found_leptons->getNumberOfElements() << " isolated letpons" << endl;
    N_Found_Isolep = found_leptons->getNumberOfElements();
    for (int i = 0; i < found_leptons->getNumberOfElements(); i++){
      EVENT::ReconstructedParticle* leptoni = (EVENT::ReconstructedParticle*) found_leptons->getElementAt(i);
      TLorentzVector LFV = TLorentzVector(leptoni->getMomentum()[0],leptoni->getMomentum()[1],leptoni->getMomentum()[2],leptoni->getEnergy());

      // is this lepton true or fake?
      RecoIsolep_Fake = 1;
      if (leptoni == recod_lepton){
	RecoIsolep_Fake = 0;
	total_true_emus++;
	//cout << "this lepton matches the reconstructed signal lepton" << endl;
      }
      else{
	total_fake_emus++;
      }

      // is this isolated lepton a tau daughter?
      for(int r = 0; r < links->getNumberOfElements(); r++){
	EVENT::LCRelation* linkr = (EVENT::LCRelation*) links->getElementAt(r);
	EVENT::ReconstructedParticle* rcr = (EVENT::ReconstructedParticle*) linkr->getFrom();
	if (rcr == leptoni){
	  EVENT::MCParticle* mcpr = (EVENT::MCParticle*) linkr->getTo();
	  if(std::find(AllMCTauDaughters.begin(), AllMCTauDaughters.end(), mcpr) != AllMCTauDaughters.end()){
	    total_emus_MCTD++;
	  }
	}
      }

      const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) leptoni->getTracks();
      if (trkvec.size()>0){
	RecoIsolep_D0 = abs(trkvec[0]->getD0());
	RecoIsolep_Z0 = abs(trkvec[0]->getZ0());
	RecoIsolep_R0 = sqrt( RecoIsolep_D0*RecoIsolep_D0 + RecoIsolep_Z0*RecoIsolep_Z0 );
      }
      else{
	RecoIsolep_D0 = 99;
	RecoIsolep_Z0 = 99;
	RecoIsolep_R0 = 99;
      }

      std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) leptoni->getClusters();
      float ecal = 0;
      float hcal = 0;
      for ( std::vector<EVENT::Cluster*>::const_iterator iCluster=clusters.begin(); iCluster!=clusters.end(); ++iCluster) {
	ecal += (*iCluster)->getSubdetectorEnergies()[0];
	hcal += (*iCluster)->getSubdetectorEnergies()[1];
      }
      RecoIsolep_CalE = ecal/(ecal+hcal);

      RecoIsolep_Energy = leptoni->getEnergy();
      RecoIsolep_Pt = LFV.Pt();
      RecoIsolep_CosTheta = cos(LFV.Theta());
      RecoIsolep_Type = leptoni->getType();

      RecoIsolep_tree->Fill();
    }

    // was the generated lepton reconstructed and found? (e/mu events only)
    if (abs(gen_lepton->getPDG())==11 || abs(gen_lepton->getPDG())==13){

      if (mc_lepton==0){
	MCIsolep_PDG = 0;
      }
      else{
	MCIsolep_PDG = mc_lepton->getPDG();
      } 

      MCIsolep_Found = 0;
      for(int r = 0; r < links->getNumberOfElements(); r++){
	EVENT::LCRelation* linkr = (EVENT::LCRelation*) links->getElementAt(r);
	EVENT::MCParticle* mcpr = (EVENT::MCParticle*) linkr->getTo();
	if (mcpr == mc_lepton){
	  EVENT::ReconstructedParticle* rcr = (EVENT::ReconstructedParticle*) linkr->getFrom();
	  // was this particle an isolep?
	  for (int i = 0; i < found_leptons->getNumberOfElements(); i++){
	    EVENT::ReconstructedParticle* leptoni = (EVENT::ReconstructedParticle*) found_leptons->getElementAt(i);
	    if(rcr == leptoni){
	      //cout << "generated lepton was reconstructed and identified as an isolep" << endl;
	      MCIsolep_Found = 1;
	    }
	  }
	}
      }
      MCIsolep_tree->Fill();
    }

    // moving up
    // get the MC taus which have daughters
    // std::vector<EVENT::MCParticle*> MCTaus;
    // MCTaus.clear();
    // MCTaus = getMCTaus(recevent);
    // //cout << "found " << MCTaus.size() << " MC taus which have daughters" << endl;
    // N_MC_Taus = MCTaus.size();
    // total_mc_taus = total_mc_taus + N_MC_Taus;

    // // get every MC tau daughter in the event
    // std::vector<EVENT::MCParticle*> AllMCTauDaughters;
    // AllMCTauDaughters.clear();
    // for(unsigned int p = 0; p < MCTaus.size(); p++){
    //   EVENT::MCParticle* MCTau = MCTaus[p];
    //   AllMCTauDaughters = getAllMCTauDaughters(MCTau, AllMCTauDaughters);
    // }
    //cout << "got a total of " << AllMCTauDaughters.size() << " MC tau daughter particles in event" << endl;

    // get the reconstructed taus
    // std::vector<EVENT::ReconstructedParticle*> RecoTaus;
    // RecoTaus.clear();
    // RecoTaus = getRecoTaus(recevent);
    // //cout << "found " << RecoTaus.size() << " reconstructed taus" << endl;
    // N_Found_Taus = RecoTaus.size();

    // // get every recontructed tau daughter in the event
    // std::vector<EVENT::ReconstructedParticle*> AllRecoTauDaughters;
    // AllRecoTauDaughters.clear();
    // for(unsigned int p = 0; p < RecoTaus.size(); p++){
    //   EVENT::ReconstructedParticle* RecoTau = RecoTaus[p];
    //   AllRecoTauDaughters = getAllRecoTauDaughters(RecoTau, AllRecoTauDaughters);
    // }
    //cout << "got a total of " << AllRecoTauDaughters.size() << " reco tau daughter particles in event" << endl;
    // end of moving up

    // Properties of MC taus
    for(unsigned int p = 0; p < MCTaus.size(); p++){
      EVENT::MCParticle* MCTau = MCTaus[p];
      TLorentzVector MCTauFV = TLorentzVector(MCTau->getMomentum()[0],MCTau->getMomentum()[1],MCTau->getMomentum()[2],MCTau->getEnergy());

      std::vector<EVENT::MCParticle*> MCTauDaughters;
      MCTauDaughters.clear();
      MCTauDaughters = getMCTauDaughters(MCTau);
      
      MCTau_Found = 0;
      MCTau_NNT = 0;
      MCTau_NCT = 0;
      MCTau_MaxDAngle = 0;
      MCTau_MaxDPt = 0;
      MCTau_Mass = MCTauFV.M();
      MCTau_NDaughters = MCTauDaughters.size();

      for (unsigned int q = 0; q < MCTauDaughters.size(); q++){
	EVENT::MCParticle* daughter_track = MCTauDaughters[q];
	TLorentzVector DFV = TLorentzVector(daughter_track->getMomentum()[0],daughter_track->getMomentum()[1],daughter_track->getMomentum()[2],daughter_track->getEnergy());

	// find the reco particle linked to this MC particle - was this daughter in any found tau?
       	for(int r = 0; r < links->getNumberOfElements(); r++){
	  EVENT::LCRelation* linkr = (EVENT::LCRelation*) links->getElementAt(r);
	  EVENT::MCParticle* mcpr = (EVENT::MCParticle*) linkr->getTo();
	  if (mcpr == daughter_track){
	    EVENT::ReconstructedParticle* rcr = (EVENT::ReconstructedParticle*) linkr->getFrom();
	    if(std::find(AllRecoTauDaughters.begin(), AllRecoTauDaughters.end(), rcr) != AllRecoTauDaughters.end()){
	      //cout << "this MC daughter was found in a tau" << endl;
	      MCTau_Found = 1;
	    }
	  }
	}

	if (DFV.Pt() > MCTau_MaxDPt){
	  MCTau_MaxDPt = DFV.Pt();
	}
	if (DFV.Angle(MCTauFV.Vect()) > MCTau_MaxDAngle){
	  MCTau_MaxDAngle = DFV.Angle(MCTauFV.Vect());
	}
	if (daughter_track->getCharge() == 0){
	  MCTau_NNT++;
	}
	else{
	  MCTau_NCT++;
	}
      }

      // loop for daughters again, now we know if mc tau is found or not
      for (unsigned int q = 0; q < MCTauDaughters.size(); q++){
	EVENT::MCParticle* daughter_track = MCTauDaughters[q];
	TLorentzVector DFV = TLorentzVector(daughter_track->getMomentum()[0],daughter_track->getMomentum()[1],daughter_track->getMomentum()[2],daughter_track->getEnergy());

	MCTauDaughter_Pt = DFV.Pt();
	MCTauDaughter_CosTheta = cos(DFV.Theta());
	MCTauDaughter_TauAngle = DFV.Angle(MCTauFV.Vect());
	MCTau_Found2 = MCTau_Found;	

	MCTauDaughter_tree->Fill();
      } 

      // energy within isolation ring about MCTau 
      MCTau_RingEnergy = 0;
      MCTau_NRingParticles = 0;
      for (int j = 0; j < pfo->getNumberOfElements(); j++) {
	EVENT::ReconstructedParticle* pfoj = (EVENT::ReconstructedParticle*) pfo->getElementAt(j);
	TLorentzVector pfoj_vec = TLorentzVector(pfoj->getMomentum()[0],pfoj->getMomentum()[1],pfoj->getMomentum()[2],pfoj->getEnergy());

	double pfoi_pfoj_angle = pfoj_vec.Angle(MCTauFV.Vect());
	if (abs(pfoi_pfoj_angle) > 0.04 && abs(pfoi_pfoj_angle) < 0.29){
	  MCTau_RingEnergy = MCTau_RingEnergy + pfoj->getEnergy();
	  MCTau_NRingParticles++;
	}
      }
      
      MCTau_tree->Fill();
    } // end of MC tau loop

    // Properties of Reco taus
    for(unsigned int p = 0; p < RecoTaus.size(); p++){
      EVENT::ReconstructedParticle* RecoTau = RecoTaus[p];
      TLorentzVector RecoTauFV = TLorentzVector(RecoTau->getMomentum()[0],RecoTau->getMomentum()[1],RecoTau->getMomentum()[2],RecoTau->getEnergy());
      //cout << "this tau FV: " << RecoTau->getMomentum()[0] << " " << RecoTau->getMomentum()[1] << " " << RecoTau->getMomentum()[2] << " " << RecoTau->getEnergy() << endl;

      std::vector<EVENT::ReconstructedParticle*> RecoTauDaughters;
      RecoTauDaughters.clear();
      RecoTauDaughters = getRecoTauDaughters(RecoTau);

      RecoTau_Fake = 1;
      RecoTau_NNT = 0;
      RecoTau_NCT = 0;
      RecoTau_MaxDPt = 0;
      RecoTau_MaxDAngle = 0;
      RecoTau_Mass = RecoTauFV.M();
      RecoTau_NDaughters = RecoTauDaughters.size();

      for (unsigned int q = 0; q < RecoTauDaughters.size(); q++){
	EVENT::ReconstructedParticle* daughter_track = RecoTauDaughters[q];
	TLorentzVector DFV = TLorentzVector(daughter_track->getMomentum()[0],daughter_track->getMomentum()[1],daughter_track->getMomentum()[2],daughter_track->getEnergy());

	// find the MC particle linked to this Reco particle -  was this daughter from any MC tau?
	// is this tau daughter an isolated lepton (an e/mu?)?
      	for(int r = 0; r < links->getNumberOfElements(); r++){
	  EVENT::LCRelation* linkr = (EVENT::LCRelation*) links->getElementAt(r);
	  EVENT::ReconstructedParticle* rcr = (EVENT::ReconstructedParticle*) linkr->getFrom();
	  if (rcr == daughter_track){
	    EVENT::MCParticle* mcpr = (EVENT::MCParticle*) linkr->getTo();
	    if(std::find(AllMCTauDaughters.begin(), AllMCTauDaughters.end(), mcpr) != AllMCTauDaughters.end()){
	      //cout << "this tau daughter was from an MC tau" << endl;
	      RecoTau_Fake = 0;
	    }
	    if (mcpr == mc_lepton){
	      //cout << "reconstructed tau daughter was an isolated e/mu" << endl;
	      total_taus_MCEMU++;
	    }
	  }
	} 

	if (DFV.Pt() > RecoTau_MaxDPt){
	  RecoTau_MaxDPt = DFV.Pt();
	}
	if (DFV.Angle(RecoTauFV.Vect()) > RecoTau_MaxDAngle){
	  RecoTau_MaxDAngle = DFV.Angle(RecoTauFV.Vect());
	}
	if (daughter_track->getCharge() == 0){
	  RecoTau_NNT++;
	}
	else{
	  RecoTau_NCT++;
	}
      }

      // loop for daughters again, now we know if reco tau is fake or not
      for (unsigned int q = 0; q < RecoTauDaughters.size(); q++){
	EVENT::ReconstructedParticle* daughter_track = RecoTauDaughters[q];
	TLorentzVector DFV = TLorentzVector(daughter_track->getMomentum()[0],daughter_track->getMomentum()[1],daughter_track->getMomentum()[2],daughter_track->getEnergy());

	const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) daughter_track->getTracks();
	if (trkvec.size()>0){
	  RecoTauDaughter_D0 = abs(trkvec[0]->getD0());
	  RecoTauDaughter_Z0 = abs(trkvec[0]->getZ0());
	  RecoTauDaughter_R0 = sqrt( RecoTauDaughter_D0*RecoTauDaughter_D0 + RecoTauDaughter_Z0*RecoTauDaughter_Z0 );
	}
	else{
	  RecoTauDaughter_D0 = 99;
	  RecoTauDaughter_Z0 = 99;
	  RecoTauDaughter_R0 = 99;
	}

	std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) daughter_track->getClusters();
	float ecal = 0;
	float hcal = 0;
	for ( std::vector<EVENT::Cluster*>::const_iterator iCluster=clusters.begin(); iCluster!=clusters.end(); ++iCluster) {
	  ecal += (*iCluster)->getSubdetectorEnergies()[0];
	  hcal += (*iCluster)->getSubdetectorEnergies()[1]; 
	} 
	RecoTauDaughter_CalE = ecal/(ecal+hcal);

	RecoTauDaughter_Pt = DFV.Pt();
	RecoTauDaughter_CosTheta = cos(DFV.Theta());
	RecoTauDaughter_TauAngle = DFV.Angle(RecoTauFV.Vect());
	RecoTau_Fake2 = RecoTau_Fake;
	
	RecoTauDaughter_tree->Fill();
      }

      // counters for output
      if (RecoTau_Fake == 0){
	total_true_taus++;
      }
      else{
	total_fake_taus++;
      }

      // energy within isolation ring about RecoTau 
      RecoTau_RingEnergy = 0;
      RecoTau_NRingParticles = 0;
      RecoTau_MinPFOAngle = 99;
      for (int j = 0; j < pfo->getNumberOfElements(); j++) {
	EVENT::ReconstructedParticle* pfoj = (EVENT::ReconstructedParticle*) pfo->getElementAt(j);
	TLorentzVector pfoj_vec = TLorentzVector(pfoj->getMomentum()[0],pfoj->getMomentum()[1],pfoj->getMomentum()[2],pfoj->getEnergy());

	// double pfoj_D0;
	// double pfoj_Z0;
	// double pfoj_R0;

	// const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) pfoj->getTracks();
	// if (trkvec.size()>0){
	//   pfoj_D0 = abs(trkvec[0]->getD0());
	//   pfoj_Z0 = abs(trkvec[0]->getZ0());
	//   pfoj_R0 = sqrt( pfoj_D0*pfoj_D0 + pfoj_Z0*pfoj_Z0 );
	// }
	// else{
	//   pfoj_D0 = 99;
	//   pfoj_Z0 = 99;
	//   pfoj_R0 = 99;
	// }


	double pfoi_pfoj_angle = pfoj_vec.Angle(RecoTauFV.Vect());
	if (abs(pfoi_pfoj_angle) > 0.04 && abs(pfoi_pfoj_angle) < 0.29){
	  RecoTau_RingEnergy = RecoTau_RingEnergy + pfoj->getEnergy();
	  RecoTau_NRingParticles++;
	}
	// is this the closest pfo to the tau?
	if ((abs(pfoi_pfoj_angle) < RecoTau_MinPFOAngle) && abs(pfoi_pfoj_angle) > 0.04){
	  RecoTau_MinPFOAngle = abs(pfoi_pfoj_angle);
	}
      }
      //cout << "ring energy: " << RecoTau_RingEnergy << " N particles: " << RecoTau_NRingParticles << endl;

      RecoTau_tree->Fill();
    } // end of Reco tau loop

    Event_tree->Fill();
  } // end of event loop

  outfile.Write();

  cout << "of " << total_mc_emus << " possible mc emus" << endl;
  cout << "found " << total_true_emus << " true emus and " << total_fake_emus << " fake emus, of which " << total_emus_MCTD << " are tau daughters" << endl;

  cout << "of " << total_mc_taus << " possible mc taus" << endl;
  cout << "found " << total_true_taus << " true taus and " << total_fake_taus << " fake taus, of which " << total_taus_MCEMU << " are e/mus" << endl;

}

Int_t leptonInvest(std::string recfileName, TString outfileName){
  // to investigate the properties of leptons from W's from top decays
  IO::LCReader* reclcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
  reclcReader->open( recfileName ) ;

  TFile outfile( outfileName,"RECREATE");

  TTree *Event_tree = new TTree("Event_tree","");
  Int_t N_MC_Taus,N_PFO_Selected,N_PFO_Tight,N_PFO_Loose;
  Event_tree->Branch("N_MC_Taus",&N_MC_Taus,"N_MC_Taus/I");
  Event_tree->Branch("N_PFO_Selected",&N_PFO_Selected,"N_PFO_Selected/I");
  Event_tree->Branch("N_PFO_Tight",&N_PFO_Tight,"N_PFO_Tight/I");
  Event_tree->Branch("N_PFO_Loose",&N_PFO_Loose,"N_PFO_Loose/I");
  
  TTree *MCTau_tree = new TTree("MCTau_tree","");
  Double_t MCTau_Mass, MCTau_MaxDPt, MCTau_MaxDAngle, MCTau_CosTheta;
  Int_t MCTau_NDaughters, MCTau_NCT, MCTau_NNT;  
  MCTau_tree->Branch("MCTau_Mass",&MCTau_Mass,"MCTau_Mass/D");
  MCTau_tree->Branch("MCTau_MaxDPt",&MCTau_MaxDPt,"MCTau_MaxDPt/D");
  MCTau_tree->Branch("MCTau_MaxDAngle",&MCTau_MaxDAngle,"MCTau_MaxDAngle/D");
  MCTau_tree->Branch("MCTau_CosTheta",&MCTau_CosTheta,"MCTau_CosTheta/D");
  MCTau_tree->Branch("MCTau_NDaughters",&MCTau_NDaughters,"MCTau_NDaughters/I");
  MCTau_tree->Branch("MCTau_NCT",&MCTau_NCT,"MCTau_NCT/I");
  MCTau_tree->Branch("MCTau_NNT",&MCTau_NNT,"MCTau_NNT/I");

  TTree *RecoTau_tree = new TTree("RecoTau_tree","");
  Double_t RecoTau_Mass, RecoTau_MaxDPt, RecoTau_MaxDAngle, RecoTau_CosTheta, RecoTau_RingEnergy, RecoTau_MinPFOAngle;
  Int_t RecoTau_NDaughters, RecoTau_NCT, RecoTau_NNT, RecoTau_NRingParticles;  
  RecoTau_tree->Branch("RecoTau_Mass",&RecoTau_Mass,"RecoTau_Mass/D");
  RecoTau_tree->Branch("RecoTau_MaxDPt",&RecoTau_MaxDPt,"RecoTau_MaxDPt/D");
  RecoTau_tree->Branch("RecoTau_MaxDAngle",&RecoTau_MaxDAngle,"RecoTau_MaxDAngle/D");
  RecoTau_tree->Branch("RecoTau_CosTheta",&RecoTau_CosTheta,"RecoTau_CosTheta/D");
  RecoTau_tree->Branch("RecoTau_RingEnergy",&RecoTau_RingEnergy,"RecoTau_RingEnergy/D");
  RecoTau_tree->Branch("RecoTau_MinPFOAngle",&RecoTau_MinPFOAngle,"RecoTau_MinPFOAngle/D");
  RecoTau_tree->Branch("RecoTau_NDaughters",&RecoTau_NDaughters,"RecoTau_NDaughters/I");
  RecoTau_tree->Branch("RecoTau_NCT",&RecoTau_NCT,"RecoTau_NCT/I");
  RecoTau_tree->Branch("RecoTau_NNT",&RecoTau_NNT,"RecoTau_NNT/I");
  RecoTau_tree->Branch("RecoTau_NRingParticles",&RecoTau_NRingParticles,"RecoTau_NRingParticles/I");

  TTree *MCTauDaughter_tree = new TTree("MCTauDaughter_tree","");
  Double_t MCTauDaughter_Energy, MCTauDaughter_Pt, MCTauDaughter_CosTheta, MCTauDaughter_TauAngle;
  Int_t MCTauDaughter_PDG;
  MCTauDaughter_tree->Branch("MCTauDaughter_Energy",&MCTauDaughter_Energy,"MCTauDaughter_Energy/D");
  MCTauDaughter_tree->Branch("MCTauDaughter_Pt",&MCTauDaughter_Pt,"MCTauDaughter_Pt/D");
  MCTauDaughter_tree->Branch("MCTauDaughter_CosTheta",&MCTauDaughter_CosTheta,"MCTauDaughter_CosTheta/D");
  MCTauDaughter_tree->Branch("MCTauDaughter_TauAngle",&MCTauDaughter_TauAngle,"MCTauDaughter_TauAngle/D");
  MCTauDaughter_tree->Branch("MCTauDaughter_PDG",&MCTauDaughter_PDG,"MCTauDaughter_PDG/I");

  TTree *RecoTauDaughter_tree = new TTree("RecoTauDaughter_tree","");
  Double_t RecoTauDaughter_Energy, RecoTauDaughter_Pt ,RecoTauDaughter_CosTheta, RecoTauDaughter_TauAngle;
  Double_t RecoTauDaughter_D0, RecoTauDaughter_Z0, RecoTauDaughter_R0, RecoTauDaughter_CalE;
  Int_t RecoTauDaughter_Type;
  RecoTauDaughter_tree->Branch("RecoTauDaughter_Energy",&RecoTauDaughter_Energy,"RecoTauDaughter_Energy/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_Pt",&RecoTauDaughter_Pt,"RecoTauDaughter_Pt/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_CosTheta",&RecoTauDaughter_CosTheta,"RecoTauDaughter_CosTheta/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_TauAngle",&RecoTauDaughter_TauAngle,"RecoTauDaughter_TauAngle/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_D0",&RecoTauDaughter_D0,"RecoTauDaughter_D0/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_Z0",&RecoTauDaughter_Z0,"RecoTauDaughter_Z0/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_R0",&RecoTauDaughter_R0,"RecoTauDaughter_R0/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_CalE",&RecoTauDaughter_CalE,"RecoTauDaughter_CalE/D");
  RecoTauDaughter_tree->Branch("RecoTauDaughter_Type",&RecoTauDaughter_Type,"RecoTauDaughter_Type/I");

  TTree *RecoIsolep_tree = new TTree("RecoIsolep_tree","");
  Int_t RecoIsolep_Type, RecoIsolep_TrueID;
  Double_t RecoIsolep_Energy, RecoIsolep_Pt, RecoIsolep_CosTheta, RecoIsolep_D0, RecoIsolep_Z0, RecoIsolep_R0, RecoIsolep_CalE, RecoIsolep_ConeE;
  RecoIsolep_tree->Branch("RecoIsolep_Type",&RecoIsolep_Type,"RecoIsolep_Type/I");
  RecoIsolep_tree->Branch("RecoIsolep_TrueID",&RecoIsolep_TrueID,"RecoIsolep_TrueID/I");
  RecoIsolep_tree->Branch("RecoIsolep_Energy",&RecoIsolep_Energy,"RecoIsolep_Energy/D");
  RecoIsolep_tree->Branch("RecoIsolep_Pt",&RecoIsolep_Pt,"RecoIsolep_Pt/D");
  RecoIsolep_tree->Branch("RecoIsolep_CosTheta",&RecoIsolep_CosTheta,"RecoIsolep_CosTheta/D");
  RecoIsolep_tree->Branch("RecoIsolep_D0",&RecoIsolep_D0,"RecoIsolep_D0/D");
  RecoIsolep_tree->Branch("RecoIsolep_Z0",&RecoIsolep_Z0,"RecoIsolep_Z0/D");
  RecoIsolep_tree->Branch("RecoIsolep_R0",&RecoIsolep_R0,"RecoIsolep_R0/D");
  RecoIsolep_tree->Branch("RecoIsolep_CalE",&RecoIsolep_CalE,"RecoIsolep_CalE/D");
  RecoIsolep_tree->Branch("RecoIsolep_ConeE",&RecoIsolep_ConeE,"RecoIsolep_ConeE/D");

  TTree *RecoPFO_tree = new TTree("RecoPFO_tree","");
  Double_t RecoPFO_Energy, RecoPFO_Pt, RecoPFO_CosTheta, RecoPFO_RingEnergy, RecoPFO_MinPFOAngle;
  Double_t RecoPFO_D0, RecoPFO_Z0, RecoPFO_R0, RecoPFO_CalE, RecoPFO_ConeE;
  Int_t RecoPFO_Type, RecoPFO_NRingParticles;
  RecoPFO_tree->Branch("RecoPFO_Energy",&RecoPFO_Energy,"RecoPFO_Energy/D");
  RecoPFO_tree->Branch("RecoPFO_Pt",&RecoPFO_Pt,"RecoPFO_Pt/D");
  RecoPFO_tree->Branch("RecoPFO_CosTheta",&RecoPFO_CosTheta,"RecoPFO_CosTheta/D");
  RecoPFO_tree->Branch("RecoPFO_RingEnergy",&RecoPFO_RingEnergy,"RecoPFO_RingEnergy/D");
  RecoPFO_tree->Branch("RecoPFO_D0",&RecoPFO_D0,"RecoPFO_D0/D");
  RecoPFO_tree->Branch("RecoPFO_Z0",&RecoPFO_Z0,"RecoPFO_Z0/D");
  RecoPFO_tree->Branch("RecoPFO_R0",&RecoPFO_R0,"RecoPFO_R0/D");
  RecoPFO_tree->Branch("RecoPFO_CalE",&RecoPFO_CalE,"RecoPFO_CalE/D");
  RecoPFO_tree->Branch("RecoPFO_MinPFOAngle",&RecoPFO_MinPFOAngle,"RecoPFO_MinPFOAngle/D");
  RecoPFO_tree->Branch("RecoPFO_ConeE",&RecoPFO_ConeE,"RecoPFO_ConeE/D");
  RecoPFO_tree->Branch("RecoPFO_Type",&RecoPFO_Type,"RecoPFO_Type/I");
  RecoPFO_tree->Branch("RecoPFO_NRingParticles",&RecoPFO_NRingParticles,"RecoPFO_NRingParticles/I");

  TTree *RecoPFOTight_tree = new TTree("RecoPFOTight_tree","");
  Double_t RecoPFOTight_Energy;
  RecoPFOTight_tree->Branch("RecoPFOTight_Energy",&RecoPFOTight_Energy,"RecoPFOTight_Energy/D");

  TTree *RecoPFOLoose_tree = new TTree("RecoPFOLoose_tree","");
  Double_t RecoPFOLoose_Energy;
  RecoPFOLoose_tree->Branch("RecoPFOLoose_Energy",&RecoPFOLoose_Energy,"RecoPFOLoose_Energy/D");

  // do a seed tree?
  // seed pt

  EVENT::LCEvent* recevent = 0;
  Int_t nEvt = 0;

  while( (recevent = reclcReader->readNextEvent()) != 0){
    nEvt++;
    //cout << "================" << endl;

    EVENT::LCCollection* links = recevent->getCollection("RecoMCTruthLink");
    EVENT::LCCollection* pfo = recevent->getCollection("SelectedPandoraPFOCollection");
    EVENT::LCCollection* pfo_loose = recevent->getCollection("LooseSelectedPandoraPFOCollection");
    EVENT::LCCollection* pfo_tight = recevent->getCollection("TightSelectedPandoraPFOCollection");

    N_PFO_Selected = pfo->getNumberOfElements();
    N_PFO_Tight = pfo_tight->getNumberOfElements();
    N_PFO_Loose = pfo_loose->getNumberOfElements();

    std::vector<EVENT::ReconstructedParticle*> RecoTauDaughters;
    RecoTauDaughters.clear();

    //per tau
    std::vector<EVENT::MCParticle*> tdvec;
    std::vector<EVENT::ReconstructedParticle*> rtdvec;
    // per event
    std::vector<EVENT::MCParticle*> evt_tdvec;
    std::vector<EVENT::ReconstructedParticle*> evt_rtdvec;


    EVENT::ReconstructedParticle* recod_lepton = getReconstructedLepton(recevent); // 11 or 13 only
    if (recod_lepton!=0){
      //cout << "found generated lepton in reconstructed particles" << endl;
      TLorentzVector LFV = TLorentzVector(recod_lepton->getMomentum()[0],recod_lepton->getMomentum()[1],recod_lepton->getMomentum()[2],recod_lepton->getEnergy());
      RecoIsolep_Type = recod_lepton->getType();
      EVENT::MCParticle* mc_lepton = getMCLepton(recevent); // used to be 11 or 13 only, not anymore
      RecoIsolep_TrueID = mc_lepton->getPDG();
      RecoIsolep_Energy = recod_lepton->getEnergy();
      RecoIsolep_Pt = LFV.Pt();
      RecoIsolep_CosTheta = cos(LFV.Theta());

      const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) recod_lepton->getTracks();
      if (trkvec.size()>0){
	RecoIsolep_D0 = abs(trkvec[0]->getD0());
	RecoIsolep_Z0 = abs(trkvec[0]->getZ0());
	RecoIsolep_R0 = sqrt( RecoIsolep_D0*RecoIsolep_D0 + RecoIsolep_Z0*RecoIsolep_Z0 );
      }
      else{
	RecoIsolep_D0 = 99;
	RecoIsolep_Z0 = 99;
	RecoIsolep_R0 = 99;
      }

      std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) recod_lepton->getClusters();
      float ecal = 0;
      float hcal = 0;
      for ( std::vector<EVENT::Cluster*>::const_iterator iCluster=clusters.begin(); iCluster!=clusters.end(); ++iCluster) {
	ecal += (*iCluster)->getSubdetectorEnergies()[0];
	hcal += (*iCluster)->getSubdetectorEnergies()[1];
      }
      RecoIsolep_CalE = ecal/(ecal+hcal);

      RecoIsolep_ConeE = 0;
      for (int j = 0; j < pfo->getNumberOfElements(); j++) {
	EVENT::ReconstructedParticle* pfoj = (EVENT::ReconstructedParticle*) pfo->getElementAt(j);
	if (pfoj == recod_lepton){
	  continue;
	}
	TLorentzVector pfoj_vec = TLorentzVector(pfoj->getMomentum()[0],pfoj->getMomentum()[1],pfoj->getMomentum()[2],pfoj->getEnergy());

	double LFV_pfoj_angle = pfoj_vec.Angle(LFV.Vect());
	if (cos(LFV_pfoj_angle) > 0.995){
	  RecoIsolep_ConeE = RecoIsolep_ConeE + pfoj->getEnergy();
	}
      }

      RecoIsolep_tree->Fill();
    }

    // get the MC taus from top decays which have daughters
    std::vector<EVENT::MCParticle*> MCTaus;
    MCTaus.clear();
    MCTaus = getMCTaus(recevent);
    //cout << "found " << MCTaus.size() << " MC taus from top decay which have daughters" << endl;
    N_MC_Taus = MCTaus.size();

    // for each MC tau, get the daughters
    for (unsigned int p = 0; p < MCTaus.size(); p++){
      EVENT::MCParticle* MCTau = MCTaus[p];
      TLorentzVector MFV = TLorentzVector(MCTau->getMomentum()[0],MCTau->getMomentum()[1],MCTau->getMomentum()[2],MCTau->getEnergy());
      TLorentzVector RecoTauFV = TLorentzVector(0,0,0,0);
      std::vector<EVENT::MCParticle*> MCTauDaughters;
      MCTauDaughters.clear();
      MCTauDaughters = getMCTauDaughters(MCTau);
      //cout << "found " << MCTauDaughters.size() << " daughters for MCTau " << p << endl;

      MCTau_Mass = MCTau->getMass();
      MCTau_NDaughters = MCTauDaughters.size();
      MCTau_MaxDPt = 0;
      MCTau_MaxDAngle = 0;
      MCTau_NCT = 0;
      MCTau_NNT = 0;
      MCTau_CosTheta = cos(MFV.Theta());

      RecoTau_NDaughters = 0;
      RecoTau_MaxDPt = 0;
      RecoTau_MaxDAngle = 0;
      RecoTau_NCT = 0;
      RecoTau_NNT = 0;

      // for each daughter, fill properties
      for (unsigned int q = 0; q < MCTauDaughters.size(); q++){
	EVENT::MCParticle* MCTauDaughter = MCTauDaughters[q];
	TLorentzVector FV = TLorentzVector(MCTauDaughter->getMomentum()[0],MCTauDaughter->getMomentum()[1],MCTauDaughter->getMomentum()[2],MCTauDaughter->getEnergy());
	MCTauDaughter_Energy = MCTauDaughter->getEnergy();
	MCTauDaughter_Pt = FV.Pt();
	MCTauDaughter_CosTheta = cos(FV.Theta());
	MCTauDaughter_PDG = MCTauDaughter->getPDG();
	MCTauDaughter_TauAngle = FV.Angle(MFV.Vect());

	if (MCTauDaughter_Pt > MCTau_MaxDPt){
	  MCTau_MaxDPt = MCTauDaughter_Pt;
	}
	if (MCTauDaughter_TauAngle > MCTau_MaxDAngle){
	  MCTau_MaxDAngle = MCTauDaughter_TauAngle;
	}

	if (MCTauDaughter->getCharge() == 0){
	  MCTau_NNT++;
	}
	else{
	  MCTau_NCT++;
	}

	MCTauDaughter_tree->Fill();

	// is there a linked reconstructed particle?
	for (int i = 0; i < links->getNumberOfElements(); i++){
	  EVENT::LCRelation* linki = (EVENT::LCRelation*) links->getElementAt(i);
	  EVENT::MCParticle* mcpi = (EVENT::MCParticle*) linki->getTo();
	  if(mcpi == MCTauDaughter) {
	    EVENT::ReconstructedParticle* rpi = (EVENT::ReconstructedParticle*) linki->getFrom();
	    RecoTauDaughters.push_back(rpi);
	    TLorentzVector RFV = TLorentzVector(rpi->getMomentum()[0],rpi->getMomentum()[1],rpi->getMomentum()[2],rpi->getEnergy());
	    RecoTauFV = RecoTauFV + RFV;

	    RecoTauDaughter_Energy = rpi->getEnergy();
	    RecoTauDaughter_Pt = RFV.Pt();
	    RecoTauDaughter_CosTheta = cos(RFV.Theta());
	    RecoTauDaughter_Type = rpi->getType();
	    RecoTauDaughter_TauAngle = FV.Angle(RFV.Vect());
	    RecoTau_NDaughters++;

	    std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) rpi->getClusters();
	    float ecal = 0;
	    float hcal = 0;
	    for ( std::vector<EVENT::Cluster*>::const_iterator iCluster=clusters.begin(); iCluster!=clusters.end(); ++iCluster) {
	      ecal += (*iCluster)->getSubdetectorEnergies()[0];
	      hcal += (*iCluster)->getSubdetectorEnergies()[1];
	    }
	    RecoTauDaughter_CalE = ecal/(ecal+hcal);

	    const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) rpi->getTracks();
	    if (trkvec.size()>0){
	      RecoTauDaughter_D0 = abs(trkvec[0]->getD0());
	      RecoTauDaughter_Z0 = abs(trkvec[0]->getZ0());
	      RecoTauDaughter_R0 = sqrt( RecoTauDaughter_D0*RecoTauDaughter_D0 + RecoTauDaughter_Z0*RecoTauDaughter_Z0 );
	    }
	    else{
	      RecoTauDaughter_D0 = 99;
	      RecoTauDaughter_Z0 = 99;
	      RecoTauDaughter_R0 = 99;
	    }

	    if (RecoTauDaughter_Pt > RecoTau_MaxDPt){
	      RecoTau_MaxDPt = RecoTauDaughter_Pt;
	    }
	    if (RecoTauDaughter_TauAngle > RecoTau_MaxDAngle){
	      RecoTau_MaxDAngle = RecoTauDaughter_TauAngle;
	    }

	    if (rpi->getCharge() == 0){
	      RecoTau_NNT++;
	    }
	    else{
	      RecoTau_NCT++;
	    }

	    RecoTauDaughter_tree->Fill();
	  }
	}
      }

      //if there was a RecoTau
      if (RecoTau_NDaughters > 0){
	RecoTau_Mass = RecoTauFV.M();
	RecoTau_CosTheta = cos(RecoTauFV.Theta());

	// energy within isolation ring about RecoTau 
	RecoTau_RingEnergy = 0;
	RecoTau_NRingParticles = 0;
	RecoTau_MinPFOAngle = 99;
	for (int j = 0; j < pfo->getNumberOfElements(); j++) {
	  EVENT::ReconstructedParticle* pfoj = (EVENT::ReconstructedParticle*) pfo->getElementAt(j);
	  TLorentzVector pfoj_vec = TLorentzVector(pfoj->getMomentum()[0],pfoj->getMomentum()[1],pfoj->getMomentum()[2],pfoj->getEnergy());
	  
	  double pfoj_D0;
	  double pfoj_Z0;
	  double pfoj_R0;
	  
	  const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) pfoj->getTracks();
	  if (trkvec.size()>0){
	    pfoj_D0 = abs(trkvec[0]->getD0());
	    pfoj_Z0 = abs(trkvec[0]->getZ0());
	    pfoj_R0 = sqrt( pfoj_D0*pfoj_D0 + pfoj_Z0*pfoj_Z0 );
	  }
	  else{
	    pfoj_D0 = 99;
	    pfoj_Z0 = 99;
	    pfoj_R0 = 99;
	  }
	  
	  // if (pfoj_vec.Pt() < 5 || pfoj_R0 > 0.02){
	  //   continue;
	  // }
	  double pfoi_pfoj_angle = pfoj_vec.Angle(RecoTauFV.Vect());
	  if (abs(pfoi_pfoj_angle) > 0.04 && abs(pfoi_pfoj_angle) < 0.29){
	    RecoTau_RingEnergy = RecoTau_RingEnergy + pfoj->getEnergy();
	    RecoTau_NRingParticles++;
	  }
	  // is this the closest pfo to the tau?
	  if ((abs(pfoi_pfoj_angle) < RecoTau_MinPFOAngle) && abs(pfoi_pfoj_angle) > 0.04){
	    RecoTau_MinPFOAngle = abs(pfoi_pfoj_angle);
	  }
	}
	RecoTau_tree->Fill();
      }

      MCTau_tree->Fill();
    }

    // properties of PFOs, not including RecoTauDaughters
    for (int i = 0; i < pfo->getNumberOfElements(); i++) {
      EVENT::ReconstructedParticle* pfoi = (EVENT::ReconstructedParticle*) pfo->getElementAt(i);
      // if the pfo is any tau daughter, continue
      if(std::find(RecoTauDaughters.begin(), RecoTauDaughters.end(), pfoi) != RecoTauDaughters.end()) {
	continue;
      }
      TLorentzVector PFV = TLorentzVector(pfoi->getMomentum()[0],pfoi->getMomentum()[1],pfoi->getMomentum()[2],pfoi->getEnergy());

      // energy within isolation ring about pfoi 
      RecoPFO_RingEnergy = 0;
      RecoPFO_ConeE = 0;
      RecoPFO_NRingParticles = 0;
      RecoPFO_MinPFOAngle = 99;
      for (int j = 0; j < pfo->getNumberOfElements(); j++) {
	EVENT::ReconstructedParticle* pfoj = (EVENT::ReconstructedParticle*) pfo->getElementAt(j);
	TLorentzVector pfoj_vec = TLorentzVector(pfoj->getMomentum()[0],pfoj->getMomentum()[1],pfoj->getMomentum()[2],pfoj->getEnergy());

	double pfoj_D0;
	double pfoj_Z0;
	double pfoj_R0;

	const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) pfoj->getTracks();
	if (trkvec.size()>0){
	  pfoj_D0 = abs(trkvec[0]->getD0());
	  pfoj_Z0 = abs(trkvec[0]->getZ0());
	  pfoj_R0 = sqrt( pfoj_D0*pfoj_D0 + pfoj_Z0*pfoj_Z0 );
	}
	else{
	  pfoj_D0 = 99;
	  pfoj_Z0 = 99;
	  pfoj_R0 = 99;
	}

	// if (pfoj_vec.Pt() < 5 || pfoj_R0 < 0.02){
	//   continue;
	// }
	double pfoi_pfoj_angle = pfoj_vec.Angle(PFV.Vect());
	if (abs(pfoi_pfoj_angle) > 0.04 && abs(pfoi_pfoj_angle) < 0.29){
	  RecoPFO_RingEnergy = RecoPFO_RingEnergy + pfoj->getEnergy();
	  RecoPFO_NRingParticles++;
	}

	// new
	if (cos(pfoi_pfoj_angle) > 0.995){
	  if (pfoj != pfoi){
	    RecoPFO_ConeE = RecoPFO_ConeE + pfoj->getEnergy();
	  }
	}


	// is this the closest pfo to the pfo?
	if ((abs(pfoi_pfoj_angle) < RecoPFO_MinPFOAngle) && abs(pfoi_pfoj_angle) > 0.04){
	  RecoPFO_MinPFOAngle = abs(pfoi_pfoj_angle);
	}
      }

      RecoPFO_Energy = pfoi->getEnergy();
      RecoPFO_Pt = PFV.Pt();
      RecoPFO_CosTheta = cos(PFV.Theta());
      RecoPFO_Type = pfoi->getType();

      std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) pfoi->getClusters();
      float ecal = 0;
      float hcal = 0;
      for ( std::vector<EVENT::Cluster*>::const_iterator iCluster=clusters.begin(); iCluster!=clusters.end(); ++iCluster) {
	ecal += (*iCluster)->getSubdetectorEnergies()[0];
	hcal += (*iCluster)->getSubdetectorEnergies()[1];
      }
      RecoPFO_CalE = ecal/(ecal+hcal);

      const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) pfoi->getTracks();
      if (trkvec.size()>0){
	RecoPFO_D0 = abs(trkvec[0]->getD0());
	RecoPFO_Z0 = abs(trkvec[0]->getZ0());
	RecoPFO_R0 = sqrt( RecoPFO_D0*RecoPFO_D0 + RecoPFO_Z0*RecoPFO_Z0 );
      }
      else{
	RecoPFO_D0 = 99;
	RecoPFO_Z0 = 99;
	RecoPFO_R0 = 99;
      }

      RecoPFO_tree->Fill();
    }

    for (int i = 0; i < pfo_loose->getNumberOfElements(); i++) {
      EVENT::ReconstructedParticle* pfoi = (EVENT::ReconstructedParticle*) pfo_loose->getElementAt(i);
      if(std::find(RecoTauDaughters.begin(), RecoTauDaughters.end(), pfoi) != RecoTauDaughters.end()) {
        continue;
      }
      RecoPFOLoose_Energy = pfoi->getEnergy();
      RecoPFOLoose_tree->Fill();
    }

    for (int i = 0; i < pfo_tight->getNumberOfElements(); i++) {
      EVENT::ReconstructedParticle* pfoi = (EVENT::ReconstructedParticle*) pfo_tight->getElementAt(i);
      if(std::find(RecoTauDaughters.begin(), RecoTauDaughters.end(), pfoi) != RecoTauDaughters.end()) {
        continue;
      }
      RecoPFOTight_Energy = pfoi->getEnergy();
      RecoPFOTight_tree->Fill();
    }

    Event_tree->Fill();

  } // end of event loop
  outfile.Write();
  cout << "Processed " << recfileName << endl;
  return nEvt;

}

// Functions copied from MCMatching
// Edit there, not here

EVENT::MCParticle* getGeneratedLepton(EVENT::LCEvent* event){
  // to return the generated signal lepton status 3
  EVENT::MCParticle* pointer_to_gen_lepton = 0;

  EVENT::LCCollection* mcp = event->getCollection("MCParticle");
  for (int i = 0; i < mcp->getNumberOfElements(); i++){
    EVENT::MCParticle* mcpi = (EVENT::MCParticle*) mcp->getElementAt(i);
    if (mcpi->getGeneratorStatus() == 3){
      if (abs(mcpi->getPDG()) == 11 || abs(mcpi->getPDG()) == 13 || abs(mcpi->getPDG()) == 15) {
        if(mcpi->getParents().size() > 0){
          if (abs(mcpi->getParents()[0]->getPDG()) == 24){
            pointer_to_gen_lepton = mcpi;
          }
        }
      }
    }
  }
  //cout << "got generated lepton " << pointer_to_gen_lepton->getPDG() << endl;
  
  return pointer_to_gen_lepton;
}

EVENT::MCParticle* getMCLepton(EVENT::LCEvent* event){
  // to return the mc signal lepton
  EVENT::MCParticle* pointer_to_mc_lepton = 0;

  EVENT::LCCollection* mcp = event->getCollection("MCParticle");
  for (int i = 0; i < mcp->getNumberOfElements(); i++){
    EVENT::MCParticle* mcpi = (EVENT::MCParticle*) mcp->getElementAt(i);
    if ((abs(mcpi->getPDG()) == 11) || (abs(mcpi->getPDG()) == 13)){
      if (mcpi->getGeneratorStatus() == 1){
        if (mcpi->getParents().size() == 1){
          if (abs(mcpi->getParents()[0]->getPDG()) == 24){
            pointer_to_mc_lepton = mcpi;
          }
        }
      }
    }
    if (abs(mcpi->getPDG()) == 15){
      if (mcpi->getParents().size() > 0){
        if (mcpi->getParents()[0]->getParents().size() > 0){
          if (abs(mcpi->getParents()[0]->getPDG()) == 24 && (abs(mcpi->getParents()[0]->getParents()[0]->getPDG()) == 94 || abs(mcpi->getParents()[0]->getParents()[0]->getPDG()) == 15 || abs(mcpi->getParents()[0]->getParents()[0]->getPDG()) == 16)){
            if (mcpi->getDaughters().size() > 0){
              pointer_to_mc_lepton = mcpi;
            }
          }
        }
      }
    }
  }
  /*
  if (pointer_to_mc_lepton){
    cout << "got mc lepton " << pointer_to_mc_lepton->getPDG() << endl;
  }
  else{
    cout << "no MC lepton" << endl;
  }
  */
  return pointer_to_mc_lepton;
}

EVENT::ReconstructedParticle* getReconstructedLepton(EVENT::LCEvent* event){
  // to return the reconstructed e/mu signal lepton
  // if a tau was generated, return nothing

  EVENT::LCRelation* link_to_rec_lepton = 0;
  EVENT::MCParticle* pointer_to_mc_lepton = 0;
  EVENT::ReconstructedParticle* pointer_to_rec_lepton = 0;
  EVENT::LCCollection* links = event->getCollection("RecoMCTruthLink");

  pointer_to_mc_lepton = getMCLepton(event);

  for (int i = 0; i < links->getNumberOfElements(); i++){
    EVENT::LCRelation* linki = (EVENT::LCRelation*) links->getElementAt(i);
    EVENT::MCParticle* mcpi = (EVENT::MCParticle*) linki->getTo();
    if (pointer_to_mc_lepton == mcpi){
      link_to_rec_lepton = linki;
    }
  }

  // follow the link
  // if no link, 0 returned
  if (link_to_rec_lepton){
    EVENT::ReconstructedParticle* rpi = (EVENT::ReconstructedParticle*) link_to_rec_lepton->getFrom();
    pointer_to_rec_lepton = rpi;
    //cout << "got linked rec lepton " << rpi->getType() << endl;
  }

  return pointer_to_rec_lepton;

}
