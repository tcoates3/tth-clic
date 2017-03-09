#define tth_analysis_cxx
#include "tth_analysis.h"
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TObjArray.h>

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

#include "jama_eig.h"
#include "tnt_math_utils.h"

using namespace std;

void eventShapes(std::string fileName, TString outfileName, std::string LeptonCollection, std::string JetCollection, std::string PFOCollection, std::string PFOsInJetsCollection) {
  // to investigate event shape variables
  IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
  lcReader->open( fileName ) ;

  TFile outfile( outfileName,"RECREATE");

  TTree *tth_tree = new TTree("tth_tree","");

  Double_t pfo_thrust, pfo_oblateness, pfo_sphericity, pfo_aplanarity;
  Double_t pfos_injets_thrust, pfos_injets_oblateness, pfos_injets_sphericity, pfos_injets_aplanarity;
  Double_t jet_thrust, jet_oblateness, jet_sphericity, jet_aplanarity;
  Double_t mcjet_thrust, mcjet_oblateness, mcjet_sphericity, mcjet_aplanarity;


  tth_tree->Branch("pfo_thrust",&pfo_thrust,"pfo_thrust/D");
  tth_tree->Branch("pfo_oblateness",&pfo_oblateness,"pfo_oblateness/D");
  tth_tree->Branch("pfo_sphericity",&pfo_sphericity,"pfo_sphericity/D");
  tth_tree->Branch("pfo_aplanarity",&pfo_aplanarity,"pfo_aplanarity/D");

  tth_tree->Branch("pfos_injets_thrust",&pfos_injets_thrust,"pfos_injets_thrust/D");
  tth_tree->Branch("pfos_injets_oblateness",&pfos_injets_oblateness,"pfos_injets_oblateness/D");
  tth_tree->Branch("pfos_injets_sphericity",&pfos_injets_sphericity,"pfos_injets_sphericity/D");
  tth_tree->Branch("pfos_injets_aplanarity",&pfos_injets_aplanarity,"pfos_injets_aplanarity/D");

  tth_tree->Branch("jet_thrust",&jet_thrust,"jet_thrust/D");
  tth_tree->Branch("jet_oblateness",&jet_oblateness,"jet_oblateness/D");
  tth_tree->Branch("jet_sphericity",&jet_sphericity,"jet_sphericity/D");
  tth_tree->Branch("jet_aplanarity",&jet_aplanarity,"jet_aplanarity/D");

  tth_tree->Branch("mcjet_thrust",&mcjet_thrust,"mcjet_thrust/D");
  tth_tree->Branch("mcjet_oblateness",&mcjet_oblateness,"mcjet_oblateness/D");
  tth_tree->Branch("mcjet_sphericity",&mcjet_sphericity,"mcjet_sphericity/D");
  tth_tree->Branch("mcjet_aplanarity",&mcjet_aplanarity,"mcjet_aplanarity/D");



  EVENT::LCEvent* event = 0;

  while( (event = lcReader->readNextEvent()) != 0 ) {
    //for(int k = 0; k < 50; k++){
    //event = lcReader->readNextEvent();

    // was one lepton found?
    EVENT::LCCollection* leptons = event->getCollection(LeptonCollection);
    Int_t nLeptons;
    nLeptons = leptons->getNumberOfElements();
    if (nLeptons != 1){
      //cout << "NLeptons not 1, continuing." << endl;
      continue;
   }

    // if the wrong lepton was found, continue
    EVENT::ReconstructedParticle* recod_lepton = getReconstructedLepton(event);
    if (recod_lepton == 0){
      //cout << "No link returned, continuing." << endl;
      continue;
    }

    EVENT::ReconstructedParticle* leptoni = (EVENT::ReconstructedParticle*) leptons->getElementAt(0);
    if (leptoni != recod_lepton){
      //cout << "Isolated lepton not signal lepton, continuing." << endl;
      continue;
    }

    TObjArray *veclist_mcjets = new TObjArray();
    veclist_mcjets->SetOwner(kTRUE);

    // find the mc jets
    EVENT::LCCollection* rec_mcp = event->getCollection("MCParticle");
    for (int i = 0; i < rec_mcp->getNumberOfElements(); i++){
      EVENT::MCParticle* rec_mcpi = (EVENT::MCParticle*) rec_mcp->getElementAt(i);
      if (rec_mcpi->getGeneratorStatus() == 3){
	if (abs(rec_mcpi->getPDG()) == 24) {
	  if (abs(rec_mcpi->getDaughters()[0]->getPDG()) < 11){
	    TVector3 vec3_tight1(rec_mcpi->getDaughters()[0]->getMomentum());
	    veclist_mcjets->Add((TVector3*)vec3_tight1.Clone());
	    TVector3 vec3_tight2(rec_mcpi->getDaughters()[1]->getMomentum());
	    veclist_mcjets->Add((TVector3*)vec3_tight2.Clone());
	  }
	}
	if (abs(rec_mcpi->getPDG()) == 6) {
	  //cout << "got a top" << endl;
	  for (int i = 0; i < rec_mcpi->getDaughters().size(); i++){
	    if (rec_mcpi->getDaughters()[i]->getPDG() == 5) {
	      TVector3 vec3_tight(rec_mcpi->getDaughters()[i]->getMomentum());
	      veclist_mcjets->Add((TVector3*)vec3_tight.Clone());
	    }
	    if (rec_mcpi->getDaughters()[i]->getPDG() == -5) {
	      TVector3 vec3_tight(rec_mcpi->getDaughters()[i]->getMomentum());
	      veclist_mcjets->Add((TVector3*)vec3_tight.Clone());
	    }
	  }
	}
      }
    }
    // find the higgs daughters
    for (int i = 0; i < rec_mcp->getNumberOfElements(); i++){
      EVENT::MCParticle* rec_mcpi = (EVENT::MCParticle*) rec_mcp->getElementAt(i);
      if (rec_mcpi->getGeneratorStatus() == 2){
	if (rec_mcpi->getParents().size() == 1){
	  if (abs(rec_mcpi->getParents()[0]->getPDG()) == 25){
	    if (rec_mcpi->getPDG() == 5) {
	      TVector3 vec3_tight(rec_mcpi->getMomentum());
	      veclist_mcjets->Add((TVector3*)vec3_tight.Clone());
	    }
	    if (rec_mcpi->getPDG() == -5) {
	      TVector3 vec3_tight(rec_mcpi->getMomentum());
	      veclist_mcjets->Add((TVector3*)vec3_tight.Clone());
	    }
	  }
	}
      }
    }

    TObjArray *veclist_pfos = new TObjArray();
    veclist_pfos->SetOwner(kTRUE);
    TObjArray *veclist_pfos_injets = new TObjArray();
    veclist_pfos_injets->SetOwner(kTRUE);
    TObjArray *veclist_jets = new TObjArray();
    veclist_jets->SetOwner(kTRUE);

    EVENT::LCCollection* pfos = event->getCollection(PFOCollection);
    for (int i = 0; i < pfos->getNumberOfElements(); i++) {
      EVENT::ReconstructedParticle* pfoi = (EVENT::ReconstructedParticle*) pfos->getElementAt(i);
      TVector3 vec3_tight(pfoi->getMomentum());
      veclist_pfos->Add((TVector3*)vec3_tight.Clone());
    }

    EVENT::LCCollection* pfos_injets = event->getCollection(PFOsInJetsCollection);
    for (int i = 0; i < pfos_injets->getNumberOfElements(); i++) {
      EVENT::ReconstructedParticle* pfoi = (EVENT::ReconstructedParticle*) pfos_injets->getElementAt(i);
      TVector3 vec3_tight(pfoi->getMomentum());
      veclist_pfos_injets->Add((TVector3*)vec3_tight.Clone());
    }

    EVENT::LCCollection* jets = event->getCollection(JetCollection);
    for (int i = 0; i < jets->getNumberOfElements(); i++) {
      EVENT::ReconstructedParticle* jeti = (EVENT::ReconstructedParticle*) jets->getElementAt(i);
      TVector3 vec3_tight(jeti->getMomentum());
      veclist_jets->Add((TVector3*)vec3_tight.Clone());
    }

    CalculateEventShapeVariables(veclist_pfos, pfo_thrust, pfo_oblateness, pfo_sphericity, pfo_aplanarity);
    //cout << pfo_thrust << pfo_oblateness << pfo_sphericity << pfo_aplanarity << endl;

    CalculateEventShapeVariables(veclist_pfos_injets, pfos_injets_thrust, pfos_injets_oblateness, pfos_injets_sphericity, pfos_injets_aplanarity);
    //cout << pfos_injets_thrust << pfos_injets_oblateness << pfos_injets_sphericity << pfos_injets_aplanarity << endl;

    CalculateEventShapeVariables(veclist_jets, jet_thrust, jet_oblateness, jet_sphericity, jet_aplanarity);
    //cout << jet_thrust << jet_oblateness << jet_sphericity << jet_aplanarity << endl;
    CalculateEventShapeVariables(veclist_mcjets, mcjet_thrust, mcjet_oblateness, mcjet_sphericity, mcjet_aplanarity);
    //cout << mcjet_thrust << mcjet_oblateness << mcjet_sphericity << mcjet_aplanarity << endl;
    
    tth_tree->Fill();

  } // end of event loop

  outfile.Write();
  cout << "Made " << outfileName << endl;

}

void CalculateEventShapeVariables(TObjArray* e, Double_t &event_thrust, Double_t &event_oblateness, Double_t &event_sphericity , Double_t &event_aplanarity ) {	

  const Int_t m_maxpart = 1000;

  // Parameters
  // ==========

  // Double_t m_dSphMomPower = 2.0;
  // PARU(41): Power of momentum dependence in sphericity finder.

  Double_t m_dDeltaThPower = 0;
  // PARU(42): Power of momentum dependence in thrust finder.   

  Int_t m_iFast = 4;
  // MSTU(44): # of initial fastest particles choosen to start search.

  Double_t m_dConv = 0.0001;
  // PARU(48): Convergence criteria for axis maximization.

  Int_t m_iGood = 2;
  // MSTU(45): # different starting configurations that must
  // converge before axis is accepted as correct.

  // Variables
  // =========

  TMatrixD m_dAxes;
  m_dAxes.ResizeTo(4,4);

  TRandom m_random;

  Double_t m_dThrust[4];
  Double_t m_dOblateness;

  //To make this look like normal physics notation the
  //zeroth element of each array, mom[i][0], will be ignored
  //and operations will be on elements 1,2,3...
  TMatrixD mom(m_maxpart,6);
  Double_t tmax = 0;
  Double_t phi = 0.;
  Double_t the = 0.;
  Double_t sgn;
  TMatrixD fast(m_iFast + 1,6);
  TMatrixD work(11,6);
  Double_t tdi[] = {0.,0.,0.,0.};
  Double_t tds;
  Double_t tpr[] = {0.,0.,0.,0.};
  Double_t thp;
  Double_t thps;
  TMatrixD temp(3,5);
  Int_t np = 0;
  Int_t numElements = e->GetEntries();
  TObject* o;
  for(Int_t elem=0;elem<numElements;elem++) {
    o = e->At(elem);
    
    if (np >= m_maxpart) { 
      printf("Too many particles input to LCDEventShape");
      return;
    }

    TString nam(o->IsA()->GetName());
    if (nam.Contains("TVector3")) {
      mom(np,1) = ((TVector3 *) o)->X();
      mom(np,2) = ((TVector3 *) o)->Y();
      mom(np,3) = ((TVector3 *) o)->Z();
      mom(np,4) = TMath::Sqrt(mom(np,1)*mom(np,1) + mom(np,2)*mom(np,2)
			      + mom(np,3)*mom(np,3));
    } else if (nam.Contains("TLorentzVector")) {
      mom(np,1) = ((TLorentzVector *) o)->X();
      mom(np,2) = ((TLorentzVector *) o)->Y();
      mom(np,3) = ((TLorentzVector *) o)->Z();
      mom(np,4) = TMath::Sqrt(mom(np,1)*mom(np,1) + mom(np,2)*mom(np,2)
			      + mom(np,3)*mom(np,3));
    } else {
      printf("LCDEventShape::SetEvent input is not a TVector3 or a TLorentzVector\n");
      continue;
    }

    if ( TMath::Abs( m_dDeltaThPower ) <= 0.001 ) {
      mom(np,5) = 1.0;
    } else {
      mom(np,5) = TMath::Power(mom(np,4),m_dDeltaThPower);
    }
    tmax = tmax + mom(np,4)*mom(np,5);
    np++;
  }
  if ( np < 2 ) {
    m_dThrust[1] = -1.0;
    m_dOblateness = -1.0;
    return;
  }
  // for pass = 1: find thrust axis.
  // for pass = 2: find major axis.
  for ( Int_t pass=1; pass < 3; pass++) {
    if ( pass == 2 ) {
      phi = ulAngle(m_dAxes(1,1), m_dAxes(1,2));
      ludbrb( &mom, 0, -phi, 0., 0., 0. );
      for ( Int_t i = 0; i < 3; i++ ) {
	for ( Int_t j = 1; j < 4; j++ ) {
	  temp(i,j) = m_dAxes(i+1,j);
	}
	temp(i,4) = 0;
      }
      ludbrb(&temp,0.,-phi,0.,0.,0.);
      for ( Int_t ib = 0; ib < 3; ib++ ) {
	for ( Int_t j = 1; j < 4; j++ ) {
	  m_dAxes(ib+1,j) = temp(ib,j);
	}
      }
      the = ulAngle( m_dAxes(1,3), m_dAxes(1,1) );
      ludbrb( &mom, -the, 0., 0., 0., 0. );
      for ( Int_t ic = 0; ic < 3; ic++ ) {
	for ( Int_t j = 1; j < 4; j++ ) {
	  temp(ic,j) = m_dAxes(ic+1,j);
	}
	temp(ic,4) = 0;
      }
      ludbrb(&temp,-the,0.,0.,0.,0.);
      for ( Int_t id = 0; id < 3; id++ ) {	
	for ( Int_t j = 1; j < 4; j++ ) {
	  m_dAxes(id+1,j) = temp(id,j);
	}
      }
    }
    for ( Int_t ifas = 0; ifas < m_iFast + 1 ; ifas++ ) {
      fast(ifas,4) = 0.;
    }
    // Find the m_iFast highest momentum particles and
    // put the highest in fast[0], next in fast[1],....fast[m_iFast-1].
    // fast[m_iFast] is just a workspace.
    for ( Int_t i = 0; i < np; i++ ) {
      if ( pass == 2 ) {
	mom(i,4) = TMath::Sqrt( mom(i,1)*mom(i,1) 
				+ mom(i,2)*mom(i,2) ); 
      }
      for ( Int_t ifas = m_iFast - 1; ifas > -1; ifas-- ) {
	if ( mom(i,4) > fast(ifas,4) ) {
	  for ( Int_t j = 1; j < 6; j++ ) {
	    fast(ifas+1,j) = fast(ifas,j);
	    if ( ifas == 0 ) fast(ifas,j) = mom(i,j);	    
	  }
	} else {
	  for ( Int_t j = 1; j < 6; j++ ) {
	    fast(ifas+1,j) = mom(i,j);
	  }
	  break;
	}
      }
    }
    // Find axis with highest thrust (case 1)/ highest major (case 2).
    for ( Int_t ie = 0; ie < work.GetNrows(); ie++ ) {
      work(ie,4) = 0.;
    }
    Int_t p = TMath::Min( m_iFast, np ) - 1;
    // Don't trust Math.pow to give right answer always.
    // Want nc = 2**p.
    Int_t nc = iPow(2,p); 
    for ( Int_t n = 0; n < nc; n++ ) {
      for ( Int_t j = 1; j < 4; j++ ) {
	tdi[j] = 0.;
      }
      for ( Int_t i = 0; i < TMath::Min(m_iFast,n); i++ ) {
	sgn = fast(i,5);
	if (iPow(2,(i+1))*((n+iPow(2,i))/iPow(2,(i+1))) >= i+1){
	  sgn = -sgn;
	}
	for ( Int_t j = 1; j < 5-pass; j++ ) {
	  tdi[j] = tdi[j] + sgn*fast(i,j);
	}
      }
      tds = tdi[1]*tdi[1] + tdi[2]*tdi[2] + tdi[3]*tdi[3];
      for ( Int_t iw = TMath::Min(n,9); iw > -1; iw--) {
	if( tds > work(iw,4) ) {
	  for ( Int_t j = 1; j < 5; j++ ) {
	    work(iw+1,j) = work(iw,j);
	    if ( iw == 0 ) {
	      if ( j < 4 ) {
		work(iw,j) = tdi[j];
	      } else {
		work(iw,j) = tds;
	      }
	    }
	  }
	} else {
	  for ( Int_t j = 1; j < 4; j++ ) {
	    work(iw+1,j) = tdi[j];
	  }
	  work(iw+1,4) = tds;
	}
      }
    }
    // Iterate direction of axis until stable maximum.
    m_dThrust[pass] = 0;
    thp = -99999.;
    Int_t nagree = 0;
    for ( Int_t iw = 0; 
	  iw < TMath::Min(nc,10) && nagree < m_iGood; iw++ ){
      thp = 0.;
      thps = -99999.;
      while ( thp > thps + m_dConv ) {
	thps = thp;
	for ( Int_t j = 1; j < 4; j++ ) {
	  if ( thp <= 1E-10 ) {
	    tdi[j] = work(iw,j);
	  } else {
	    tdi[j] = tpr[j];
	    tpr[j] = 0;
	  }
	}
	for ( Int_t i = 0; i < np; i++ ) {
	  sgn = sign(mom(i,5), 
		     tdi[1]*mom(i,1) + 
		     tdi[2]*mom(i,2) + 
		     tdi[3]*mom(i,3));
	  for ( Int_t j = 1; j < 5 - pass; j++ ){
	    tpr[j] = tpr[j] + sgn*mom(i,j);
	  }
	}
	thp = TMath::Sqrt(tpr[1]*tpr[1] 
			  + tpr[2]*tpr[2] 
			  + tpr[3]*tpr[3])/tmax;
      }
      // Save good axis. Try new initial axis until enough
      // tries agree.
      if ( thp < m_dThrust[pass] - m_dConv ) {
	break;
      }
      if ( thp > m_dThrust[pass] + m_dConv ) {
	nagree = 0;
	sgn = iPow( -1, (Int_t)TMath::Nint(m_random.Rndm()) );
	for ( Int_t j = 1; j < 4; j++ ) {
	  m_dAxes(pass,j) = sgn*tpr[j]/(tmax*thp);
	}
	m_dThrust[pass] = thp;
      }
      nagree = nagree + 1;
    }
  }
  // Find minor axis and value by orthogonality.
  sgn = iPow( -1, (Int_t)TMath::Nint(m_random.Rndm()));
  m_dAxes(3,1) = -sgn*m_dAxes(2,2);
  m_dAxes(3,2) = sgn*m_dAxes(2,1);
  m_dAxes(3,3) = 0.;
  thp = 0.;
  for ( Int_t i = 0; i < np; i++ ) {
    thp += mom(i,5)*TMath::Abs(m_dAxes(3,1)*mom(i,1) + 
			       m_dAxes(3,2)*mom(i,2));
  }
  m_dThrust[3] = thp/tmax;
  // Rotate back to original coordinate system.
  for ( Int_t i6 = 0; i6 < 3; i6++ ) {
    for ( Int_t j = 1; j < 4; j++ ) {
      temp(i6,j) = m_dAxes(i6+1,j);
    }
    temp(i6,4) = 0;
  }
  ludbrb(&temp,the,phi,0.,0.,0.);
  for ( Int_t i7 = 0; i7 < 3; i7++ ) {
    for ( Int_t j = 1; j < 4; j++ ) {
      m_dAxes(i7+1,j) = temp(i7,j);
    }
  }
  m_dOblateness = m_dThrust[2] - m_dThrust[3];

  event_thrust = m_dThrust[1];
  event_oblateness = m_dOblateness;
  

  // sphericity and aplanarity
  double Pduzi; 
  float sp[3][3];
  float norm=0.0;
  double dvojka=2.0;
    
  for(int i=0;i<3;i++){
    for (int j=0;j<3;j++)       
      {     
	sp[j][i]=0.0;
      }
  }
  
  for(Int_t elem=0;elem<numElements;elem++) {
    o = e->At(elem);
    double pp [] = {((TVector3 *) o)->X(),((TVector3 *) o)->Y(),((TVector3 *) o)->Z()};

    Pduzi=sqrt(pow(pp[0],dvojka)+pow(pp[1],dvojka)+pow(pp[2],dvojka));
      
    norm=norm+pow(Pduzi,(double) 2.0); // _r: exponent in sphericity tensor use 2.0 for classical 1.0 for C,D
    for(int j=0;j<3;j++){
      for(int i=0;i<3;i++){
	sp[j][i]=sp[j][i]+pp[i]*pp[j]*pow(Pduzi,(2.0-dvojka));
      }	
    }    
  }
  
  for(int j=0;j<3;j++){
    for(int i=0;i<3;i++){
      sp[j][i]=sp[j][i]/norm;
      //  cout << "sp tenzor "<< sp[j][i] << endl;
    }	
  }
  float lre[3];
  JAMMA::Eigenvalue test(sp);
  test.getRealEigenvalues(lre);
  // cout << lre[0] << "   " << lre[1] << "  " << lre[2] << endl;  
  float sphericity;
  float aplanarity;
  sphericity=1.5*(lre[0]+lre[1]);
  aplanarity=1.5*lre[0];

  event_sphericity = sphericity;
  event_aplanarity = aplanarity;

}

Double_t ulAngle(Double_t x, Double_t y)
{
  Double_t ulangl = 0;
  Double_t r = TMath::Sqrt(x*x + y*y);
  if ( r < 1.0E-20 ) {
    return ulangl; 
  }
  if ( TMath::Abs(x)/r < 0.8 ) {
    ulangl = sign(TMath::ACos(x/r),y);
  } else {
    ulangl = TMath::ASin(y/r);
    if ( x < 0. && ulangl >= 0. ) {
      ulangl = TMath::Pi() - ulangl;
    } else if ( x < 0. ) {
      ulangl = -TMath::Pi() - ulangl;
    }
  }
  return ulangl;
}

Double_t sign(Double_t a, Double_t b) {
  if ( b < 0 ) {
    return -TMath::Abs(a);
  } else {
    return TMath::Abs(a);
  }
}

void ludbrb(TMatrixD* mom, Double_t the, Double_t phi, Double_t bx, Double_t by, Double_t bz) {
  // Ignore "zeroth" elements in rot,pr,dp.
  // Trying to use physics-like notation.
  TMatrixD rot(4,4);
  Double_t pr[4];
  Double_t dp[5];
  Int_t np = mom->GetNrows();
  if ( the*the + phi*phi > 1.0E-20 )
    {
      rot(1,1) = TMath::Cos(the)*TMath::Cos(phi);
      rot(1,2) = -TMath::Sin(phi);
      rot(1,3) = TMath::Sin(the)*TMath::Cos(phi);
      rot(2,1) = TMath::Cos(the)*TMath::Sin(phi);
      rot(2,2) = TMath::Cos(phi);
      rot(2,3) = TMath::Sin(the)*TMath::Sin(phi);
      rot(3,1) = -TMath::Sin(the);
      rot(3,2) = 0.0;
      rot(3,3) = TMath::Cos(the);
      for ( Int_t i = 0; i < np; i++ ) {
        for ( Int_t j = 1; j < 4; j++ ) {
          pr[j] = (*mom)(i,j);
          (*mom)(i,j) = 0;
        }
        for ( Int_t jb = 1; jb < 4; jb++) {
          for ( Int_t k = 1; k < 4; k++) {
            (*mom)(i,jb) = (*mom)(i,jb) + rot(jb,k)*pr[k];
          }
        }
      }
      Double_t beta = TMath::Sqrt( bx*bx + by*by + bz*bz );
      if ( beta*beta > 1.0E-20 ) {
        if ( beta >  0.99999999 ) {
          //send message: boost too large, resetting to <~1.0.
          bx = bx*(0.99999999/beta);
          by = by*(0.99999999/beta);
          bz = bz*(0.99999999/beta);
          beta =   0.99999999;
        }
        Double_t gamma = 1.0/TMath::Sqrt(1.0 - beta*beta);
        for ( Int_t i = 0; i < np; i++ ) {
          for ( Int_t j = 1; j < 5; j++ ) {
            dp[j] = (*mom)(i,j);
          }
          Double_t bp = bx*dp[1] + by*dp[2] + bz*dp[3];
          Double_t gbp = gamma*(gamma*bp/(1.0 + gamma) + dp[4]);
          (*mom)(i,1) = dp[1] + gbp*bx;
          (*mom)(i,2) = dp[2] + gbp*by;
          (*mom)(i,3) = dp[3] + gbp*bz;
          (*mom)(i,4) = gamma*(dp[4] + bp);
        }
      }
    }
  return;
}

Int_t iPow(Int_t man, Int_t exp) {
  Int_t ans = 1;
  for( Int_t k = 0; k < exp; k++) {
    ans = ans*man;
  }
  return ans;
}

EVENT::ReconstructedParticle* getReconstructedLepton(EVENT::LCEvent* event){
  // to return the reconstructed signal lepton
  Int_t generator_lepton_id = 0;
  Int_t mc_lepton_id = 0;
  Int_t linked_mc_lepton_id = 0;

  TLorentzVector generator_lepton;
  TLorentzVector mc_lepton;
  TLorentzVector linked_mc_lepton;
  TLorentzVector linked_rec_lepton;

  EVENT::LCRelation* link_to_rec_lepton;
  EVENT::MCParticle* pointer_to_mc_lepton = 0;
  EVENT::ReconstructedParticle* pointer_to_rec_lepton = 0;

  // find generated lepton status 3
  EVENT::LCCollection* mcp = event->getCollection("MCParticle");
  for (int i = 0; i < mcp->getNumberOfElements(); i++){
    EVENT::MCParticle* mcpi = (EVENT::MCParticle*) mcp->getElementAt(i);
    if (mcpi->getGeneratorStatus() == 3){
      if (abs(mcpi->getPDG()) == 11 || abs(mcpi->getPDG()) == 13) {
        if(mcpi->getParents().size() > 0){
          if (abs(mcpi->getParents()[0]->getPDG()) == 24){
            generator_lepton = TLorentzVector(mcpi->getMomentum());
            generator_lepton_id = mcpi->getPDG();
          }
        }
      }
    }
  }
  //cout << "got generated lepton " << generator_lepton_id << " " << generator_lepton.Theta() << " " << generator_lepton.Pt() << endl;
  
  // find generated lepton status 1 with lepton ID and W parent
  for (int i = 0; i < mcp->getNumberOfElements(); i++){
    EVENT::MCParticle* mcpi = (EVENT::MCParticle*) mcp->getElementAt(i);
    if (mcpi->getPDG() == generator_lepton_id){
      if (mcpi->getGeneratorStatus() == 1){
        if (mcpi->getParents().size() == 1){
          if (abs(mcpi->getParents()[0]->getPDG()) == 24){
            pointer_to_mc_lepton = mcpi;
            mc_lepton = TLorentzVector(mcpi->getMomentum());
            mc_lepton_id = mcpi->getPDG();
          }
        }
      }
    }
  }
  //cout << "got mc lepton " << mc_lepton_id << " " << mc_lepton.Theta() << " " << mc_lepton.Pt() << endl;

  // find link which points to MC lepton
  EVENT::LCCollection* links = event->getCollection("RecoMCTruthLink");

  for (int i = 0; i < links->getNumberOfElements(); i++){
    EVENT::LCRelation* linki = (EVENT::LCRelation*) links->getElementAt(i);
    EVENT::MCParticle* mcpi = (EVENT::MCParticle*) linki->getTo();
    if (pointer_to_mc_lepton == mcpi){
      linked_mc_lepton = TLorentzVector(mcpi->getMomentum());
      linked_mc_lepton_id = mcpi->getPDG();
      link_to_rec_lepton = linki;
    }
  }
  //cout << "got linked mc lepton " << linked_mc_lepton_id << " " << linked_mc_lepton.Theta() << " " << linked_mc_lepton.Pt() << endl;

  // follow the link
  // if no link, 0 returned
  if (linked_mc_lepton.P() != 0){
    EVENT::ReconstructedParticle* rpi = (EVENT::ReconstructedParticle*) link_to_rec_lepton->getFrom();
    pointer_to_rec_lepton = rpi;
    linked_rec_lepton = TLorentzVector(rpi->getMomentum());
    //cout << "got linked rec lepton " << rpi->getType() << " " << linked_rec_lepton.Theta() << " " << linked_rec_lepton.Pt() << endl;
  }

  return pointer_to_rec_lepton;
}
