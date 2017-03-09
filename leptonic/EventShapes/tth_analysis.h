
#include <TFile.h>
#include <TMatrixD.h>
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCEvent.h"

void eventShapes(std::string fileName, TString outfileName, std::string LeptonCollection, std::string JetCollection, std::string PFOCollection, std::string PFOsInJetsCollection);

void CalculateEventShapeVariables(TObjArray* e, Double_t &event_thrust, Double_t &event_oblateness, Double_t &event_sphericity , Double_t &event_aplanarity );
Double_t ulAngle(Double_t x, Double_t y);
Double_t sign(Double_t a, Double_t b);
void ludbrb(TMatrixD* mom, Double_t the, Double_t phi, Double_t bx, Double_t by, Double_t bz);
Int_t iPow(Int_t man, Int_t exp);

EVENT::ReconstructedParticle* getReconstructedLepton(EVENT::LCEvent* event);


