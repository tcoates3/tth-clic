
#include <TFile.h>
#include <TMatrixD.h>
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCEvent.h"

Int_t treeMaker(std::string fileName, TString outfileName);
Int_t treeMaker_hadronic(std::string fileName, TString outfileName);

EVENT::MCParticle* getGeneratedLepton(EVENT::LCEvent* event);
EVENT::MCParticle* getMCLepton(EVENT::LCEvent* event);
EVENT::ReconstructedParticle* getReconstructedLepton(EVENT::LCEvent* event);

void CalculateEventShapeVariables(TObjArray* e, Double_t &event_thrust, Double_t &event_oblateness, Double_t &event_sphericity , Double_t &event_aplanarity );
Double_t ulAngle(Double_t x, Double_t y);
Double_t sign(Double_t a, Double_t b);
void ludbrb(TMatrixD* mom, Double_t the, Double_t phi, Double_t bx, Double_t by, Double_t bz);
Int_t iPow(Int_t man, Int_t exp);
bool compare(const std::pair<Double_t,Double_t>&i, const std::pair<Double_t,Double_t>&j);


