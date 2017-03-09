
#include <TFile.h>
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCEvent.h"

std::vector<EVENT::MCParticle*> getMCTaus(EVENT::LCEvent* event);
std::vector<EVENT::MCParticle*> getMCTauDaughters(EVENT::MCParticle* mcpi);
std::vector<EVENT::MCParticle*>  getAllMCTauDaughters(EVENT::MCParticle* mcpi, std::vector<EVENT::MCParticle*> daughters);

std::vector<EVENT::ReconstructedParticle*> getRecoTaus(EVENT::LCEvent* event);
std::vector<EVENT::ReconstructedParticle*> getRecoTauDaughters(EVENT::ReconstructedParticle* taui);
std::vector<EVENT::ReconstructedParticle*> getAllRecoTauDaughters(EVENT::ReconstructedParticle* taui, std::vector<EVENT::ReconstructedParticle*> daughters);

Int_t leptonInvest(std::string recfileName, TString outfileName);

void postLeptonFinder(std::string recfileName, TString outfileName);

EVENT::MCParticle* getGeneratedLepton(EVENT::LCEvent* event);
EVENT::MCParticle* getMCLepton(EVENT::LCEvent* event);
EVENT::ReconstructedParticle* getReconstructedLepton(EVENT::LCEvent* event);
