#pragma once


#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/core/include/Hists.h"

#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

#include <vector>

/**
 *   Example class for booking and filling histograms, the new version using AnalysisModule mechanisms.
 */

struct singleHists
{
  TH1F *pt, *eta, *phi, *mass, *charge;
  TH2F* pt_eta;
};


struct GenParticleHists
{ 
  TH1F* number, *decay_mom, *decay_daughter, *deltaRmin, *nextParticle, *deltaR_daughters;
  std::vector<singleHists> stdHists;
  

};


class VLQGenHists: public uhh2::Hists {
 public:
  // use the same constructor arguments as Hists for forwarding:
  VLQGenHists(uhh2::Context & ctx, const std::string & dirname);
  
  virtual void fill(const uhh2::Event & ev) override;
  virtual ~VLQGenHists();
  
 private:
  
  GenParticleHists histoBooker(const std::string& HistName, double minMass, double maxMass);
  void histoFiller(std::vector<GenParticle> & particles, int partNumber, double weight);
  void decayFiller(double weight, int partNumber, int mother1, int mother2, int daughter1, int daughter2);
  void fill_nearest(int position, double weight, double deltaRmin, double pdgId);
  void fillDaughterDistance(int position, double deltaR, double weight);

  int positionHelper(const std::string& Name);

  std::vector<GenParticleHists> m_Hists;
  std::vector<std::string> PartNames;
  std::vector<int> PartPdgId;

  TH1F* particles_noMother, *particles_noMother_pT, *particles_noMother_eta, *particles_noMother_phi;
  TH1F* particles_noMotherNoDaughter, *particles_noMotherNoDaughter_pT, *particles_noMotherNoDaughter_eta, *particles_noMotherNoDaughter_phi;

  TH1F* VLQ_mother;
  TH2F* VLQ_mother1_mother2;

  TH1F* deltaR_w, *deltaPhi_w;
};
