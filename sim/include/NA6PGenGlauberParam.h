// NA6PCCopyright

#ifndef NA6P_GEN_GLAUBER_PARAM_H
#define NA6P_GEN_GLAUBER_PARAM_H

#include <TF1.h>
#include "NA6PGenerator.h"
#include "NA6PGlauber.h"
#include "NA6PGenParam.h"

class NA6PGenGlauberParam : public NA6PGenerator
{
 public:
  using NA6PGenerator::NA6PGenerator;
  NA6PGenGlauberParam() = default;
  NA6PGenGlauberParam(const NA6PGlauber& gl, float bmin = 0., float bmax = 20.);
  ~NA6PGenGlauberParam() override = default;
  void setKeepOnlyEventsWithInteraction(bool opt) { mKeepOnlyInteracting = opt; }
  void generate() override;

  void init() override;

 protected:
  static constexpr int kMaxTrials = 1000;
  NA6PGlauber mGlauberMC;           // Glauber helper class
  float mImpParMin = 0.f;           // minimum impact parameter
  float mImpParMax = 20.f;          // maximum impact parameter
  float mPtMin = 0.f;               // minimum pt
  float mPtMax = 10.f;              // maximum pt
  float mYMin = 0.5f;               // minimum y
  float mYMax = 5.f;                // maximum y
  float mNParticipantsRef = 350.f;  // n. of participants in reference events
  bool mKeepOnlyInteracting = true; // keep only events with >=1 NN collision
  NA6PGenParam mGenPiMinus;         // parametrized generator for pi-
  NA6PGenParam mGenPiPlus;          // parametrized generator for pi+
  NA6PGenParam mGenKMinus;          // parametrized generator for K-
  NA6PGenParam mGenKPlus;           // parametrized generator for K+
  NA6PGenParam mGenProtons;         // parametrized generator for p

  bool configureParametrizations();
  int generateSpecies(const NA6PGenParam& genp, int nParticipants);
  void generateProjectileSpectators(int nSpec, int zProj, int aProj);
  ClassDefOverride(NA6PGenGlauberParam, 1);
};

#endif
