// NA6PCCopyright

#include "NA6PGenGlauberParam.h"
#include "NA6PBeamParam.h"
#include "NA6PGenCutParam.h"
#include "NA6PMCStack.h"
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <cmath>

NA6PGenGlauberParam::NA6PGenGlauberParam(const NA6PGlauber& gl, float bmin, float bmax) : NA6PGenerator("GenGlauber"),
                                                                                          mGlauberMC(gl),
                                                                                          mImpParMin(bmin),
                                                                                          mImpParMax(bmax)
{
}

void NA6PGenGlauberParam::init()
{
  if (isInitDone()) {
    return;
  }
  configureParametrizations();
  mGenPiMinus.initMultiplicity();
  mGenPiPlus.initMultiplicity();
  mGenKMinus.initMultiplicity();
  mGenKPlus.initMultiplicity();
  mGenProtons.initMultiplicity();
  NA6PGenerator::init();
}

bool NA6PGenGlauberParam::configureParametrizations()
{
  float T_piM, y0_piM, ysig_piM, dNdY_piM;
  float T_piP, y0_piP, ysig_piP, dNdY_piP;
  float T_KM, y0_KM, ysig_KM, dNdY_KM;
  float T_KP, y0_KP, ysig_KP, dNdY_KP;
  float T_Prot, y0_Prot, ysig_Prot, dNdY_Prot;

  float en = NA6PBeamParam::Instance().energyPerNucleon;
  // pions and Kaons from  NA49 nucl-ex/0205002

  if (std::abs(en - 40.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 40. GeV/c");
    T_piM = T_piP = 0.169;
    T_KM = 0.232;
    T_KP = 0.226;
    T_Prot = 0.25;

    y0_piM = y0_piP = 0.666; // pions 2-gaussian poles half distance
    y0_KM = 0.694;           // K- 2-gaussian poles half distance
    y0_KP = 0.694;           // K+ 2-gaussian poles half distance
    y0_Prot = 0.907;         // ???

    ysig_piM = ysig_piP = 0.872;
    ysig_KM = 0.635;
    ysig_KP = 0.725;
    ysig_Prot = 0.798;

    dNdY_piM = 106.1;
    dNdY_piP = 96.6;
    dNdY_KM = 7.58;
    dNdY_KP = 20.1;
    dNdY_Prot = 39.9;
  } else if (std::abs(en - 158.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 158. GeV/c");
    T_piM = T_piP = 0.18;
    T_KM = 0.226;
    T_KP = 0.232;
    T_Prot = 0.31;

    y0_piM = y0_piP = 0.72; // pions 2-gaussian poles half distance
    y0_KM = 0.727;          // K- 2-gaussian poles half distance
    y0_KP = 0.839;          // K+ 2-gaussian poles half distance
    y0_Prot = 39.8;         // ???

    ysig_piM = ysig_piP = 1.18;
    ysig_KM = 0.81;
    ysig_KP = 0.88;
    ysig_Prot = 8.07;

    dNdY_piM = 175.4;
    dNdY_piP = 170.1;
    dNdY_KM = 16.8;
    dNdY_KP = 29.6;
    dNdY_Prot = 27.;
  } else {
    LOGP(fatal, "No parameters defined for beam energy {}", en);
    return false;
  }
  mNParticipantsRef = 350.; // Pb-Pb 0-10%
  const float ycm = NA6PBeamParam::Instance().getYCM();

  std::string dndptFun = "x*exp(-sqrt(x*x+[0]*[0])/[1])";
  std::string dndyFun = "exp(-0.5*pow((x-[0]-[1])/[2],2))+exp(-0.5*pow((x-[0]+[1])/[2],2))";

  mGenPiMinus.setPDGCode(-211);
  mGenPiMinus.setMultiplicity(0);
  mGenPiMinus.setPoisson(true);
  mGenPiMinus.setLongIsY(true);
  mGenPiMinus.setTransIsPt(true);
  mGenPiMinus.setFunTrans(dndptFun, mPtMin, mPtMax);
  mGenPiMinus.setFunLong(dndyFun, mYMin, mYMax);
  mGenPiMinus.setParametersTrans({0.1396f, T_piM});
  mGenPiMinus.setParametersLong({ycm, y0_piM, ysig_piM});
  mGenPiMinus.setdNdY(dNdY_piM);

  mGenPiPlus.setPDGCode(211);
  mGenPiPlus.setMultiplicity(0);
  mGenPiPlus.setPoisson(true);
  mGenPiPlus.setLongIsY(true);
  mGenPiPlus.setTransIsPt(true);
  mGenPiPlus.setFunTrans(dndptFun, mPtMin, mPtMax);
  mGenPiPlus.setFunLong(dndyFun, mYMin, mYMax);
  mGenPiPlus.setParametersTrans({0.1396f, T_piP});
  mGenPiPlus.setParametersLong({ycm, y0_piP, ysig_piP});
  mGenPiPlus.setdNdY(dNdY_piP);

  mGenKMinus.setPDGCode(-321);
  mGenKMinus.setMultiplicity(0);
  mGenKMinus.setPoisson(true);
  mGenKMinus.setLongIsY(true);
  mGenKMinus.setTransIsPt(true);
  mGenKMinus.setFunTrans(dndptFun, mPtMin, mPtMax);
  mGenKMinus.setFunLong(dndyFun, mYMin, mYMax);
  mGenKMinus.setParametersTrans({0.4937f, T_KM});
  mGenKMinus.setParametersLong({ycm, y0_KM, ysig_KM});
  mGenKMinus.setdNdY(dNdY_KM);

  mGenKPlus.setPDGCode(321);
  mGenKPlus.setMultiplicity(0);
  mGenKPlus.setPoisson(true);
  mGenKPlus.setLongIsY(true);
  mGenKPlus.setTransIsPt(true);
  mGenKPlus.setFunTrans(dndptFun, mPtMin, mPtMax);
  mGenKPlus.setFunLong(dndyFun, mYMin, mYMax);
  mGenKPlus.setParametersTrans({0.4937f, T_KP});
  mGenKPlus.setParametersLong({ycm, y0_KP, ysig_KP});
  mGenKPlus.setdNdY(dNdY_KP);

  mGenProtons.setPDGCode(2212);
  mGenProtons.setMultiplicity(0);
  mGenProtons.setPoisson(true);
  mGenProtons.setLongIsY(true);
  mGenProtons.setTransIsPt(true);
  mGenProtons.setFunTrans(dndptFun, mPtMin, mPtMax);
  mGenProtons.setFunLong(dndyFun, mYMin, mYMax);
  mGenProtons.setParametersTrans({0.9383f, T_Prot});
  mGenProtons.setParametersLong({ycm, y0_Prot, ysig_Prot});
  mGenProtons.setdNdY(dNdY_Prot);

  return true;
}

void NA6PGenGlauberParam::generate()
{
  generatePrimaryVertex(); // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack.
  NA6PGlauber::GlauberEvent ev = mGlauberMC.generateEvent(mImpParMin, mImpParMax);
  int nCollisions = ev.nCollisions;
  int ntry = 1;
  if (mKeepOnlyInteracting && nCollisions == 0) {
    while (nCollisions == 0 && ntry < kMaxTrials) {
      ev = mGlauberMC.generateEvent(mImpParMin, mImpParMax);
      ++ntry;
      nCollisions = ev.nCollisions;
    }
    if (nCollisions == 0) {
      LOGP(error, "No interacting event in {} < b < {} after {} trials", mImpParMin, mImpParMax, kMaxTrials);
    }
  }
  float b = ev.impactParameter;
  int nParticipants = ev.nParticipants;
  int nSpectA = ev.nSpectatorsA;
  int nGenerated = 0;
  nGenerated += generateSpecies(mGenPiMinus, nParticipants);
  nGenerated += generateSpecies(mGenPiPlus, nParticipants);
  nGenerated += generateSpecies(mGenKMinus, nParticipants);
  nGenerated += generateSpecies(mGenKPlus, nParticipants);
  nGenerated += generateSpecies(mGenProtons, nParticipants);
  int nSpecGenerated = 0;
  if (nCollisions > 0) {
    generateProjectileSpectators(nSpectA, mGlauberMC.getNucleusA().getZ(), mGlauberMC.getNucleusA().getA());
    nSpecGenerated = nSpectA;
  }
  // register generator header in the MCHeader
  auto mcHead = getStack()->getEventHeader();
  mcHead->getGenHeaders().emplace_back(nGenerated, 0, mcHead->getNPrimaries(), 0, "GenGlauber_ProducedParticles");
  mcHead->incNPrimaries(nGenerated);
  mcHead->getGenHeaders().emplace_back(nSpecGenerated, 0, mcHead->getNPrimaries(), 0, "GenGlauber_Spectators");
  mcHead->incNPrimaries(nSpecGenerated);
  getStack()->setImpPar(b);
  getStack()->setNPart(nParticipants);
  getStack()->setNColl(nCollisions);
}

int NA6PGenGlauberParam::generateSpecies(const NA6PGenParam& genp, int nParticipants)
{
  // mostly duplicated from NA6PGenParam
  const auto& param = NA6PGenCutParam::Instance();
  int pdgCode = genp.getPDGCode();
  auto m = TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass();
  float nParticles = genp.getMultiplicity() * (float)nParticipants / mNParticipantsRef;
  int nToGenerate = genp.isPoisson() ? gRandom->Poisson(nParticles) : TMath::Nint(nParticles);
  LOGP(info, "Generate {} particles of species {} in event with {} participant nucleons", nToGenerate, pdgCode, nParticipants);
  for (int i = 0; i < nToGenerate; i++) {
    int ntrial = 0;
    do {
      double vTrans, vLong, phi, pt2, pt, mt2, mt, pZ, en;
      try {
        vTrans = genp.getFunTrans()->GetRandom();
        vLong = genp.getFunLong()->GetRandom();
        phi = gRandom->Uniform(0, TMath::TwoPi());
        pt2 = vTrans * vTrans;
        mt2 = pt2;
        pt = vTrans;
        mt = pt;
        if (!genp.isTransPt()) { // mT
          pt2 -= m * m;
          pt = std::sqrt(pt2);
        } else {
          mt2 += m * m;
          mt = std::sqrt(mt2);
        }
        if (genp.isLongY()) {
          pZ = mt * std::sinh(vLong);
        } else { // eta
          auto tht = 2. * std::atan(std::exp(-vLong));
          pZ = pt / std::tan(tht);
        }
        en = std::sqrt(mt2 + pZ * pZ);
      } catch (std::exception const& e) {
        LOGP(info, "Exception at trial {} of generating particle PDG={}, VTrans={}, vLong={}, error: {}", ntrial, pdgCode, vTrans, vLong, e.what());
        if (ntrial++ > param.maxTrialsPerParticle) {
          LOGP(fatal, "Failed to generate track PDG:{} after max allowed {} trials, check generator", pdgCode, ntrial);
        }
        continue;
      }
      ntrial = 0;
      int dummy = 0;
      float sn = std::sin(phi), cs = std::cos(phi);
      getStack()->PushTrack(true, -1, pdgCode, pt * cs, pt * sn, pZ, en, getOriginX(), getOriginY(), getOriginZ(),
                            0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., 0);
    } while (ntrial != 0);
  }
  return nToGenerate;
}

void NA6PGenGlauberParam::generateProjectileSpectators(int nSpec, int zProj, int aProj)
{
  double mP = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  double mN = TDatabasePDG::Instance()->GetParticle(2112)->Mass();
  double eBeam = NA6PBeamParam::Instance().energyPerNucleon;
  int zRemaining = zProj, aRemaining = aProj;
  int nProtons = 0, nNeutrons = 0;
  for (int i = 0; i < nSpec; ++i) {
    bool isProton = gRandom->Rndm() * aRemaining < zRemaining;
    --aRemaining;
    if (isProton)
      --zRemaining;
    int pdgCode = isProton ? 2212 : 2112;
    double m = isProton ? mP : mN;
    double pZ = std::sqrt(eBeam * eBeam - m * m);
    double en = eBeam;
    int dummy = 0;
    getStack()->PushTrack(true, -1, pdgCode, 0., 0., pZ, en, getOriginX(), getOriginY(), getOriginZ(),
                          0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., 0);
    if (isProton)
      ++nProtons;
    else
      ++nNeutrons;
  }
  LOGP(info, "Add {} spectator protons and {} spectator neutrons (total = {})", nProtons, nNeutrons, nSpec);
}
