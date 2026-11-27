#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenGlauberParam.h"
#include <fairlogger/Logger.h>
#endif

bool setNuclearShape(NA6PNuclearShape& nuc, std::string name)
{
  if (strcmp(name.c_str(), "Pb") == 0) {
    nuc.setNucleusPb();
  } else if (strcmp(name.c_str(), "C") == 0) {
    nuc.setNucleusC();
  } else if (strcmp(name.c_str(), "O") == 0) {
    nuc.setNucleusO();
  } else if (strcmp(name.c_str(), "Al") == 0) {
    nuc.setNucleusAl();
  } else if (strcmp(name.c_str(), "In") == 0) {
    nuc.setNucleusIn();
  } else if (strcmp(name.c_str(), "p") == 0) {
    nuc.setNucleusFromParametrization(1, 1);
  } else {
    LOGP(fatal, "Nucleus {} not supported", name.c_str());
    return false;
  }
  return true;
}

NA6PGenerator* genGlauberNA49(std::string nucProj = "Pb",
                              std::string nucTarg = "Pb",
                              double bmin = 0.,
                              double bmax = 18.)
{
  NA6PNuclearShape nucA, nucB;
  bool setA = setNuclearShape(nucA, nucProj);
  bool setB = setNuclearShape(nucB, nucTarg);
  if (!setA || !setB)
    return nullptr;
  float sigmb = 30.;
  NA6PGlauber gl(nucA, nucB, sigmb);
  LOGP(info, "Projectile A = {} Z = {}", nucA.getA(), nucA.getZ());
  LOGP(info, "Target     A = {} Z = {}", nucB.getA(), nucB.getZ());
  LOGP(info, "Impact parameter range: [{}, {}] fm", bmin, bmax);
  NA6PGenGlauberParam* genG = new NA6PGenGlauberParam(gl, bmin, bmax);
  return genG;
}
