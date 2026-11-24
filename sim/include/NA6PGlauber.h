// NA6PCCopyright

// Based on:
// Copyright 2019-2030 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_GLAUBER_H
#define NA6P_GLAUBER_H

#include "NA6PNuclearShape.h"

class NA6PGlauber : public TObject
{
 public:
  struct GlauberEvent {
    float impactParameter;
    float phiImpPar;
    int nParticipants;
    int nSpectatorsA;
    int nSpectatorsB;
    int nCollisions;
    float eccentricity;
    std::vector<float> xPart, yPart;
    std::vector<float> xColl, yColl;
    void clear()
    {
      xPart.clear();
      yPart.clear();
      xColl.clear();
      yColl.clear();
    }
  };

  NA6PGlauber();
  NA6PGlauber(const NA6PNuclearShape& nucA, const NA6PNuclearShape& nucB, double sigmb);
  virtual ~NA6PGlauber() = default;

  void setSigmamb(double sigmb) { mSigmamb = sigmb; }
  void setSigmafm2(double sigfm2) { mSigmamb = 10. * sigfm2; } // conversion from mb(=10^-31 m2) to fm2 (=10-30 m2) -> 10 mb=1 fm2
  void setNucleusA(const NA6PNuclearShape& nuc) { mNuclA = nuc; }
  void setNucleusB(const NA6PNuclearShape& nuc) { mNuclB = nuc; }
  void setImpactParameter(double b) { mImpPar = b; }
  void setIntegrationSteps(double dx, double dy, double dz);

  double getSigmamb() const { return mSigmamb; }
  double getSigmafm2() const { return mSigmamb / 10.; }
  double getImpactParam() const { return mImpPar; }
  const NA6PNuclearShape& getNucleusA() const { return mNuclA; }
  const NA6PNuclearShape& getNucleusB() const { return mNuclB; }
  double getTAB() const { return mTAB; }
  double getTa() const { return mTa; }
  double getTb() const { return mTb; }
  double getTaTb() const { return mTaTb; }
  double getInteractionProbability() const
  {
    double pnin = std::pow((1. - mTAB * getSigmafm2()), mNuclA.getA() * mNuclB.getA());
    return 1. - pnin;
  }
  double getMeanX() const { return mX / mTaTb; }
  double getMeanY() const { return mY / mTaTb; }
  double getVarX() const
  {
    double mX = getMeanX();
    return mXX / mTaTb - mX * mX;
  }
  double getVarY() const
  {
    double mY = getMeanY();
    return mYY / mTaTb - mY * mY;
  }
  double getEccentricity() const;
  double getAveragePartDens() const
  {
    if (mNuclA.getA() == 1 || mNuclB.getA() == 1 || mTaTb <= 0.)
      return -999.;
    return mPartDens / mTaTb;
  }
  double getMaximumPartDens() const
  {
    if (mNuclA.getA() == 1 || mNuclB.getA() == 1 || mTaTb <= 0.)
      return -999.;
    return mPartDensMax;
  }
  double getAverageCollDens() const
  {
    if (mNuclA.getA() == 1 || mNuclB.getA() == 1 || mTaTb <= 0.)
      return -999.;
    return mCollDens / mTaTb;
  }
  double getMaximumCollDens() const
  {
    if (mNuclA.getA() == 1 || mNuclB.getA() == 1 || mTaTb <= 0.)
      return -999.;
    return mCollDensMax;
  }
  double getAverageNPart() const { return mAverNPart; }
  double getAverageNPartInInteractions() const
  {
    if (getInteractionProbability() > 0.)
      return mAverNPart / getInteractionProbability();
    else
      return 0.;
  }
  double getAverageNColl() const { return mNuclA.getA() * mNuclB.getA() * getSigmafm2() * mTAB; }
  double getRhoL() const { return mRhoL / mTaTb; }

  void doOpticalIntegrals(double b);
  double calcPartDens(double tA, double tB) const;
  double calcPartDensFromCoord(double x, double y, double b) const;
  double calcCollDens(double tA, double tB) const;
  double calcCollDensFromCoord(double x, double y, double b) const;
  bool isInside(double x, double y, double b) const;

  GlauberEvent generateEvent(double bmin = -1., double bmax = -1., bool storePositions = false);

 private:
  NA6PNuclearShape mNuclA;
  NA6PNuclearShape mNuclB;
  double mSigmamb = 30.;                  // NN c.s. in mb
  double mIntegStep[3] = {0.3, 0.3, 0.3}; // Integration step in x,y,z
  double mImpPar = 9999.;                 // impact parameter (fm) (directed along y)
  double mPartDens = 0.;
  double mPartDensMax = 0.;
  double mCollDens = 0.;
  double mCollDensMax = 0.;
  double mTa = 0.;
  double mTb = 0.;
  double mX = 0.;
  double mY = 0.;
  double mXX = 0.;
  double mYY = 0.;
  double mTAB = 0.;
  double mTaTb = 0.;
  double mAverNPart = 0.;
  double mRhoL = 0.;

  void doOpticalIntegralsPA(double b);

  ClassDef(NA6PGlauber, 1)
};

#endif
