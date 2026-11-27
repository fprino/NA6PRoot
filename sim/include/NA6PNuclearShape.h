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

#ifndef NA6PNUCLEARSHAPE_H
#define NA6PNUCLEARSHAPE_H

#include "TObject.h"
#include <span>

class TF1;

class NA6PNuclearShape : public TObject
{
 public:
  static constexpr int kMaxNucleons = 250;
  static constexpr double kRMin = 0.;
  static constexpr double kRMax = 20.;
  static constexpr int kMaxAttempts = 30;
  static constexpr double kPrecision = 1.e-8;

  NA6PNuclearShape();
  NA6PNuclearShape(int Z, int A);
  NA6PNuclearShape(int Z, int A, double r0, double wid, double w, int e);
  ~NA6PNuclearShape() = default;

  void setZ(int Z) { mZ = Z; }
  void setA(int A);
  void setShapeParameters(double r0, double wid, double w, int e);
  void setMinimumDistanceBetweenNucleons(double d) { mMinNucleonDist = d; }
  void setNucleusFromParametrization(int Z, int A);
  void setNucleusC();     // C-12      3-parameter Fermi model
  void setNucleusO();     // O-16      3-parameter Fermi model
  void setNucleusAl();    // Al-27     2-parameter Fermi model
  void setNucleusS();     // S-32      3-parameter Fermi model
  void setNucleusAr2pF(); // Ar-40     2-parameter Fermi model (Default)
  void setNucleusAr3pF(); // Ar-40     3-parameter Fermi model
  void setNucleusFe();    // Fe-56     2-parameter Fermi model
  void setNucleusCu2pF(); // Cu-63     2-parameter Fermi model (Default)
  void setNucleusCu3pG(); // Cu-63     3-parameter Gaussian model
  void setNucleusIn();    // In-115    2-parameter Fermi model
  void setNucleusAu();    // Au-197    2-parameter Fermi model
  void setNucleusPb();    // Pb-208    2-parameter Fermi model
  void setNucleusU();     // U-238     2-parameter Fermi model

  int getZ() const { return mZ; }
  int getA() const { return mA; }
  double getRho0() const { return mRho0; }
  double getRadius() const { return mRadius; }
  double getSkinDepth() const { return mSkinDepth; }
  double getW() const { return mW; }
  int getExponent() const { return mExpon; }

  std::span<const double> getXpos() const
  {
    return {mXpos.data(), static_cast<std::size_t>(mA)};
  }
  std::span<const double> getYpos() const
  {
    return {mYpos.data(), static_cast<std::size_t>(mA)};
  }
  std::span<const double> getZpos() const
  {
    return {mZpos.data(), static_cast<std::size_t>(mA)};
  }

  double calcRho(double r) const;
  double calcSurfNuclDens(double x, double y, double zstep) const;
  double calcSurfNuclDens(double x, double y, double zstep, double zm) const;
  double calcSurfNuclDens2(double x, double y, double zstep, double zm) const;
  double radialDensity(double r) const
  {
    return r * r * calcRho(r);
  }
  double integrateRho() const;
  double getMaxRadius() const { return mMaxRadius; }
  void normalizeRho0();

  TF1 getDensity() const;
  void resetNucleonPositions()
  {
    mXpos.fill(0.);
    mYpos.fill(0.);
    mZpos.fill(0.);
  }
  double sampleRadius(double rhoMax, double rMin, double rMax) const;
  void setNucleonPositions(double shiftx, double shifty);

 private:
  int mZ = 82;                  // atomic number
  int mA = 208;                 // mass number
  double mRho0 = 0.160;         // density at r=0
  double mRadius = 6.624;       // Wood-Saxon radius (r0)
  double mSkinDepth = 0.549;    // Wood-Saxon width (C)
  double mW = 0.;               // 3rd parameter
  int mExpon = 1;               // 1->3pF param., 2->3pG param
  double mRadDensMax = 0.;      // maximum radial density
  double mMaxRadius = 20.;      // maximum radius
  double mMinNucleonDist = 0.4; // minimum distance between nucleons
  std::array<double, kMaxNucleons> mXpos = {0};
  std::array<double, kMaxNucleons> mYpos = {0};
  std::array<double, kMaxNucleons> mZpos = {0};

  double estimateRadialMax(double rMin, double rMax, int samples = 1000);
  bool isFarEnough(int n, double x, double y, double z, double d2) const;
  double computeMaxRadius() const;

  ClassDef(NA6PNuclearShape, 1)
};

#endif
