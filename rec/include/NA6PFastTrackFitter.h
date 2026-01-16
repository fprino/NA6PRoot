// NA6PCCopyright

// Based on:
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_FAST_TRACK_FITTER_H
#define NA6P_FAST_TRACK_FITTER_H

#include <string>
#include <Rtypes.h>
#include "NA6PBaseCluster.h"

// Fast track fit based on Kalman filter

class TGeoManager;
class NA6PTrack;

class NA6PFastTrackFitter
{

 public:
  enum { kTwoPointSeed = 0,
         kThreePointSeed = 1 };

  NA6PFastTrackFitter();
  ~NA6PFastTrackFitter(){};

  // setters for configurable parameters
  void setMaxChi2Cl(double v = 10) { mMaxChi2Cl = v; }
  void setNLayers(int n)
  {
    mNLayers = n;
    mClusters.clear();
    mClusters.resize(n);
  }
  void setParticleHypothesis(int pdg);
  void enableMaterialCorrections() { mCorrectForMaterial = true; }
  void disableMaterialCorrections() { mCorrectForMaterial = false; }
  void setPropagateToPrimaryVertex(bool opt = true) { mPropagateToPrimVert = opt; }
  void setPrimaryVertexZ(double zvert)
  {
    mPrimVertZ = zvert;
    mIsPrimVertSet = true;
  }
  void setSeed(const double* pos, const double* mom, int charge = 1);
  void setCharge(int ch)
  {
    if (ch != 0)
      mCharge = ch;
  };
  void unsetSeed() { mIsSeedSet = false; }
  void setSeedFromTwoOutermostHits() { mSeedOption = kTwoPointSeed; }
  void setSeedFromThreeOutermostHits() { mSeedOption = kThreePointSeed; }
  void computeSeed(int dir = -1);
  void computeSeedOuter() { computeSeed(-1); }
  void computeSeedInner() { computeSeed(1); }
  void printSeed() const;

  void addCluster(int jLay, const NA6PBaseCluster& cl);
  void resetClusters()
  {
    std::fill(mClusters.begin(), mClusters.end(), nullptr);
  }
  void cleanupAndStartFit()
  {
    resetClusters();
    unsetSeed();
  }

  bool loadGeometry(const char* filename = "geometry.root", const char* geoname = "NA6P");

  NA6PTrack* fitTrackPoints(int dir = -1, NA6PTrack* seed = nullptr);
  NA6PTrack* fitTrackPointsInward() { return fitTrackPoints(-1); }
  NA6PTrack* fitTrackPointsOutward(NA6PTrack* seed = nullptr) { return fitTrackPoints(1, seed); }
  bool updateTrack(NA6PTrack* trc, const NA6PBaseCluster* cl) const;
  int propagateToZ(NA6PTrack* trc, double zFrom, double zTo, int dir) const;
  int propagateToZ(NA6PTrack* trc, double zTo) const;
  void getMeanMaterialBudgetFromGeom(double* start, double* end, double* mparam) const;

  static const Double_t kMassP;
  static const Double_t kMassK;
  static const Double_t kMassPi;
  static const Double_t kMassMu;
  static const Double_t kMassE;
  static const Double_t kAlmostZero;

 protected:
  int mNLayers = 5;                  // number of active
  double mMaxChi2Cl = 10.;           // max cluster-track chi2
  bool mIsSeedSet = false;           // flag for set seed
  int mSeedOption = kThreePointSeed; // seed option (see enum)
  double mSeedPos[3];                // seed for track position
  double mSeedMom[3];                // seed for track momentum
  int mCharge = 1;                   // track charge for seed
  double mMass = kMassPi;            // mass hypothesis for particle
  bool mPropagateToPrimVert = false; // flag for propagation to primary vertex
  double mPrimVertZ = 0.0;           // primary vertex z
  bool mIsPrimVertSet = false;       // flag for presence of prim vert z
  bool mCorrectForMaterial = true;   // flag for material corrections

  std::vector<const NA6PBaseCluster*> mClusters; // array with clusters

 protected:
  int getNumberOfClusters() const
  {
    int nClus = 0;
    for (int jLay = 0; jLay < mNLayers; ++jLay)
      if (mClusters[jLay])
        nClus++;
    return nClus;
  }

  ClassDefNV(NA6PFastTrackFitter, 1);
};

#endif
