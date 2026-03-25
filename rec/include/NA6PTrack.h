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

#ifndef NA6P_TRACK_H
#define NA6P_TRACK_H

#include <string>
#include <Rtypes.h>
#include "ExtTrackPar.h"

// Basic track class

class NA6PBaseCluster;

class NA6PTrack
{
 public:
  enum { kNDOF = 5,
         kMaxLr = 16 };
  enum { kY2 = 0,
         kZ2 = 2,
         kSnp2 = 5,
         kTgl2 = 9,
         kPtI2 = 14 };
  enum { kY,
         kZ,
         kSnp,
         kTgl,
         kPtI };

  NA6PTrack();
  NA6PTrack(const double* xyz, const double* pxyz, int sign, double errLoose = -1);
  NA6PTrack(const NA6PTrack&) = default;
  NA6PTrack& operator=(const NA6PTrack&) = default;
  virtual ~NA6PTrack() {}

  bool init(const double* xyz, const double* pxyz, int sign, double errLoose = -1);
  void imposeKinematics(const double* xyzLab, const double* cosinesLab, double en, double mass, int charge);
  void reset();
  void resetCovariance(float err = -1.);
  void setParam(const ExtTrackPar& p) { mExtTrack = p; }
  void setOuterParam(const ExtTrackPar& p) { mOuter = p; }
  void setVertexConstrainedParam(const ExtTrackPar& p) { mConstrained = p; }
  void setStatusRefitInward(bool status) { mStatusRefitInward = status; }
  void setStatusConstrained(bool status) { mStatusConstrained = status; }

  bool getStatusRefitInward() const { return mStatusRefitInward; }
  bool getStatusConstrained() const { return mStatusConstrained; }
  double getMass() const { return mMass; }
  double getChi2() const { return mChi2VT + mChi2MS; }
  double getChi2VT() const { return mChi2VT; }
  double getChi2MS() const { return mChi2MS; }
  double getMatchChi2() const { return mMatchChi2; }
  const ExtTrackPar& getTrackExtParam() const { return mExtTrack; }
  const double* getCovariance() const { return mExtTrack.getCovariance(); }
  double getCurvature(float by) const { return mExtTrack.getC(by); }
  std::array<float, 3> getCircleParams(float by) const;

  const ExtTrackPar& getOuterParam() const { return mOuter; }
  ExtTrackPar& getOuterParam() { return mOuter; }
  const ExtTrackPar& getVertexConstrainedParam() const { return mConstrained; }
  ExtTrackPar& getVertexConstrainedParam() { return mConstrained; }
  int getNVTHits() const { return mNClustersVT; }
  int getNMSHits() const { return mNClustersMS; }
  int getNTRHits() const { return mNClustersTR; }
  int getNHits() const { return mNClusters; }
  uint32_t getClusterMap() const { return mClusterMap; }
  int getClusterIndex(int lr) const { return lr < kMaxLr ? mClusterIndices[lr] : -1; }
  int getParticleLabel(int lr) const { return lr < kMaxLr ? mClusterPartID[lr] : -2; }
  int getParticleID() const { return mParticleID; }
  int getCAIteration() const { return mCAIteration; }

  double getAlpha() const { return mExtTrack.getAlpha(); }
  double getCharge() const { return mExtTrack.Charge(); }
  bool negDir() const { return std::abs(mExtTrack.getAlpha()) > M_PI / 2.; }
  double getXLab() const { return mExtTrack.getY(); }
  double getYLab() const { return mExtTrack.getZ(); }
  double getZLab() const { return negDir() ? -mExtTrack.getX() : mExtTrack.getX(); }
  double getX() const { return getXLab(); }
  double getY() const { return getYLab(); }
  double getZ() const { return getZLab(); }
  double getZLabOuter() const { return negDir() ? -mOuter.getX() : mOuter.getX(); }
  double getZLabVertexConstrained() const { return negDir() ? -mConstrained.getX() : mConstrained.getX(); }

  double getR() const
  {
    double x = getXLab(), y = getYLab(), r = x * x + y * y;
    return r > 0 ? std::sqrt(r) : 0;
  }
  double getXTF() const { return mExtTrack.getX(); }
  double getYTF() const { return mExtTrack.getY(); }
  double getZTF() const { return mExtTrack.getZ(); }
  void getXYZ(double* xyz) const;
  void getPXYZ(double* pxyz) const;
  template <typename T = float>
  std::array<T, 3> getXYZ() const
  {
    return std::array<T, 3>{static_cast<float>(getXLab()), static_cast<float>(getYLab()), static_cast<float>(getZLab())};
  }

  template <typename T = float>
  std::array<T, 3> getPXYZ() const
  {
    double pxyz[3];
    getPXYZ(pxyz);
    return {static_cast<float>(pxyz[0]), static_cast<float>(pxyz[1]), static_cast<float>(pxyz[2])};
  }
  double getP() const { return mExtTrack.getP(); }
  void getXYZOuter(double* xyz) const;
  void getPXYZOuter(double* pxyz) const;
  double getPOuter() const { return mOuter.getP(); }
  void getXYZVertexConstrained(double* xyz) const;
  void getPXYZVertexConstrained(double* pxyz) const;
  double getPVertexConstrained() const { return mConstrained.getP(); }

  double getSigmaX2() const { return mExtTrack.getSigmaY2(); }
  double getSigmaY2() const { return mExtTrack.getSigmaZ2(); }
  double getSigmaXY() const { return mExtTrack.getSigmaZY(); }
  double getSigmaYX() const { return getSigmaXY(); }
  double getSigmaP2() const;
  double getSigmaPX2() const;
  double getSigmaPY2() const;
  double getSigmaPZ2() const;
  double getChi2VTOuter() const { return mChi2VTOuter; }
  double getChi2MSOuter() const { return mChi2MSOuter; }
  double getChi2VTRefit() const { return mChi2VTRefit; }
  double getChi2MSRefit() const { return mChi2MSRefit; }
  double getNormChi2() const { return mNClusters < 3 ? 0 : (mChi2VT + mChi2MS) / ((mNClusters << 1) - kNDOF); }
  double getPredictedChi2(double* p, double* cov) const { return mExtTrack.getPredictedChi2(p, cov); }

  void setMass(double m) { mMass = m; }
  void setMatchChi2(double chi2) { mMatchChi2 = chi2; }
  void setChi2VT(double chi2) { mChi2VT = chi2; }
  void setChi2MS(double chi2) { mChi2MS = chi2; }
  void setChi2VTOuter(double chi2) { mChi2VTOuter = chi2; }
  void setChi2MSOuter(double chi2) { mChi2MSOuter = chi2; }
  void setChi2VTRefit(double chi2) { mChi2VTRefit = chi2; }
  void setChi2MSRefit(double chi2) { mChi2MSRefit = chi2; }
  void setParticleLabel(int idx, int lr)
  {
    if (lr < kMaxLr)
      mClusterPartID[lr] = idx;
  }
  void setClusterIndex(int idx, int lr)
  {
    if (lr < kMaxLr)
      mClusterIndices[lr] = idx;
  }
  void setParticleID(int idx) { mParticleID = idx; }
  void setCAIteration(int iter) { mCAIteration = iter; }

  static void lab2trk(const double* vLab, double* vTrk);
  static void trk2lab(const double* vTrk, double* vLab);

  template <typename ClusterType>
  void addCluster(const ClusterType* clu, int cluIndex, double chi2);
  bool correctForMeanMaterial(double xOverX0, double xTimesRho, double density = 0.f, double atomicZ = 0.f, double zOverA = 0.f, bool anglecorr = false)
  {
    return mExtTrack.correctForMeanMaterialGeneral(xOverX0, xTimesRho, mMass, anglecorr, density, atomicZ, zOverA);
  }
  bool propagateToZBxByBz(double z, double maxDZ = 1.0, double xOverX0 = 0., double xTimesRho = 0., bool outer = false, double rho = 0.f, double atomicZ = 0.f, double atomicZoverA = 0.f);
  bool propagateToDCA(NA6PTrack* partner);
  bool propagateToDCABeamAxis(float beamX, float beamY, float maxDCA);
  bool update(double p[2], double cov[3]) { return mExtTrack.Update(p, cov); }

  virtual void print() const;
  std::string asString() const;

 protected:
  double mMass = 0.140;                      // particle mass
  double mMatchChi2 = 0.f;                   // total chi2
  double mChi2VT = 0.f;                      // total chi2 VT
  double mChi2MS = 0.f;                      // total chi2 MS
  double mChi2VTOuter = 0.f;                 // total chi2 VT outward fit
  double mChi2MSOuter = 0.f;                 // total chi2 MS outward fit
  double mChi2VTRefit = 0.f;                 // total chi2 VT inward refit
  double mChi2MSRefit = 0.f;                 // total chi2 MS inward refit
  uint32_t mClusterMap = 0;                  // pattern of clusters per layer
  int mNClusters = 0;                        // total hits
  int mNClustersVT = 0;                      // total VT hits
  int mNClustersMS = 0;                      // total MS hits
  int mNClustersTR = 0;                      // total TR hits
  std::array<int, kMaxLr> mClusterIndices{}; // cluster indices
  std::array<int, kMaxLr> mClusterPartID{};  // particle ID (per cluster)
  int mParticleID = -1;                      // particle ID (MC truth)
  int mCAIteration = -1;                     //! CA iteration (for debug)
  ExtTrackPar mExtTrack;                     // track params
  ExtTrackPar mOuter{};                      // parametrization for outward fit
  ExtTrackPar mConstrained{};                // parametrization with vertex constrain
  bool mStatusRefitInward = false;           // boolean for status of inward refit
  bool mStatusConstrained = false;           // boolean for status of contrained track

 private:
  bool propagateToZBxByBz(double z, const double* bxyz) { return mExtTrack.propagateToBxByBz(z, bxyz); }
  bool propagateToZBxByBzOuter(double z, const double* bxyz) { return mOuter.propagateToBxByBz(z, bxyz); }

  ClassDefNV(NA6PTrack, 4)
};

//_______________________________________________________________________
inline void NA6PTrack::lab2trk(const double* vLab, double* vTrk)
{
  // convert alice coordinates to modified
  vTrk[0] = vLab[2];
  vTrk[1] = vLab[0];
  vTrk[2] = vLab[1];
}

//_______________________________________________________________________
inline void NA6PTrack::trk2lab(const double* vTrk, double* vLab)
{
  // convert modified coordinates to Lab ones
  vLab[0] = vTrk[1];
  vLab[1] = vTrk[2];
  vLab[2] = vTrk[0];
}

//_______________________________________________________________________
inline void NA6PTrack::getXYZ(double* xyz) const
{
  // track position in Lab coordinates
  double xyzTrk[3];
  mExtTrack.getXYZ(xyzTrk);
  trk2lab(xyzTrk, xyz);
}

//_______________________________________________________________________
inline void NA6PTrack::getPXYZ(double* pxyz) const
{
  // track position in Lab coordinates
  double pxyzTrk[3];
  mExtTrack.getPxPyPz(pxyzTrk);
  trk2lab(pxyzTrk, pxyz);
}

//_______________________________________________________________________
inline void NA6PTrack::getXYZOuter(double* xyz) const
{
  // track position in Lab coordinates
  double xyzTrk[3];
  mOuter.getXYZ(xyzTrk);
  trk2lab(xyzTrk, xyz);
}

//_______________________________________________________________________
inline void NA6PTrack::getPXYZOuter(double* pxyz) const
{
  // track position in Lab coordinates
  double pxyzTrk[3];
  mOuter.getPxPyPz(pxyzTrk);
  trk2lab(pxyzTrk, pxyz);
}
//_______________________________________________________________________
inline void NA6PTrack::getXYZVertexConstrained(double* xyz) const
{
  // track position in Lab coordinates
  double xyzTrk[3];
  mConstrained.getXYZ(xyzTrk);
  trk2lab(xyzTrk, xyz);
}

//_______________________________________________________________________
inline void NA6PTrack::getPXYZVertexConstrained(double* pxyz) const
{
  // track position in Lab coordinates
  double pxyzTrk[3];
  mConstrained.getPxPyPz(pxyzTrk);
  trk2lab(pxyzTrk, pxyz);
}

#endif
