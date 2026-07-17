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
#include "NA6PMCComposedLabel.h"
#include "NA6PTrackParCov.h"

class NA6PBaseCluster;

class NA6PTrack : public NA6PTrackParCov
{
 public:
  static constexpr int kMaxLr = 16;

  NA6PTrack();
  NA6PTrack(const float* xyz, const float* pxyz, int sign, float errLoose = -2);
  NA6PTrack& operator=(const NA6PTrack&) = default;
  NA6PTrack(const NA6PTrack&) = default;
  ~NA6PTrack() = default;

  void reset();
  void setInwardParam(const NA6PTrackParCov& p) { (*(NA6PTrackParCov*)this) = p; }
  void setOuterParam(const NA6PTrackParCov& p) { mOuter = p; }
  void setVertexConstrainedParam(const NA6PTrackParCov& p) { mConstrained = p; }
  bool getStatusRefitInward() const { return getInwardParam().isValid(); }
  bool getStatusConstrained() const { return getVertexConstrainedParam().isValid(); }
  bool getStatusRefitOutward() const { return getOuterParam().isValid(); }

  void setChi2(float c) { mChi2 = c; }
  float getChi2() const { return mChi2; }
  float getChi2Norm(int ndof = 5) const { return mNClusters > 2 ? mChi2 / (2 * mNClusters - ndof) : -1; }

  void setChi2Out(float c) { mChi2Outer = c; }
  float getChi2Outer() const { return mChi2Outer; }
  float getChi2OuterNorm(int ndof = 5) const { return mNClusters > 2 ? mChi2Outer / (2 * mNClusters - ndof) : -1; }

  const NA6PTrackParCov& getOuterParam() const { return mOuter; }
  NA6PTrackParCov& getOuterParam() { return mOuter; }

  const NA6PTrackParCov& getInwardParam() const { return *this; }
  NA6PTrackParCov& getInwardParam() { return *this; }

  const NA6PTrackParCov& getVertexConstrainedParam() const { return mConstrained; }
  NA6PTrackParCov& getVertexConstrainedParam() { return mConstrained; }

  int getNHits() const { return mNClusters; }

  uint32_t getClusterMap() const { return mClusterMap; }
  int getClusterIndex(int lr) const { return lr < kMaxLr ? mClusterIndices[lr] : -1; }
  int getParticleLabel(int lr) const { return lr < kMaxLr ? mClusterPartID[lr] : -2; }
  int getParticleID() const { return mParticleID; }
  NA6PMCComposedLabel getMCLabel() const { return mMCCompLab; }
  int getCAIteration() const { return mCAIteration; }

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
  void setMCLabel(NA6PMCComposedLabel& lab) { mMCCompLab = lab; }
  void setCAIteration(int iter) { mCAIteration = iter; }

  template <typename ClusterType>
  void addCluster(const ClusterType* clu);

  void print() const;
  std::string asString() const;

 protected:
  NA6PTrackParCov mOuter{};                  // parametrization for outward fit
  NA6PTrackParCov mConstrained{};            // parametrization with vertex constrain
  std::array<int, kMaxLr> mClusterIndices{}; // cluster indices
  std::array<int, kMaxLr> mClusterPartID{};  // particle ID (per cluster) // RSTOD why this is needed? This info must be available from the cluster indices

  float mChi2 = 0.f;                         // total chi2
  float mChi2Outer = 0.f;                    // total chi2 outward fit
  uint32_t mClusterMap = 0;                  // pattern of clusters per layer
  int mNClusters = 0;                        // total hits
  int mParticleID = -1;                      // particle ID (MC truth)
  NA6PMCComposedLabel mMCCompLab;            // MC composed label (will replace particle ID)
  int mCAIteration = -1;                     //! CA iteration (for debug)

 private:
  ClassDefNV(NA6PTrack, 3)
};

#endif
