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

#ifndef NA6P_TRACKER_CA_H
#define NA6P_TRACKER_CA_H

// #define _CHI2_TUNING_MODE_

#include <string>
#include <Rtypes.h>
#include "NA6PFastTrackFitter.h"
#include "NA6PTrack.h"
#include "NA6PTreeStreamRedirector.h"
#include "NA6PMCTruthContainer.h"

// Cellular Automaton track finder

class NA6PVertex;

// structures for temporary objects used in track finding
struct TrackletCandidate {
  int startingLayer;
  int firstClusterIndex;
  int secondClusterIndex;
  float tanL;
  float phi;
  float pxpz;
  float pypz;
};

struct CellCandidate {
  int startingLayer;
  int firstTrackletIndex;
  int secondTrackletIndex;
  std::array<int, 3> cluIDs;
  NA6PTrack trackFitFast;
  //    genfit::Track trackFit;
};

struct TrackCandidate {
  int innerLayer;
  int innerCellIndex;
  int outerLayer;
  int outerCellIndex;
  std::vector<int> cluIDs;
  int getCluIDsWOGaps(std::vector<int>& dest) const // fill vectors of cluster IDs w/o gaps. RSTOD: check if cluIDs has actually gaps
  {
    dest.clear();
    for (const auto cid : cluIDs) {
      if (cid >= 0) {
        dest.push_back(cid);
      }
    }
    return dest.size();
  }
};

struct TrackFitted {
  int innerLayer;
  int outerLayer;
  int nClus;
  std::vector<int> cluIDs;
  NA6PTrack trackFitFast;
  //    genfit::Track trackFit;
  float chi2ndf;
};

enum class ExtendDirection { kInward,
                             kOutward };

class NA6PTrackerCA
{
 public:
  enum RecoType {
    NONE = 0,
    VT,
    MS,
    NRecoType
  };
  static constexpr int kMaxIterationsCA = 10;

  NA6PTrackerCA();
  ~NA6PTrackerCA();

  // setters for configurable parameters
  void setNLayers(int n) { mNLayers = n; }
  void setStartLayer(int start) { mLayerStart = start; }
  void setMaxNumberOfSharedClusters(int n) { mMaxSharedClusters = n; }
  void setPropagateTracksToPrimaryVertex(bool opt = true) { mPropagateTracksToPrimaryVertex = opt; }
  void setDoOutwardPropagation(bool opt = true) { mDoOutwardPropagation = opt; }
  void setZForOutwardPropagation(float zout) { mZOutProp = zout; }
  void setDoInwardRefit(bool opt = true) { mDoInwardRefit = opt; }
  void setDoTrackConstrainedToPrimVert(bool opt = true) { mDoTrackConstrainedToPrimVert = opt; }
  void setNumberOfIterations(int nIter);
  void setIterationParams(int iter,
                          float maxDeltaThetaTracklets,
                          float maxDeltaPhiTracklets,
                          float maxDeltaTanLCells,
                          float maxDeltaPhiCells,
                          float maxDeltaPxPzCells,
                          float maxDeltaPyPzCells,
                          float maxChi2TrClCells,
                          float maxChi2ndfCells,
                          float maxChi2ndfTracks,
                          int minNClusTracks);
  void setPID(PID pid) { mPID = pid; }

  void setMaxPropagationStep(float step) { mTrackFitter->setMaxPropagationStep(step); }
  void configureFromRecoParamVT();
  void configureFromRecoParamMS();
  void setVerbosity(bool opt = true) { mVerbose = opt; }
  void setUseLinRef(bool v) { mUseLinRef = v; }

  void printConfiguration() const;
  int getNIterations() const { return mNIterationsCA; }
  bool loadGeometry(const char* filename, const char* geoname = "NA6P");

  template <typename ClusterType>
  void findTracks(std::vector<ClusterType>& cluArr, const NA6PMCTruthContainer& mcCluLabels, const NA6PVertex* primVert, bool useMCCont = true);
  template <typename ClusterType>
  void findTracks(std::vector<ClusterType>& cluArr, const NA6PVertex* primVert)
  {
    NA6PMCTruthContainer trCont; // dummy object
    findTracks(cluArr, trCont, primVert, false);
  }
  std::vector<NA6PTrack> getTracks();
  template <typename ClusterType>
  std::vector<std::pair<ClusterType, ClusterType>> findTracklets(int jFirstLay, int jLastLay, std::vector<ClusterType>& cluArr, const NA6PVertex* primVert);
  NA6PFastTrackFitter* getTrackFitter() { return mTrackFitter.get(); }

 protected:
  // methods used in tracking
  template <typename ClusterType>
  void sortClustersByLayerAndEta(std::vector<ClusterType>& cluArr,
                                 std::vector<int>& firstIndex,
                                 std::vector<int>& lastIndex);
  void sortTrackletsByLayerAndIndex(std::vector<TrackletCandidate>& tracklets,
                                    std::vector<int>& firstIndex,
                                    std::vector<int>& lastIndex);
  void sortCellsByLayerAndIndex(std::vector<CellCandidate>& cells,
                                std::vector<int>& firstIndex,
                                std::vector<int>& lastIndex);
  template <typename ClusterType>
  void computeLayerTracklets(const std::vector<ClusterType>& cluArr,
                             const std::vector<int>& firstIndex,
                             const std::vector<int>& lastIndex,
                             std::vector<TrackletCandidate>& tracklets,
                             float deltaThetaMax,
                             float deltaPhiMax);
  template <typename ClusterType>
  void computeLayerCells(const std::vector<TrackletCandidate>& tracklets,
                         const std::vector<int>& firstIndex,
                         const std::vector<int>& lastIndex,
                         const std::vector<ClusterType>& cluArr,
                         std::vector<CellCandidate>& cells,
                         float deltaTanLMax,
                         float deltaPhiMax,
                         float deltaPxPzMax,
                         float deltaPyPzMax,
                         float maxChi2TrClu,
                         float maxChi2NDF);
  template <typename ClusterType>
  float computeTrackToClusterChi2(const NA6PTrackParCov& track,
                                  const ClusterType& clu);
  template <typename ClusterType>
  float fitTrackPointsFast(const std::vector<int>& cluIDs,
                           const std::vector<ClusterType>& cluArr,
                           NA6PTrackParCov& fitTrack,
                           float maxChi2TrClu,
                           float maxChi2NDF);
  template <typename ClusterType>
  void findCellsNeighbours(const std::vector<CellCandidate>& cells,
                           const std::vector<int>& firstIndex,
                           const std::vector<int>& lastIndex,
                           std::vector<std::pair<int, int>>& cneigh,
                           const std::vector<ClusterType>& cluArr,
                           float maxChi2TrClu);
  template <typename ClusterType>
  std::vector<TrackCandidate> prolongSeed(const TrackCandidate& seed,
                                          const std::vector<CellCandidate>& cells,
                                          const std::vector<int>& firstIndex,
                                          const std::vector<int>& lastIndex,
                                          const std::vector<ClusterType>& cluArr,
                                          float maxChi2TrClu,
                                          ExtendDirection dir);
  template <typename ClusterType>
  void findRoads(const std::vector<std::pair<int, int>>& cneigh,
                 const std::vector<CellCandidate>& cells,
                 const std::vector<int>& firstIndex,
                 const std::vector<int>& lastIndex,
                 const std::vector<TrackletCandidate>& tracklets,
                 const std::vector<ClusterType>& cluArr,
                 std::vector<TrackCandidate>& trackCands,
                 float maxChi2TrClu);
  template <typename ClusterType>
  void fitAndSelectTracks(const std::vector<TrackCandidate>& trackCands,
                          const std::vector<ClusterType>& cluArr,
                          const NA6PMCTruthContainer& mcCluLabels,
                          std::vector<TrackFitted>& tracks,
                          const NA6PVertex* primVert,
                          float maxChi2TrClu,
                          int minNClu,
                          float maxChi2NDF,
                          bool useMCCont = false);
  template <typename T, typename ClusterType>
  void printStats(const std::vector<T>& candidates,
                  const std::vector<ClusterType>& cluArr,
                  const std::vector<CellCandidate>& cells,
                  const std::string& label,
                  int requiredClus = -1);

  template <typename T, typename ClusterType>
  std::pair<int, bool> getTrackMCTruthStatus(const T& clIDs, const std::vector<ClusterType>& cluArr) const;

 private:
  static constexpr std::array<std::string_view, RecoType::NRecoType> mRecoTypeNames{"NONE", "VerTel", "MS"};
  RecoType mRecoType = RecoType::NONE;
  PID mPID{};
  int mNLayers = 5;
  int mLayerStart = 0;
  float mPrimVertPos[3] = {};
  std::unique_ptr<NA6PFastTrackFitter> mTrackFitter;
  std::vector<bool> mIsClusterUsed = {};
  std::vector<TrackFitted> mFinalTracks = {};
  int mMaxSharedClusters = 0;
  bool mUseLinRef = true;
  bool mVerbose = false;
  bool mPropagateTracksToPrimaryVertex = false;
  bool mDoOutwardPropagation = false;
  bool mDoInwardRefit = false;
  bool mDoTrackConstrainedToPrimVert = false;
  float mZOutProp = 38.1175;
  int mNIterationsCA = 2;
  int mCurIteration = -1;
  int mNDOF = 5;
  float mMaxDeltaThetaTrackletsCA[kMaxIterationsCA] = {0.04, 0.1, 0.15, 0.3, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float mMaxDeltaPhiTrackletsCA[kMaxIterationsCA] = {0.1, 0.2, 0.25, 0.5, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float mMaxDeltaTanLCellsCA[kMaxIterationsCA] = {4., 9., 18., 40., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float mMaxDeltaPhiCellsCA[kMaxIterationsCA] = {0.4, 0.6, 1.2, 2.4, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float mMaxDeltaPxPzCellsCA[kMaxIterationsCA] = {0.02, 0.05, 0.1, 0.2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float mMaxDeltaPyPzCellsCA[kMaxIterationsCA] = {2e-3, 7.5e-3, 1.5e-2, 3.e-2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float mMaxChi2TrClCellsCA[kMaxIterationsCA] = {100., 500., 999., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float mMaxChi2ndfCellsCA[kMaxIterationsCA] = {100., 500., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float mMaxChi2ndfTracksCA[kMaxIterationsCA] = {100., 500., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  int mMinNClusTracksCA[kMaxIterationsCA] = {5, 3, 0, 0, 0, 0, 0, 0, 0, 0};

  std::unique_ptr<NA6PTreeStreamRedirector> dbgStream;

  ClassDefNV(NA6PTrackerCA, 1);
};

template <typename T, typename ClusterType>
std::pair<int, bool> NA6PTrackerCA::getTrackMCTruthStatus(const T& clIDs, const std::vector<ClusterType>& cluArr) const
{
  std::vector<std::pair<int, int>> occurrences;
  int ncl = 0;
  for (auto id : clIDs) {
    if (id < 0) {
      continue;
    }
    ncl++;
    const int pid = cluArr[id].getParticleID();
    bool found{false};
    for (int i = 0; i < (int)occurrences.size(); i++) {
      if (pid == occurrences[i].first) {
        occurrences[i].second++;
        found = true;
        break;
      }
    }
    if (!found) {
      occurrences.emplace_back(pid, 1);
    }
  }
  std::sort(std::begin(occurrences), std::end(occurrences), [](auto e1, auto e2) { return e1.second > e2.second; });
  return {occurrences.front().second == ncl, occurrences.front().first};
}

#endif
