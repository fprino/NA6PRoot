// NA6PCCopyright

#include <fmt/format.h>
#include <iostream>
#include <fairlogger/Logger.h>
#include <TGeoManager.h>
#include <TSystem.h>
#include "Propagator.h"
#include "MagneticField.h"
#include "NA6PRecoParam.h"
#include "NA6PVertex.h"
#include "NA6PTrackerCA.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PVerTelCluster.h"

static constexpr float kThetaMax = M_PI_2 * 0.99f;

namespace
{
std::string formatClusterIDs(const std::vector<int>& ids)
{
  std::string out;
  for (size_t i = 0; i < ids.size(); ++i) {
    if (i) {
      out += ",";
    }
    out += std::to_string(ids[i]);
  }
  return out;
}
} // namespace

//______________________________________________________________________
NA6PTrackerCA::NA6PTrackerCA()
{
  // RSTODO create init method for the fitter
  mTrackFitter = std::make_unique<NA6PFastTrackFitter>();
  mTrackFitter->setPropagateToPrimaryVertex(false);
}

//______________________________________________________________________
NA6PTrackerCA::~NA6PTrackerCA()
{
#ifdef _CHI2_TUNING_MODE_
  if (dbgStream) {
    dbgStream->Close();
  }
#endif
}

//______________________________________________________________________
void NA6PTrackerCA::setNumberOfIterations(int nIter)
{
  if (nIter >= 0 && nIter <= kMaxIterationsCA)
    mNIterationsCA = nIter;
  else
    LOGP(error, "Number of iterations should be <= {}", kMaxIterationsCA);
}

void NA6PTrackerCA::setIterationParams(int iter,
                                       float maxDeltaThetaTracklets,
                                       float maxDeltaPhiTracklets,
                                       float maxDeltaTanLCells,
                                       float maxDeltaPhiCells,
                                       float maxDeltaPxPzCells,
                                       float maxDeltaPyPzCells,
                                       float maxChi2TrClCells,
                                       float maxChi2ndfCells,
                                       float maxChi2ndfTracks,
                                       int minNClusTracks)
{
  if (iter < 0 || iter >= mNIterationsCA) {
    LOGP(error, "Iteration {} out of allowed range [0, {}]", iter, mNIterationsCA - 1);
    return;
  }
  mMaxDeltaThetaTrackletsCA[iter] = maxDeltaThetaTracklets;
  mMaxDeltaPhiTrackletsCA[iter] = maxDeltaPhiTracklets;
  mMaxDeltaTanLCellsCA[iter] = maxDeltaTanLCells;
  mMaxDeltaPhiCellsCA[iter] = maxDeltaPhiCells;
  mMaxDeltaPxPzCellsCA[iter] = maxDeltaPxPzCells;
  mMaxDeltaPyPzCellsCA[iter] = maxDeltaPyPzCells;
  mMaxChi2TrClCellsCA[iter] = maxChi2TrClCells;
  mMaxChi2ndfCellsCA[iter] = maxChi2ndfCells;
  mMaxChi2ndfTracksCA[iter] = maxChi2ndfTracks;
  mMinNClusTracksCA[iter] = minNClusTracks;
}

void NA6PTrackerCA::configureFromRecoParamVT()
{
  mRecoType = RecoType::VT;
  setPID(PID::Pion);
  const auto& param = NA6PRecoParam::Instance();
  setDoOutwardPropagation(param.vtDoOutwardPropagation);
  setZForOutwardPropagation(param.vtZOutProp);
  setDoInwardRefit(param.vtDoInwardRefit);
  setPropagateTracksToPrimaryVertex(param.vtPropagateTracksToPV);
  setDoTrackConstrainedToPrimVert(param.vtDoConstrainedTrack);
  setNumberOfIterations(param.vtNIterationsTrackerCA);
  setNLayers(param.vtNLayers);
  setMaxPropagationStep(param.maxPropagationStep);
  setUseLinRef(param.useLinRefVT);
  mTrackFitter->setSeedImprovePrec(param.seedImprovePrecVT);
  mTrackFitter->setPID(mPID);

  for (int jIter = 0; jIter < mNIterationsCA; ++jIter) {
    setIterationParams(jIter,
                       param.vtMaxDeltaThetaTrackletsCA[jIter],
                       param.vtMaxDeltaPhiTrackletsCA[jIter],
                       param.vtMaxDeltaTanLCellsCA[jIter],
                       param.vtMaxDeltaPhiCellsCA[jIter],
                       param.vtMaxDeltaPxPzCellsCA[jIter],
                       param.vtMaxDeltaPyPzCellsCA[jIter],
                       param.vtMaxChi2TrClCellsCA[jIter],
                       param.vtMaxChi2ndfCellsCA[jIter],
                       param.vtMaxChi2ndfTracksCA[jIter],
                       param.vtMinNClusTracksCA[jIter]);
  }
#ifdef _CHI2_TUNING_MODE_
  dbgStream = std::make_unique<NA6PTreeStreamRedirector>(fmt::format("trackerCA_{}_dbg.root", mRecoTypeNames[mRecoType]).c_str());
#endif
}

void NA6PTrackerCA::configureFromRecoParamMS()
{
  mRecoType = RecoType::MS;
  setPID(PID::Muon);
  const auto& param = NA6PRecoParam::Instance();
  setNumberOfIterations(param.msNIterationsTrackerCA);
  setNLayers(param.msNLayers);
  setStartLayer(param.vtNLayers);
  setDoTrackConstrainedToPrimVert(param.msDoConstrainedTrack);
  setMaxPropagationStep(param.maxPropagationStep);
  setUseLinRef(param.useLinRefMS);
  mTrackFitter->setSeedImprovePrec(param.seedImprovePrecMS);
  mTrackFitter->setPID(mPID);

  for (int jIter = 0; jIter < mNIterationsCA; ++jIter) {
    setIterationParams(jIter,
                       param.msMaxDeltaThetaTrackletsCA[jIter],
                       param.msMaxDeltaPhiTrackletsCA[jIter],
                       param.msMaxDeltaTanLCellsCA[jIter],
                       param.msMaxDeltaPhiCellsCA[jIter],
                       param.msMaxDeltaPxPzCellsCA[jIter],
                       param.msMaxDeltaPyPzCellsCA[jIter],
                       param.msMaxChi2TrClCellsCA[jIter],
                       param.msMaxChi2ndfCellsCA[jIter],
                       param.msMaxChi2ndfTracksCA[jIter],
                       param.msMinNClusTracksCA[jIter]);
  }
#ifdef _CHI2_TUNING_MODE_
  dbgStream = std::make_unique<NA6PTreeStreamRedirector>(fmt::format("trackerCA_{}_dbg.root", mRecoTypeNames[mRecoType]).c_str());
#endif
}

void NA6PTrackerCA::printConfiguration() const
{
  std::cout << "=== Tracker CA Configuration ===\n";
  std::cout << "Number of iterations: " << mNIterationsCA << "\n\n";

  std::cout << std::setw(5) << "Iter"
            << std::setw(10) << "dThetaTrk"
            << std::setw(10) << "dPhiTrk"
            << std::setw(10) << "dTanLCell"
            << std::setw(10) << "dPhiCell"
            << std::setw(10) << "dPxPz"
            << std::setw(10) << "dPyPz"
            << std::setw(10) << "Chi2Cell"
            << std::setw(10) << "Chi2CellNdf"
            << std::setw(10) << "Chi2TrackNdf"
            << std::setw(8) << "MinNClu"
            << "\n";

  for (int i = 0; i < mNIterationsCA; ++i) {
    std::cout << std::setw(5) << i
              << std::setw(10) << mMaxDeltaThetaTrackletsCA[i]
              << std::setw(10) << mMaxDeltaPhiTrackletsCA[i]
              << std::setw(10) << mMaxDeltaTanLCellsCA[i]
              << std::setw(10) << mMaxDeltaPhiCellsCA[i]
              << std::setw(10) << mMaxDeltaPxPzCellsCA[i]
              << std::setw(10) << mMaxDeltaPyPzCellsCA[i]
              << std::setw(10) << mMaxChi2TrClCellsCA[i]
              << std::setw(10) << mMaxChi2ndfCellsCA[i]
              << std::setw(10) << mMaxChi2ndfTracksCA[i]
              << std::setw(8) << mMinNClusTracksCA[i]
              << "\n";
  }
  std::cout << "==============================\n";
}

//______________________________________________________________________

bool NA6PTrackerCA::loadGeometry(const char* filename, const char* geoname)
{
  return mTrackFitter->loadGeometry(filename, geoname);
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::sortClustersByLayerAndEta(std::vector<ClusterType>& cluArr,
                                              std::vector<int>& firstIndex,
                                              std::vector<int>& lastIndex)
{

  // Count clusters per layer for all layers
  int maxLayers = 25;
  std::vector<int> count(maxLayers, 0);
  for (const auto& clu : cluArr) {
    int jLay = clu.getLayer();
    count[jLay]++;
  }
  // Compute starting offset for each layer
  std::vector<int> firstAll(maxLayers);
  firstAll[0] = 0;
  for (int i = 1; i < maxLayers; ++i) {
    firstAll[i] = firstAll[i - 1] + count[i - 1];
  }

  // Fill firstIndex and lastIndex only for selected layers
  firstIndex.resize(mNLayers);
  lastIndex.resize(mNLayers);
  for (int j = 0; j < mNLayers; ++j) {
    int jLay = mLayerStart + j;
    firstIndex[j] = firstAll[jLay];
    lastIndex[j] = firstAll[jLay] + count[jLay];
  }

  // Create a reordered vector based on layer grouping
  std::vector<int> countReord = firstAll;
  std::vector<ClusterType> reordered(cluArr.size());
  for (const auto& clu : cluArr) {
    int jLay = clu.getLayer();
    if (jLay >= 0 && jLay < maxLayers)
      reordered[countReord[jLay]++] = clu;
  }
  cluArr = std::move(reordered);
  // sort by theta within each layer (use z^2/r^2 as a proxy of theta to avoid sqrt and atan)
  const float pvx = mPrimVertPos[0];
  const float pvy = mPrimVertPos[1];
  const float pvz = mPrimVertPos[2];
  for (int jLay = 0; jLay < maxLayers; jLay++) {
    if (count[jLay] == 0)
      continue;
    auto first = cluArr.begin() + firstAll[jLay];
    auto last = cluArr.begin() + firstAll[jLay] + count[jLay];
    std::sort(first, last, [pvx, pvy, pvz](const ClusterType& a, const ClusterType& b) {
      float xa = a.getX() - pvx;
      float ya = a.getY() - pvy;
      float za = a.getZ() - pvz;
      float xb = b.getX() - pvx;
      float yb = b.getY() - pvy;
      float zb = b.getZ() - pvz;
      float r2a = xa * xa + ya * ya;
      float r2b = xb * xb + yb * yb;
      return za * za * r2b < zb * zb * r2a;
    });
  }
}

//______________________________________________________________________

void NA6PTrackerCA::sortTrackletsByLayerAndIndex(std::vector<TrackletCandidate>& tracklets,
                                                 std::vector<int>& firstIndex,
                                                 std::vector<int>& lastIndex)
{

  // count clusters per layer
  std::vector<int> count(mNLayers - 1, 0);
  for (const auto& trkl : tracklets) {
    int jLay = trkl.startingLayer;
    if (jLay >= 0 && jLay < mNLayers - 1) {
      count[jLay]++;
    }
  }
  // starting offset for each layer
  firstIndex.resize(mNLayers - 1);
  lastIndex.resize(mNLayers - 1);
  firstIndex[0] = 0;
  lastIndex[0] = count[0];
  for (int jLay = 1; jLay < mNLayers - 1; jLay++) {
    firstIndex[jLay] = firstIndex[jLay - 1] + count[jLay - 1];
    lastIndex[jLay] = firstIndex[jLay] + count[jLay];
  }
  // Create a reordered vector based on layer grouping
  std::vector<int> countReord = firstIndex;
  std::vector<TrackletCandidate> reordered(tracklets.size());
  for (const auto& trkl : tracklets) {
    int jLay = trkl.startingLayer;
    if (jLay >= 0 && jLay < mNLayers - 1)
      reordered[countReord[jLay]++] = trkl;
  }
  tracklets = std::move(reordered);
  // sort by cluster index within each layer
  for (int jLay = 0; jLay < mNLayers - 1; jLay++) {
    auto first = tracklets.begin() + firstIndex[jLay];
    auto last = tracklets.begin() + lastIndex[jLay];
    std::sort(first, last, [](const TrackletCandidate& a, const TrackletCandidate& b) {
      return a.firstClusterIndex < b.firstClusterIndex || (a.firstClusterIndex == b.firstClusterIndex && a.secondClusterIndex < b.secondClusterIndex);
    });
  }
}

//______________________________________________________________________

void NA6PTrackerCA::sortCellsByLayerAndIndex(std::vector<CellCandidate>& cells,
                                             std::vector<int>& firstIndex,
                                             std::vector<int>& lastIndex)
{

  // count cells per layer
  std::vector<int> count(mNLayers - 2, 0);
  for (const auto& cell : cells) {
    int jLay = cell.startingLayer;
    if (jLay >= 0 && jLay < mNLayers - 2) {
      count[jLay]++;
    }
  }
  // starting offset for each layer
  firstIndex.resize(mNLayers - 2);
  lastIndex.resize(mNLayers - 2);
  firstIndex[0] = 0;
  lastIndex[0] = count[0];
  for (int jLay = 1; jLay < mNLayers - 2; jLay++) {
    firstIndex[jLay] = firstIndex[jLay - 1] + count[jLay - 1];
    lastIndex[jLay] = firstIndex[jLay] + count[jLay];
  }
  // Create a reordered vector based on layer grouping
  std::vector<int> countReord = firstIndex;
  std::vector<CellCandidate> reordered(cells.size());
  for (const auto& cell : cells) {
    int jLay = cell.startingLayer;
    if (jLay >= 0 && jLay < mNLayers - 2)
      reordered[countReord[jLay]++] = cell;
  }
  cells = std::move(reordered);
  // sort by tracklet index within each layer
  for (int jLay = 0; jLay < mNLayers - 2; jLay++) {
    auto first = cells.begin() + firstIndex[jLay];
    auto last = cells.begin() + lastIndex[jLay];
    std::sort(first, last, [](const CellCandidate& a, const CellCandidate& b) {
      return a.firstTrackletIndex < b.firstTrackletIndex || (a.firstTrackletIndex == b.firstTrackletIndex && a.secondTrackletIndex < b.secondTrackletIndex);
    });
  }
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::computeLayerTracklets(const std::vector<ClusterType>& cluArr,
                                          const std::vector<int>& firstIndex,
                                          const std::vector<int>& lastIndex,
                                          std::vector<TrackletCandidate>& tracklets,
                                          float deltaThetaMax,
                                          float deltaPhiMax)
{

  tracklets.clear();
  const float pvx = mPrimVertPos[0];
  const float pvy = mPrimVertPos[1];
  const float pvz = mPrimVertPos[2];

  for (int iLayer = 0; iLayer < mNLayers - 1; ++iLayer) {
    auto layerBegin = cluArr.begin() + firstIndex[iLayer + 1];
    auto layerEnd = cluArr.begin() + lastIndex[iLayer + 1];
    for (int jClu1 = firstIndex[iLayer]; jClu1 < lastIndex[iLayer]; ++jClu1) {
      if (mIsClusterUsed[jClu1])
        continue;
      const ClusterType& clu1 = cluArr[jClu1];
      float x1 = clu1.getX() - pvx;
      float y1 = clu1.getY() - pvy;
      float z1 = clu1.getZ() - pvz;
      float r1 = std::sqrt(x1 * x1 + y1 * y1);
      float theta1 = std::atan2(z1, r1);
      float phi1 = std::atan2(y1, x1);
      float tanth2Min = std::max(0.f, std::tan(theta1 - 1.2f * deltaThetaMax)); // 1.2 is a safety margin, the tan(theta) range is limited at zero to protect for the case in which theta1 - 1.2 * deltaThetaMax is <0, which would give a negative tangent. The minimum acceptable value for theta is 0
      float tanth2Max = std::tan(std::min(kThetaMax, theta1 + 1.2f * deltaThetaMax));
      auto lower = std::partition_point(layerBegin, layerEnd,
                                        [pvx, pvy, pvz, tanth2Min](const ClusterType& clu) {
                                          float x = clu.getX() - pvx;
                                          float y = clu.getY() - pvy;
                                          float z = clu.getZ() - pvz;
                                          float r2 = x * x + y * y;
                                          float tan2 = z * z / r2;
                                          return tan2 < tanth2Min * tanth2Min;
                                        });

      auto upper = std::partition_point(layerBegin, layerEnd,
                                        [pvx, pvy, pvz, tanth2Max](const ClusterType& clu) {
                                          float x = clu.getX() - pvx;
                                          float y = clu.getY() - pvy;
                                          float z = clu.getZ() - pvz;
                                          float r2 = x * x + y * y;
                                          float tan2 = z * z / r2;
                                          return tan2 <= tanth2Max * tanth2Max;
                                        });
      int lowerIdx = std::distance(cluArr.begin(), lower);
      int upperIdx = std::distance(cluArr.begin(), upper);
      for (int jClu2 = lowerIdx; jClu2 < upperIdx; ++jClu2) {
        if (mIsClusterUsed[jClu2])
          continue;
        const ClusterType& clu2 = cluArr[jClu2];
        float x2 = clu2.getX() - pvx;
        float y2 = clu2.getY() - pvy;
        float z2 = clu2.getZ() - pvz;
        float r2 = std::sqrt(x2 * x2 + y2 * y2);
        float theta2 = std::atan2(z2, r2);
        float phi2 = std::atan2(y2, x2);
        float dphi = phi2 - phi1;
        if (dphi > M_PI)
          dphi -= 2 * M_PI;
        else if (dphi < -M_PI)
          dphi += 2 * M_PI;
        if (std::abs(theta2 - theta1) < deltaThetaMax && std::abs(dphi) < deltaPhiMax) {
          float phi = std::atan2(y2 - y1, x2 - x1);
          float tanL = (z2 - z1) / (r2 - r1);
          float pxpz = (x2 - x1) / (z2 - z1);
          float pypz = (y2 - y1) / (z2 - z1);
          tracklets.emplace_back(iLayer, jClu1, jClu2, tanL, phi, pxpz, pypz);
          // do not assign the clusters as used, it will be done when the tracklets are used into tracks
          // mIsClusterUsed[jClu1] = true;
          // mIsClusterUsed[jClu2] = true;
        }
      }
    }
  }
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::computeLayerCells(const std::vector<TrackletCandidate>& tracklets,
                                      const std::vector<int>& firstIndex,
                                      const std::vector<int>& lastIndex,
                                      const std::vector<ClusterType>& cluArr,
                                      std::vector<CellCandidate>& cells,
                                      float deltaTanLMax,
                                      float deltaPhiMax,
                                      float deltaPxPzMax,
                                      float deltaPyPzMax,
                                      float maxChi2TrClu,
                                      float maxChi2NDF)
{

  cells.clear();
  for (int iLayer = 0; iLayer < mNLayers - 2; ++iLayer) {
    auto layerBegin = tracklets.begin() + firstIndex[iLayer + 1];
    auto layerEnd = tracklets.begin() + lastIndex[iLayer + 1];
    for (int jTrkl1 = firstIndex[iLayer]; jTrkl1 < lastIndex[iLayer]; ++jTrkl1) {
      const TrackletCandidate& trkl1 = tracklets[jTrkl1];
      const int nextLayerClusterIndex = trkl1.secondClusterIndex;
      auto lower = std::partition_point(layerBegin, layerEnd,
                                        [nextLayerClusterIndex](const TrackletCandidate& trkl) {
                                          return trkl.firstClusterIndex < nextLayerClusterIndex;
                                        });
      auto upper = std::partition_point(layerBegin, layerEnd,
                                        [nextLayerClusterIndex](const TrackletCandidate& trkl) {
                                          return trkl.firstClusterIndex <= nextLayerClusterIndex;
                                        });
      for (auto it = lower; it != upper; ++it) {
        int jTrkl2 = it - tracklets.begin();
        const TrackletCandidate& trkl2 = *it;
        if (trkl2.firstClusterIndex != nextLayerClusterIndex)
          continue;
        const float deltaTanLambda = std::abs(trkl2.tanL - trkl1.tanL);
        float dphi = trkl2.phi - trkl1.phi;
        if (dphi > M_PI)
          dphi -= 2 * M_PI;
        else if (dphi < -M_PI)
          dphi += 2 * M_PI;
        float deltapxpz = std::abs(trkl2.pxpz - trkl1.pxpz);
        float deltapypz = std::abs(trkl2.pypz - trkl1.pypz);
        if (deltapypz < deltaPyPzMax && deltapxpz < deltaPxPzMax && deltaTanLambda < deltaTanLMax && std::abs(dphi) < deltaPhiMax) {
          std::array<int, 3> cluIDs = {trkl1.firstClusterIndex, trkl2.firstClusterIndex, trkl2.secondClusterIndex};
          NA6PTrack fitTrackFast;
          //   genfit::Track fitTrack;
          float chi = fitTrackPointsFast(std::vector<int>(cluIDs.begin(), cluIDs.end()), cluArr, fitTrackFast, maxChi2TrClu, maxChi2NDF);
          if (chi >= 0) {
#ifdef _CHI2_TUNING_MODE_
            auto mcTruth = mTrackFitter->getMCTruthStatus();
            if (mcTruth.second) {
              (*dbgStream) << "chiCells"
                           << "iter=" << mCurIteration << "trc=" << ((NA6PTrackParCov&)fitTrackFast) << "chi2vec=" << mTrackFitter->getChi2Buffer() << "chi2Tot=" << chi << "\n";
            }
#endif
            if (mVerbose) {
              LOGP(info, "DBG cell accept lay={} trkl=({},{}) clu=({},{},{}) chi={}",
                   iLayer, jTrkl1, jTrkl2, cluIDs[0], cluIDs[1], cluIDs[2], chi);
            }
            cells.emplace_back(iLayer, jTrkl1, jTrkl2, std::move(cluIDs), fitTrackFast);
          } else if (mVerbose) {
            LOGP(info, "DBG cell reject lay={} trkl=({},{}) clu=({},{},{})",
                 iLayer, jTrkl1, jTrkl2, cluIDs[0], cluIDs[1], cluIDs[2]);
          }
        }
      }
    }
  }
}

//______________________________________________________________________
template <typename ClusterType>
float NA6PTrackerCA::computeTrackToClusterChi2(const NA6PTrackParCov& track, const ClusterType& clu)
{
  NA6PTrackParCov t{track};
  if (!Propagator::Instance()->propagateToZ(t, clu.getZ())) {
    if (mVerbose) {
      LOGP(info, "DBG trk-cl chi2 fail: propagate to z={} state={}", clu.getZ(), track.asString());
    }
    return 1e99;
  }
  return t.getPredictedChi2(clu);
}

//______________________________________________________________________

template <typename ClusterType>
float NA6PTrackerCA::fitTrackPointsFast(const std::vector<int>& cluIDs,
                                        const std::vector<ClusterType>& cluArr,
                                        NA6PTrackParCov& fitTrack,
                                        float maxChi2TrClu,
                                        float maxChi2NDF)
{

  int nClus = cluIDs.size();
  mTrackFitter->cleanupAndStartFit();
  mTrackFitter->setMaxChi2Cl(maxChi2TrClu);
  for (int jClu = 0; jClu < nClus; jClu++) {
    mTrackFitter->addCluster(cluArr[cluIDs[jClu]]);
  }
  if (!mTrackFitter->computeSeed(-1, &fitTrack)) {
    if (mVerbose) {
      LOGP(info, "DBG fit fail: seed nClus={} cluIDs={}", nClus, formatClusterIDs(cluIDs));
    }
    return -1.f;
  }
  auto chiTot = mTrackFitter->fitSeedInward(fitTrack, true, mUseLinRef);
  if (chiTot < 0.f) {
    if (mVerbose) {
      LOGP(info, "DBG fit fail: fitSeedInward nClus={} cluIDs={} seed={}",
           nClus, formatClusterIDs(cluIDs), fitTrack.asString());
    }
    return -1.f;
  }
  float chi2ndf = nClus < 3 ? 0.f : chiTot / (2 * nClus - mNDOF);
  if (chi2ndf > maxChi2NDF) {
    if (mVerbose) {
      LOGP(info, "DBG fit fail: chi2ndf={} max={} nClus={} cluIDs={}",
           chi2ndf, maxChi2NDF, nClus, formatClusterIDs(cluIDs));
    }
    return -1.f;
  }
  return chiTot;
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::findCellsNeighbours(const std::vector<CellCandidate>& cells,
                                        const std::vector<int>& firstIndex,
                                        const std::vector<int>& lastIndex,
                                        std::vector<std::pair<int, int>>& cneigh,
                                        const std::vector<ClusterType>& cluArr,
                                        float maxChi2TrClu)
{

  cneigh.clear();

  for (int iLayer = 0; iLayer < mNLayers - 3; ++iLayer) {
    auto layerBegin = cells.begin() + firstIndex[iLayer + 1];
    auto layerEnd = cells.begin() + lastIndex[iLayer + 1];
    for (int jCe1 = firstIndex[iLayer]; jCe1 < lastIndex[iLayer]; ++jCe1) {
      const CellCandidate& cell1 = cells[jCe1];
      const int nextLayerTrackletIndex = cell1.secondTrackletIndex;
      auto lower = std::partition_point(layerBegin, layerEnd,
                                        [&](const CellCandidate& c) { return c.firstTrackletIndex < nextLayerTrackletIndex; });
      auto upper = std::partition_point(layerBegin, layerEnd,
                                        [&](const CellCandidate& c) { return c.firstTrackletIndex <= nextLayerTrackletIndex; });
      for (auto it = lower; it != upper; ++it) {
        int jCe2 = it - cells.begin();
        const CellCandidate& cell2 = *it;
        if (cell2.firstTrackletIndex != nextLayerTrackletIndex)
          continue;
        if (cell1.cluIDs[1] != cell2.cluIDs[0] || cell1.cluIDs[2] != cell2.cluIDs[1]) {
          LOGP(error, "mismatch in cluIDs");
          continue;
        }

        // Check compatibility of cells via track-cluster chi2
        // Check done only for the inward propagation of the outer cell
        //   because we have the track parameterization at the innermost point of each cell
        // int jClu2 = cell2.cluIDs[2];
        // const auto& clu2 = cluArr[jClu2];
        // float cluchi2a = computeTrackToClusterChi2(cell1.trackFitFast, clu2);
        // if (cluchi2a > maxChi2TrClu)
        //   continue;

        int jClu0 = cell1.cluIDs[0];
        const auto& clu0 = cluArr[jClu0];
        float cluchi2b = computeTrackToClusterChi2(cell2.trackFitFast, clu0);
#ifdef _CHI2_TUNING_MODE_
        auto mcTruth = getTrackMCTruthStatus(cell2.cluIDs, cluArr);
        if (mcTruth.second && mcTruth.first == clu0.getParticleID()) {
          (*dbgStream) << "chiNb"
                       << "iter=" << mCurIteration << "trc=" << ((NA6PTrackParCov&)cell2.trackFitFast) << "trackNcl=" << cell2.trackFitFast.getNHits() << "trackChi2=" << cell2.trackFitFast.getChi2() << "chi2cl=" << cluchi2b << "\n";
        }
#endif
        if (cluchi2b > maxChi2TrClu) {
          if (mVerbose) {
            LOGP(info, "DBG neigh reject: lay={} c1={} c2={} addClu={} chi2={} max={}",
                 iLayer, jCe1, jCe2, jClu0, cluchi2b, maxChi2TrClu);
          }
          continue;
        }

        cneigh.push_back(std::make_pair(jCe1, jCe2));
      }
    }
  }
}

//______________________________________________________________________

template <typename ClusterType>
std::vector<TrackCandidate> NA6PTrackerCA::prolongSeed(const TrackCandidate& seed,
                                                       const std::vector<CellCandidate>& cells,
                                                       const std::vector<int>& firstIndex,
                                                       const std::vector<int>& lastIndex,
                                                       const std::vector<ClusterType>& cluArr,
                                                       float maxChi2TrClu,
                                                       ExtendDirection dir)
{

  std::vector<TrackCandidate> current = {seed};
  std::vector<TrackCandidate> next;
  std::vector<TrackCandidate> result;

  int step = (dir == ExtendDirection::kInward) ? -1 : 1;
  int startLayer = (dir == ExtendDirection::kInward) ? seed.innerLayer - 1 : seed.outerLayer + 1;
  int endLayer = (dir == ExtendDirection::kInward) ? -1 : mNLayers;
  int maxValidLayerToSearch = firstIndex.size();

  for (int iLayer = startLayer; iLayer != endLayer; iLayer += step) {
    next.clear();
    bool foundAnyProlongation = false;
    for (auto& cand : current) {
      const CellCandidate& refCell = (dir == ExtendDirection::kInward) ? cells[cand.innerCellIndex] : cells[cand.outerCellIndex];
      const int nextLayerTrackletIndex = (dir == ExtendDirection::kInward) ? refCell.firstTrackletIndex : refCell.secondTrackletIndex;
      // in case of outward prolongation, pick up firstIndex and lastIndex at iLayer-2 (the cell contains 3 clus in layers: iLayer-2, iLayer-1 and iLayer)
      int layToSearch = (dir == ExtendDirection::kInward) ? iLayer : iLayer - 2;
      if (layToSearch < 0 || layToSearch >= maxValidLayerToSearch) {
        LOGP(error, "in prolongSeed direction {}: layToSearch out of range {}", step, layToSearch);
        continue;
      }
      for (int jCe = firstIndex[layToSearch]; jCe < lastIndex[layToSearch]; ++jCe) {
        const CellCandidate& ccNext = cells[jCe];
        if ((dir == ExtendDirection::kInward && ccNext.secondTrackletIndex != nextLayerTrackletIndex) ||
            (dir == ExtendDirection::kOutward && ccNext.firstTrackletIndex != nextLayerTrackletIndex))
          continue;
        // --- CluID continuity checks ---
        if (dir == ExtendDirection::kInward) {
          if (ccNext.cluIDs[1] != refCell.cluIDs[0] || ccNext.cluIDs[2] != refCell.cluIDs[1]) {
            LOGP(error, "mismatch in CluIDs in prolongSeed in inward direction");
            continue;
          }
        } else {
          if (ccNext.cluIDs[0] != refCell.cluIDs[1] || ccNext.cluIDs[1] != refCell.cluIDs[2]) {
            LOGP(error, "mismatch in CluIDs in prolongSeed in outward direction");
            continue;
          }
        }
        // --- Chi2 checks ---
        // Compute chi2 only for the innermost cluster of the cells to be connected
        const auto& cellTst = (dir == ExtendDirection::kInward) ? refCell : ccNext;
        const auto& fitCurr = cellTst.trackFitFast;
        const auto& cluToAdd = cluArr[cellTst.cluIDs[0]];
        float chi2 = computeTrackToClusterChi2(fitCurr, cluToAdd);
#ifdef _CHI2_TUNING_MODE_
        auto mcTruth = getTrackMCTruthStatus(cellTst.cluIDs, cluArr);
        if (mcTruth.second && mcTruth.first == cluToAdd.getParticleID()) {
          (*dbgStream) << "chiProl"
                       << "iter=" << mCurIteration << "trc=" << ((NA6PTrackParCov&)fitCurr) << "trackNcl=" << fitCurr.getNHits() << "trackChi2=" << fitCurr.getChi2() << "chi2cl=" << chi2 << "\n";
        }
#endif
        if (chi2 > maxChi2TrClu) {
          if (mVerbose) {
            LOGP(info, "DBG prolong reject dir={} seedLayers=({},{}) cell={} nextClu={} chi2={} max={}",
                 dir == ExtendDirection::kInward ? "in" : "out",
                 seed.innerLayer, seed.outerLayer, jCe,
                 dir == ExtendDirection::kInward ? ccNext.cluIDs[0] : refCell.cluIDs[0],
                 chi2, maxChi2TrClu);
          }
          continue;
        }

        // create new prolonged track
        TrackCandidate extended = cand;
        if (dir == ExtendDirection::kInward) {
          extended.innerCellIndex = jCe;
          extended.innerLayer = ccNext.startingLayer;
          int nLay = cluArr[ccNext.cluIDs[0]].getLayer() - mLayerStart;
          extended.cluIDs[nLay] = ccNext.cluIDs[0];
        } else {
          extended.outerCellIndex = jCe;
          extended.outerLayer = ccNext.startingLayer;
          int nLay = cluArr[ccNext.cluIDs[2]].getLayer() - mLayerStart;
          extended.cluIDs[nLay] = ccNext.cluIDs[2];
        }
        next.push_back(std::move(extended));
        foundAnyProlongation = true;
      }
    }
    if (!foundAnyProlongation)
      break;
    current.swap(next);
  }
  // add remaining candidates from last layer
  for (auto& cand : current)
    result.push_back(cand);
  return result;
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::findRoads(const std::vector<std::pair<int, int>>& cneigh,
                              const std::vector<CellCandidate>& cells,
                              const std::vector<int>& firstIndex,
                              const std::vector<int>& lastIndex,
                              const std::vector<TrackletCandidate>& tracklets,
                              const std::vector<ClusterType>& cluArr,
                              std::vector<TrackCandidate>& trackCands,
                              float maxChi2TrClu)
{

  trackCands.clear();
  int nCellPairs = cneigh.size();
  for (int jPair = 0; jPair < nCellPairs; jPair++) {
    const auto& thisPair = cneigh[jPair];
    int jCe1 = thisPair.first;
    int jCe2 = thisPair.second;
    const CellCandidate& cc1 = cells[jCe1];
    const CellCandidate& cc2 = cells[jCe2];
    bool cc1Inner = (cc1.startingLayer < cc2.startingLayer);
    const int innerIndex = cc1Inner ? jCe1 : jCe2;
    const int outerIndex = cc1Inner ? jCe2 : jCe1;
    const CellCandidate& inner = cells[innerIndex];
    const CellCandidate& outer = cells[outerIndex];

    // build initial candidate
    std::vector<int> cluIDsFull(mNLayers, -1);
    int nClusIn = inner.cluIDs.size();
    int innerLay = 99;
    for (int jClu = 0; jClu < nClusIn; jClu++) {
      int cluID = inner.cluIDs[jClu];
      int nLay = cluArr[cluID].getLayer() - mLayerStart;
      cluIDsFull[nLay] = cluID;
      if (nLay < innerLay)
        innerLay = nLay;
    }
    int nClusOut = outer.cluIDs.size();
    int outerLay = 0;
    for (int jClu = 0; jClu < nClusOut; jClu++) {
      int cluID = outer.cluIDs[jClu];
      int nLay = cluArr[cluID].getLayer() - mLayerStart;
      if (cluIDsFull[nLay] != -1 && cluIDsFull[nLay] != cluID) {
        LOGP(error, "clu mismatch in layer {}", nLay);
        continue;
      }
      cluIDsFull[nLay] = cluID;
      if (nLay > outerLay)
        outerLay = nLay;
    }
    // consistency checks
    if (innerLay != inner.startingLayer) {
      LOGP(error, "mismatch on inner layer");
      continue;
    }
    if (outerLay != outer.startingLayer + 2) {
      LOGP(error, "mismatch on outer layer");
      continue;
    }
    TrackCandidate seed{innerLay, innerIndex, outerLay, outerIndex, cluIDsFull};
    auto inwardTracks = prolongSeed(seed, cells, firstIndex, lastIndex, cluArr, maxChi2TrClu, ExtendDirection::kInward);
    for (const auto& track : inwardTracks) {
      auto fullTracks = prolongSeed(track, cells, firstIndex, lastIndex, cluArr, maxChi2TrClu, ExtendDirection::kOutward);
      for (auto& finalTrack : fullTracks) {
        trackCands.push_back(std::move(finalTrack));
      }
    }
  }
  // remove duplicated roads
  std::sort(trackCands.begin(), trackCands.end(), [](const TrackCandidate& a, const TrackCandidate& b) {
    return a.cluIDs < b.cluIDs;
  });
  auto last = std::unique(trackCands.begin(), trackCands.end(), [](const TrackCandidate& a, const TrackCandidate& b) {
    return a.cluIDs == b.cluIDs;
  });
  trackCands.erase(last, trackCands.end());
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::fitAndSelectTracks(const std::vector<TrackCandidate>& trackCands,
                                       const std::vector<ClusterType>& cluArr,
                                       const NA6PMCTruthContainer& mcCluLabels,
                                       std::vector<TrackFitted>& tracks,
                                       const NA6PVertex* primVert,
                                       float maxChi2TrClu,
                                       int minNClu,
                                       float maxChi2NDF)
{

  std::vector<TrackFitted> fittedTracks;
  fittedTracks.reserve(trackCands.size());
  auto prop = Propagator::Instance();

  // fit all track candidates
  std::vector<int> cluIDsForfit;
  cluIDsForfit.reserve(mNLayers);
  for (const auto& cand : trackCands) {
    int nClus = cand.getCluIDsWOGaps(cluIDsForfit);
    if (nClus < minNClu) { // RSTODO do we really have such candidates at this stage?
      if (mVerbose) {
        LOGP(info, "DBG final skip: nClus={} min={} cluIDs={}", nClus, minNClu, formatClusterIDs(cand.cluIDs));
      }
      continue;
    }
    // genfit::Track fitTrack;
    NA6PTrack fitTrackFast;
    float chi2Inw = fitTrackPointsFast(cluIDsForfit, cluArr, fitTrackFast, maxChi2TrClu, maxChi2NDF);
#ifdef _CHI2_TUNING_MODE_
    auto mcTruth = mTrackFitter->getMCTruthStatus();
    if (mcTruth.second) {
      (*dbgStream) << "chiTracks"
                   << "iter=" << mCurIteration << "trc=" << ((NA6PTrackParCov&)fitTrackFast) << "chi2vec=" << mTrackFitter->getChi2Buffer() << "chi2Tot=" << chi2Inw << "\n";
    }
#endif
    if (chi2Inw < 0.f) {
      if (mVerbose) {
        LOGP(info, "DBG final reject: fit failed nClus={} cluIDs={}", nClus, formatClusterIDs(cluIDsForfit));
      }
      continue; // skip if fit failed
    }
    if (mDoOutwardPropagation) {
      // Outward state is auxiliary for matching; keep the inward track if it fails.
      auto& outer = fitTrackFast.getOuterParam();
      outer = fitTrackFast; // use inward as a seed
      float chi2Out = mTrackFitter->fitSeedOutward(outer, true, mUseLinRef);
      if (chi2Out >= 0.f) {
        if (mDoInwardRefit) {
          NA6PTrackParCov inwRefit = outer; // use outward fit as a seed
          float chi2Refit = mTrackFitter->fitSeedInward(inwRefit, true, mUseLinRef);
          if (chi2Refit >= 0.f) {
            static_cast<NA6PTrackParCov&>(fitTrackFast) = inwRefit;
            chi2Inw = chi2Refit;
          } else {
            fitTrackFast.getInwardParam().invalidate();
          }
        }
        if (prop->propagateToZ(outer, mZOutProp, mTrackFitter->getPropOpt())) {
          fitTrackFast.setChi2Out(chi2Out);
        } else {
          outer.invalidate();
        }
      } else {
        outer.invalidate();
      }
    } else {
      fitTrackFast.getOuterParam().invalidate();
    }
    if (mPropagateTracksToPrimaryVertex) {
      NA6PTrackParCov pvProp = fitTrackFast;
      if (prop->propagateToZ(pvProp, mPrimVertPos[2], mTrackFitter->getPropOpt())) {
        static_cast<NA6PTrackParCov&>(fitTrackFast) = pvProp;
      }
    }
    fitTrackFast.setChi2(chi2Inw);
    if (mDoTrackConstrainedToPrimVert && primVert) {
      mTrackFitter->constrainTrackToVertex(fitTrackFast, *primVert); // no exit on failure!
    } else {
      fitTrackFast.getVertexConstrainedParam().invalidate();
    }
    float chi2NDF = nClus < 3 ? 0.f : chi2Inw / (2 * nClus - mNDOF);
    fittedTracks.emplace_back(cand.innerLayer, cand.outerLayer, nClus, cand.cluIDs, std::move(fitTrackFast), chi2NDF);
  }
  // sort by chi2 in view of selection
  std::sort(fittedTracks.begin(), fittedTracks.end(), [](const TrackFitted& a, const TrackFitted& b) { // RSTODO do we need really in-place sorting or just reindexing
    return a.chi2ndf < b.chi2ndf;
  });
  // treat tracks with shared clus: in case of sharing keep the track with lower chi2
  for (auto& track : fittedTracks) {
    int nShared = 0;
    for (int cluID : track.cluIDs) {
      if (cluID >= 0 && mIsClusterUsed[cluID]) {
        nShared++;
      }
    }
    if (nShared > mMaxSharedClusters) { // tracks should not be stored in the list of selected tracks
      continue;
    }
    // mark the clusters of the selected track as used
    // and register clusters in the final track. RSTODO currently there is a redundancy in clusters indices of the TrackFitted and those of the NA6PTrack
    for (int cluID : track.cluIDs) {
      if (cluID >= 0) {
        mIsClusterUsed[cluID] = true;
        track.trackFitFast.addCluster(&cluArr[cluID]);
      }
    }
    // assign overall MC label to the track
    int nClus = track.cluIDs.size();
    std::vector<std::pair<int, int>> counts;
    for (int jClu = 0; jClu < nClus; jClu++) {
      int cluID = track.cluIDs[jClu];
      if (cluID >= 0) {
        const auto& clu = cluArr[cluID];
        int idPartClu = clu.getParticleID();
        bool found = false;
        for (auto& p : counts) {
          if (p.first == idPartClu) {
            ++p.second;
            found = true;
            break;
          }
        }
        if (!found) {
          counts.push_back({idPartClu, 1});
        }
      }
    }
    int idPartTrack = -9999999;
    int maxCount = 0;
    for (const auto& p : counts) {
      if (p.second > maxCount) {
        maxCount = p.second;
        idPartTrack = p.first;
      }
    }
    // assign negative label to tracks with misassociations
    if (counts.size() > 1 && idPartTrack > 0) {
      idPartTrack *= -1;
    }
    track.trackFitFast.setParticleID(idPartTrack);
    // assign label using MCTruthContainer
    std::vector<std::pair<NA6PMCComposedLabel, int>> countLabs;
    for (int jClu = 0; jClu < nClus; jClu++) {
      int cluID = track.cluIDs[jClu];
      if (cluID >= 0) {
        const auto& clu = cluArr[cluID];
        int originalID = clu.getClusterIndex();
        std::span labels = mcCluLabels.getLabels(originalID);
        int nLabels = labels.size();
        for (int jLab = 0; jLab < nLabels; jLab++) {
          NA6PMCComposedLabel lbl = labels[jLab];
          bool found = false;
          for (auto& p : countLabs) {
            if (p.first == lbl) {
              ++p.second;
              found = true;
              break;
            }
          }
          if (!found) {
            countLabs.push_back({lbl, 1});
          }
        }
      }
    }
    NA6PMCComposedLabel lblTrack;
    lblTrack.unset();
    int maxCountLabs = 0;
    for (const auto& p : countLabs) {
      if (p.second > maxCountLabs) {
        maxCountLabs = p.second;
        lblTrack = p.first;
      }
    }
    // assign fake label to tracks with misassociations
    if (lblTrack.isSet()) {
      for (int jClu = 0; jClu < nClus; jClu++) {
        int cluID = track.cluIDs[jClu];
        if (cluID >= 0) {
          const auto& clu = cluArr[cluID];
          int originalID = clu.getClusterIndex();
          std::span labels = mcCluLabels.getLabels(originalID);
          int nLabels = labels.size();
          bool found = false;
          for (int jLab = 0; jLab < nLabels; jLab++) {
            NA6PMCComposedLabel lbl = labels[jLab];
            if (lbl == lblTrack) {
              found = true;
              break;
            }
          }
          if (!found) {
            lblTrack.setFakeFlag();
            break;
          }
        }
      }
    }
    track.trackFitFast.setMCLabel(lblTrack);
    tracks.push_back(std::move(track));
  }
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::findTracks(std::vector<ClusterType>& cluArr,
                               const NA6PMCTruthContainer& mcCluLabels,
                               const NA6PVertex* primVert)
{
  mNDOF = Propagator::Instance()->getNDOFTrack();
  mFinalTracks.clear();
  uint nClus = cluArr.size();
  NA6PVertex vertWithDefaultCov;
  if (primVert) {
    mPrimVertPos[0] = primVert->getX();
    mPrimVertPos[1] = primVert->getY();
    mPrimVertPos[2] = primVert->getZ();
    if (primVert->getSigmaX2() <= 0.f ||
        primVert->getSigmaY2() <= 0.f ||
        primVert->getSigmaZ2() <= 0.f) {
      if (mDoTrackConstrainedToPrimVert)
        LOGP(warn, "Covariance matrix of vertex not set. Use default values");
      // Assign default values to covariance matrix (based on beam width and target thickness)
      float defaultSigmaXY = 0.02f;
      float defaultSigmaZ = 0.15f;
      vertWithDefaultCov = *primVert;
      vertWithDefaultCov.setSigmaX(defaultSigmaXY);
      vertWithDefaultCov.setSigmaY(defaultSigmaXY);
      vertWithDefaultCov.setSigmaZ(defaultSigmaZ);
      vertWithDefaultCov.setSigmaXY(0.f);
      vertWithDefaultCov.setSigmaXZ(0.f);
      vertWithDefaultCov.setSigmaYZ(0.f);
      primVert = &vertWithDefaultCov;
    }
  } else {
    mPrimVertPos[0] = mPrimVertPos[1] = mPrimVertPos[2] = 0.0f;
  }

  LOGP(info, "Process event with nClusters {}, primary vertex in z = {} cm", nClus, mPrimVertPos[2]);
  mIsClusterUsed.resize(nClus);
  for (uint jClu = 0; jClu < nClus; jClu++)
    mIsClusterUsed[jClu] = false;
  std::vector<int> firstCluPerLay;
  std::vector<int> lastCluPerLay;
  sortClustersByLayerAndEta(cluArr, firstCluPerLay, lastCluPerLay);

  std::vector<TrackletCandidate> foundTracklets;
  std::vector<CellCandidate> foundCells;
  std::vector<std::pair<int, int>> cellsNeighbours;
  std::vector<TrackCandidate> trackCandidates;
  std::vector<TrackFitted> iterationTracks;
  for (int jIteration = 0; jIteration < mNIterationsCA; ++jIteration) {
    mCurIteration = jIteration;
    if (mVerbose) {
      LOGP(info, " -> Iteration {} <-", jIteration);
    }
    foundTracklets.clear();
    foundCells.clear();
    cellsNeighbours.clear();
    trackCandidates.clear();
    iterationTracks.clear();
    computeLayerTracklets(cluArr, firstCluPerLay, lastCluPerLay, foundTracklets, mMaxDeltaThetaTrackletsCA[jIteration], mMaxDeltaPhiTrackletsCA[jIteration]);
    std::vector<int> firstTrklPerLay;
    std::vector<int> lastTrklPerLay;
    sortTrackletsByLayerAndIndex(foundTracklets, firstTrklPerLay, lastTrklPerLay);
    if (mVerbose)
      printStats(foundTracklets, cluArr, foundCells, "tracklets");
    //
    computeLayerCells(foundTracklets, firstTrklPerLay, lastTrklPerLay, cluArr, foundCells, mMaxDeltaTanLCellsCA[jIteration], mMaxDeltaPhiCellsCA[jIteration], mMaxDeltaPxPzCellsCA[jIteration], mMaxDeltaPyPzCellsCA[jIteration], mMaxChi2TrClCellsCA[jIteration], mMaxChi2ndfCellsCA[jIteration]);
    std::vector<int> firstCellPerLay;
    std::vector<int> lastCellPerLay;
    sortCellsByLayerAndIndex(foundCells, firstCellPerLay, lastCellPerLay);
    if (mVerbose)
      printStats(foundCells, cluArr, foundCells, "cells");
    //
    findCellsNeighbours(foundCells, firstCellPerLay, lastCellPerLay, cellsNeighbours, cluArr, mMaxChi2TrClCellsCA[jIteration]);
    if (mVerbose)
      printStats(cellsNeighbours, cluArr, foundCells, "cell pairs");
    //
    findRoads(cellsNeighbours, foundCells, firstCellPerLay, lastCellPerLay, foundTracklets, cluArr, trackCandidates, mMaxChi2TrClCellsCA[jIteration]);
    if (mVerbose)
      printStats(trackCandidates, cluArr, foundCells, "track candidates");
    //
    fitAndSelectTracks(trackCandidates, cluArr, mcCluLabels, iterationTracks, primVert, mMaxChi2TrClCellsCA[jIteration], mMinNClusTracksCA[jIteration], mMaxChi2ndfTracksCA[jIteration]);
    if (mVerbose) {
      printStats(iterationTracks, cluArr, foundCells, "selected tracks");
      printStats(iterationTracks, cluArr, foundCells, "selected tracks", 5);
      printStats(iterationTracks, cluArr, foundCells, "selected tracks", 4);
    }
    for (auto& track : iterationTracks) {
      track.trackFitFast.setCAIteration(jIteration);
    }
    mFinalTracks.insert(mFinalTracks.end(), iterationTracks.begin(), iterationTracks.end());
  }
}

//______________________________________________________________________
std::vector<NA6PTrack> NA6PTrackerCA::getTracks()
{
  std::vector<NA6PTrack> trackArr;
  trackArr.reserve(mFinalTracks.size());
  for (auto& track : mFinalTracks)
    trackArr.push_back(std::move(track.trackFitFast));
  return trackArr;
}

//______________________________________________________________________
template <typename ClusterType>
std::vector<std::pair<ClusterType, ClusterType>> NA6PTrackerCA::findTracklets(int jFirstLay, int jLastLay, std::vector<ClusterType>& cluArr, const NA6PVertex* primVert)
{
  if (primVert) {
    mPrimVertPos[0] = primVert->getX();
    mPrimVertPos[1] = primVert->getY();
    mPrimVertPos[2] = primVert->getZ();
  } else {
    mPrimVertPos[0] = mPrimVertPos[1] = mPrimVertPos[2] = 0.0f;
  }
  int backupStart = mLayerStart;
  int backupLay = mNLayers;
  mLayerStart = jFirstLay;
  mNLayers = jLastLay - jFirstLay + 1;
  std::vector<int> firstCluPerLay;
  std::vector<int> lastCluPerLay;
  sortClustersByLayerAndEta(cluArr, firstCluPerLay, lastCluPerLay);
  std::vector<TrackletCandidate> foundTracklets;
  float cutDeltaTheta = mMaxDeltaThetaTrackletsCA[0];
  float cutDeltaPhi = mMaxDeltaPhiTrackletsCA[0];
  for (int jIteration = 0; jIteration < mNIterationsCA; ++jIteration) {
    if (mMaxDeltaThetaTrackletsCA[jIteration] > cutDeltaTheta)
      cutDeltaTheta = mMaxDeltaThetaTrackletsCA[jIteration];
    if (mMaxDeltaPhiTrackletsCA[jIteration] > cutDeltaPhi)
      cutDeltaPhi = mMaxDeltaPhiTrackletsCA[jIteration];
  }
  computeLayerTracklets(cluArr, firstCluPerLay, lastCluPerLay, foundTracklets, cutDeltaTheta, cutDeltaPhi);
  const int nTracklets = foundTracklets.size();
  std::vector<std::pair<ClusterType, ClusterType>> trackletLines;
  trackletLines.reserve(nTracklets);
  for (auto& trkl : foundTracklets) {
    auto clu0 = cluArr[trkl.firstClusterIndex];
    auto clu1 = cluArr[trkl.secondClusterIndex];
    trackletLines.emplace_back(std::make_pair(clu0, clu1));
  }
  mLayerStart = backupStart;
  mNLayers = backupLay;
  return trackletLines;
}

//______________________________________________________________________
template <typename T, typename ClusterType>
void NA6PTrackerCA::printStats(const std::vector<T>& candidates,
                               const std::vector<ClusterType>& cluArr,
                               const std::vector<CellCandidate>& cells,
                               const std::string& label,
                               int requiredClus)
{

  int nFound = candidates.size();
  if (requiredClus < 0)
    LOGP(info, "Number of {} = {}", label.c_str(), nFound);
  int nGood = 0;
  int nSelected = 0;
  float aveClus = 0;
  for (int j = 0; j < nFound; ++j) {
    const auto& tr = candidates[j];
    int nClus = -1;
    int idClus[10];
    int startLay = -1;
    if constexpr (std::is_same_v<T, TrackletCandidate>) {
      nClus = 2;
      idClus[0] = tr.firstClusterIndex;
      idClus[1] = tr.secondClusterIndex;
      startLay = tr.startingLayer;
    } else if constexpr (std::is_same_v<T, std::pair<int, int>>) {
      int jCe1 = tr.first;
      int jCe2 = tr.second;
      CellCandidate cc1 = cells[jCe1];
      CellCandidate cc2 = cells[jCe2];
      nClus = cc1.cluIDs.size() + 1;
      for (int jClu = 0; jClu < nClus; jClu++) {
        if (jClu < nClus - 1)
          idClus[jClu] = cc1.cluIDs[jClu];
        else
          idClus[jClu] = cc2.cluIDs[2];
      }
      startLay = std::min(cc1.startingLayer, cc2.startingLayer);
    } else if constexpr (std::is_same_v<T, TrackCandidate> || std::is_same_v<T, TrackFitted>) {
      nClus = tr.cluIDs.size();
      for (int jClu = 0; jClu < nClus; jClu++)
        idClus[jClu] = tr.cluIDs[jClu];
      startLay = tr.innerLayer;
    } else {
      nClus = tr.cluIDs.size();
      for (int jClu = 0; jClu < nClus; jClu++)
        idClus[jClu] = tr.cluIDs[jClu];
      startLay = tr.startingLayer;
    }
    int idPartTrack = -9999999;
    int nTrueClus = 0;
    for (int jClu = 0; jClu < nClus; jClu++) {
      int cluID = idClus[jClu];
      if (cluID >= 0) {
        ++nTrueClus;
        const auto& clu = cluArr[cluID];
        int jLay = clu.getLayer() - mLayerStart;
        if (jClu == 0 && jLay != startLay)
          LOGP(error, "mismatch in {} layers: {} {}", label.c_str(), jLay, startLay);
        int idPartClu = clu.getParticleID();
        if (jClu == 0)
          idPartTrack = idPartClu;
        else if (idPartClu != idPartTrack && idPartTrack > 0)
          idPartTrack *= (-1);
      }
    }
    if (requiredClus > 0 && nTrueClus != requiredClus)
      continue;
    nSelected++;
    if (idPartTrack > 0)
      nGood++;
    aveClus += nTrueClus;
  }
  if (requiredClus > 0) {
    if (nSelected > 0) {
      LOGP(info, "{} with {} clus: Fraction of good = {} / {} = {}  --- average clus = {}",
           label.c_str(), requiredClus, nGood, nSelected,
           static_cast<float>(nGood) / static_cast<float>(nSelected),
           aveClus / static_cast<float>(nSelected));
    } else {
      LOGP(info, "No {} having {} clus", label.c_str(), requiredClus);
    }
  } else {
    if (nFound > 0)
      LOGP(info, "Fraction of good {} = {} / {} = {}  --- average clus = {}",
           label.c_str(), nGood, nFound,
           static_cast<float>(nGood) / static_cast<float>(nFound),
           aveClus / static_cast<float>(nFound));
  }
}

template void NA6PTrackerCA::findTracks<NA6PMuonSpecCluster>(std::vector<NA6PMuonSpecCluster>&, const NA6PMCTruthContainer&, const NA6PVertex*);
template void NA6PTrackerCA::findTracks<NA6PVerTelCluster>(std::vector<NA6PVerTelCluster>&, const NA6PMCTruthContainer&, const NA6PVertex*);
template std::vector<std::pair<NA6PMuonSpecCluster, NA6PMuonSpecCluster>> NA6PTrackerCA::findTracklets<NA6PMuonSpecCluster>(int, int, std::vector<NA6PMuonSpecCluster>&, const NA6PVertex*);
template std::vector<std::pair<NA6PVerTelCluster, NA6PVerTelCluster>> NA6PTrackerCA::findTracklets<NA6PVerTelCluster>(int, int, std::vector<NA6PVerTelCluster>&, const NA6PVertex*);
