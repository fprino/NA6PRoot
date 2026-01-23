// NA6PCCopyright

#include <fmt/format.h>
#include <iostream>
#include <fairlogger/Logger.h>
#include <TGeoManager.h>
#include <TSystem.h>
#include "MagneticField.h"
#include "NA6PFastTrackFitter.h"
#include "NA6PRecoParam.h"
#include "NA6PTrackerCA.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PVerTelCluster.h"

ClassImp(NA6PTrackerCA)

  //______________________________________________________________________

  //______________________________________________________________________

  NA6PTrackerCA::NA6PTrackerCA() : mLayerStart{0},
                                   mNLayers{5},
                                   mPrimVertPos{0., 0., 0.},
                                   mTrackFitter{nullptr},
                                   mIsClusterUsed{false},
                                   mMaxSharedClusters{0},
                                   mVerbose{false},
                                   mNIterationsCA{2}
{
  mTrackFitter = new NA6PFastTrackFitter();
  mTrackFitter->setNLayers(mNLayers);
  mTrackFitter->setPropagateToPrimaryVertex(false);
}
//______________________________________________________________________

NA6PTrackerCA::~NA6PTrackerCA()
{
  if (mTrackFitter)
    delete mTrackFitter;
  mTrackFitter = nullptr;
}
//______________________________________________________________________
void NA6PTrackerCA::setNLayers(int n)
{
  mNLayers = n;
  if (mTrackFitter)
    mTrackFitter->setNLayers(n);
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
                                       double maxDeltaThetaTracklets,
                                       double maxDeltaPhiTracklets,
                                       double maxDeltaTanLCells,
                                       double maxDeltaPhiCells,
                                       double maxDeltaPxPzCells,
                                       double maxDeltaPyPzCells,
                                       double maxChi2TrClCells,
                                       double maxChi2ndfCells,
                                       double maxChi2ndfTracks,
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

void NA6PTrackerCA::configureFromRecoParam(const std::string& filename)
{
  if (filename != "") {
    na6p::conf::ConfigurableParamHelper<NA6PRecoParam>::updateFromFile(filename);
  }
  const auto& param = NA6PRecoParam::Instance();
  setNumberOfIterations(param.nIterationsTrackerCA);
  setNLayers(param.nLayers);
  for (int jIter = 0; jIter < mNIterationsCA; ++jIter) {
    setIterationParams(jIter,
                       param.maxDeltaThetaTrackletsCA[jIter],
                       param.maxDeltaPhiTrackletsCA[jIter],
                       param.maxDeltaTanLCellsCA[jIter],
                       param.maxDeltaPhiCellsCA[jIter],
                       param.maxDeltaPxPzCellsCA[jIter],
                       param.maxDeltaPyPzCellsCA[jIter],
                       param.maxChi2TrClCellsCA[jIter],
                       param.maxChi2ndfCellsCA[jIter],
                       param.maxChi2ndfTracksCA[jIter],
                       param.minNClusTracksCA[jIter]);
  }
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

  if (!mTrackFitter) {
    LOGP(error, "Fitter not yet instantiated");
    return false;
  }
  return mTrackFitter->loadGeometry(filename, geoname);
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::sortClustersByLayerAndEta(std::vector<ClusterType>& cluArr,
                                              std::vector<int>& firstIndex,
                                              std::vector<int>& lastIndex)
{

  // count clus per layer
  std::vector<int> count(mNLayers, 0);
  for (const auto& clu : cluArr) {
    int jDet = clu.getDetectorID();
    int jLay = clu.getLayer() - mLayerStart;
    if (jLay >= 0 && jLay < mNLayers) {
      count[jLay]++;
    }
  }
  // starting offset for each layer
  firstIndex.resize(mNLayers);
  lastIndex.resize(mNLayers);
  firstIndex[0] = 0;
  lastIndex[0] = count[0];
  for (int jLay = 1; jLay < mNLayers; jLay++) {
    firstIndex[jLay] = firstIndex[jLay - 1] + count[jLay - 1];
    lastIndex[jLay] = firstIndex[jLay] + count[jLay];
  }
  // Create a reordered vector based on layer grouping
  std::vector<int> countReord = firstIndex;
  std::vector<ClusterType> reordered(cluArr.size());
  for (const auto& clu : cluArr) {
    int jLay = clu.getLayer() - mLayerStart;
    if (jLay >= 0 && jLay < mNLayers)
      reordered[countReord[jLay]++] = clu;
  }
  cluArr = std::move(reordered);
  // sort by theta within each layer (use z^2/r^2 as a proxy of theta to avoid sqrt and atan)
  const double pvx = mPrimVertPos[0];
  const double pvy = mPrimVertPos[1];
  const double pvz = mPrimVertPos[2];
  for (int jLay = 0; jLay < mNLayers; jLay++) {
    auto first = cluArr.begin() + firstIndex[jLay];
    auto last = cluArr.begin() + lastIndex[jLay];
    std::sort(first, last, [pvx, pvy, pvz](const ClusterType& a, const ClusterType& b) {
      double xa = a.getX() - pvx;
      double ya = a.getY() - pvy;
      double za = a.getZ() - pvz;
      double xb = b.getX() - pvx;
      double yb = b.getY() - pvy;
      double zb = b.getZ() - pvz;
      double r2a = xa * xa + ya * ya;
      double r2b = xb * xb + yb * yb;
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
                                          double deltaThetaMax,
                                          double deltaPhiMax)
{

  tracklets.clear();
  const double pvx = mPrimVertPos[0];
  const double pvy = mPrimVertPos[1];
  const double pvz = mPrimVertPos[2];

  for (int iLayer = 0; iLayer < mNLayers - 1; ++iLayer) {
    auto layerBegin = cluArr.begin() + firstIndex[iLayer + 1];
    auto layerEnd = cluArr.begin() + lastIndex[iLayer + 1];
    for (int jClu1 = firstIndex[iLayer]; jClu1 < lastIndex[iLayer]; ++jClu1) {
      if (mIsClusterUsed[jClu1])
        continue;
      const ClusterType& clu1 = cluArr[jClu1];
      double x1 = clu1.getX() - pvx;
      double y1 = clu1.getY() - pvy;
      double z1 = clu1.getZ() - pvz;
      double r1 = std::sqrt(x1 * x1 + y1 * y1);
      double theta1 = std::atan2(z1, r1);
      double phi1 = std::atan2(y1, x1);
      double tanth2Min = std::max(0., std::tan(theta1 - 1.2 * deltaThetaMax)); // 1.2 is a safety margin, the tan(theta) range is limited at zero to protect for the case in which theta1 - 1.2 * deltaThetaMax is <0, which would give a negative tangent. The minimum acceptable value for theta is 0
      double tanth2Max = std::tan(theta1 + 1.2 * deltaThetaMax);
      auto lower = std::partition_point(layerBegin, layerEnd,
                                        [pvx, pvy, pvz, tanth2Min](const ClusterType& clu) {
                                          double x = clu.getX() - pvx;
                                          double y = clu.getY() - pvy;
                                          double z = clu.getZ() - pvz;
                                          double r2 = x * x + y * y;
                                          double tan2 = z * z / r2;
                                          return tan2 < tanth2Min * tanth2Min;
                                        });

      auto upper = std::partition_point(layerBegin, layerEnd,
                                        [pvx, pvy, pvz, tanth2Max](const ClusterType& clu) {
                                          double x = clu.getX() - pvx;
                                          double y = clu.getY() - pvy;
                                          double z = clu.getZ() - pvz;
                                          double r2 = x * x + y * y;
                                          double tan2 = z * z / r2;
                                          return tan2 <= tanth2Max * tanth2Max;
                                        });
      int lowerIdx = std::distance(cluArr.begin(), lower);
      int upperIdx = std::distance(cluArr.begin(), upper);
      for (int jClu2 = lowerIdx; jClu2 < upperIdx; ++jClu2) {
        if (mIsClusterUsed[jClu2])
          continue;
        const ClusterType& clu2 = cluArr[jClu2];
        double x2 = clu2.getX() - pvx;
        double y2 = clu2.getY() - pvy;
        double z2 = clu2.getZ() - pvz;
        double r2 = std::sqrt(x2 * x2 + y2 * y2);
        double theta2 = std::atan2(z2, r2);
        double phi2 = std::atan2(y2, x2);
        double dphi = phi2 - phi1;
        if (dphi > M_PI)
          dphi -= 2 * M_PI;
        else if (dphi < -M_PI)
          dphi += 2 * M_PI;
        if (std::abs(theta2 - theta1) < deltaThetaMax && std::abs(dphi) < deltaPhiMax) {
          double phi = std::atan2(y2 - y1, x2 - x1);
          double tanL = (z2 - z1) / (r2 - r1);
          double pxpz = (x2 - x1) / (z2 - z1);
          double pypz = (y2 - y1) / (z2 - z1);
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
                                      double deltaTanLMax,
                                      double deltaPhiMax,
                                      double deltaPxPzMax,
                                      double deltaPyPzMax,
                                      double maxChi2TrClu,
                                      double maxChi2NDF)
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
        const double deltaTanLambda = std::abs(trkl2.tanL - trkl1.tanL);
        double dphi = trkl2.phi - trkl1.phi;
        if (dphi > M_PI)
          dphi -= 2 * M_PI;
        else if (dphi < -M_PI)
          dphi += 2 * M_PI;
        double deltapxpz = std::abs(trkl2.pxpz - trkl1.pxpz);
        double deltapypz = std::abs(trkl2.pypz - trkl1.pypz);
        if (deltapypz < deltaPyPzMax && deltapxpz < deltaPxPzMax && deltaTanLambda < deltaTanLMax && std::abs(dphi) < deltaPhiMax) {
          std::array<int, 3> cluIDs = {trkl1.firstClusterIndex, trkl2.firstClusterIndex, trkl2.secondClusterIndex};
          NA6PTrack fitTrackFast;
          //	  genfit::Track fitTrack;
          if (fitTrackPointsFast(std::vector<int>(cluIDs.begin(), cluIDs.end()), cluArr, fitTrackFast, maxChi2TrClu, maxChi2NDF)) {
            cells.emplace_back(iLayer, jTrkl1, jTrkl2, std::move(cluIDs), fitTrackFast);
          }
        }
      }
    }
  }
}

//______________________________________________________________________

template <typename ClusterType>
double NA6PTrackerCA::computeTrackToClusterChi2(const NA6PTrack& track,
                                                const ClusterType& clu)
{

  double meas[2] = {clu.getYTF(), clu.getZTF()}; // ideal cluster coordinate, tracking (AliExtTrParam frame)
  double measErr2[3] = {clu.getSigYY(), clu.getSigYZ(), clu.getSigZZ()};
  NA6PTrack copyToProp = track;
  copyToProp.propagateToZBxByBz(clu.getZ()); // no material correction temporarily
  double cluchi2 = copyToProp.getTrackExtParam().getPredictedChi2(meas, measErr2);
  return cluchi2;
}

//______________________________________________________________________

template <typename ClusterType>
bool NA6PTrackerCA::fitTrackPointsFast(const std::vector<int>& cluIDs,
                                       const std::vector<ClusterType>& cluArr,
                                       NA6PTrack& fitTrack,
                                       double maxChi2TrClu,
                                       double maxChi2NDF)
{

  int nClus = cluIDs.size();
  mTrackFitter->cleanupAndStartFit();
  mTrackFitter->setMaxChi2Cl(maxChi2TrClu);
  std::vector<const ClusterType*> clusters;
  clusters.reserve(nClus);
  for (int jClu = 0; jClu < nClus; jClu++) {
    int cluID = cluIDs[jClu];
    const auto& clu = cluArr[cluID];
    int nLay = clu.getLayer() - mLayerStart;
    clusters.push_back(&clu); // store the pointer for later use
    mTrackFitter->addCluster(nLay, clu);
  }

  std::unique_ptr<NA6PTrack> fitTrackPtr(mTrackFitter->fitTrackPoints());
  if (!fitTrackPtr)
    return false;

  double chi2ndf = fitTrackPtr->getNormChi2();
  if (chi2ndf > maxChi2NDF)
    return false;

  fitTrack = *fitTrackPtr;
  return true;
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::findCellsNeighbours(const std::vector<CellCandidate>& cells,
                                        const std::vector<int>& firstIndex,
                                        const std::vector<int>& lastIndex,
                                        std::vector<std::pair<int, int>>& cneigh,
                                        const std::vector<ClusterType>& cluArr,
                                        double maxChi2TrClu)
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

        if (mChi2TrackCluOpt == kBothWoOutwFit) {
          int jClu2 = cell2.cluIDs[2];
          const auto& clu2 = cluArr[jClu2];
          double cluchi2a = computeTrackToClusterChi2(cell1.trackFitFast, clu2);
          if (cluchi2a > maxChi2TrClu)
            continue;
        }

        int jClu0 = cell1.cluIDs[0];
        const auto& clu0 = cluArr[jClu0];
        double cluchi2b = computeTrackToClusterChi2(cell2.trackFitFast, clu0);
        if (cluchi2b > maxChi2TrClu)
          continue;

        if (mChi2TrackCluOpt == kBothWithOutwFit) {
          int nClus = cell1.cluIDs.size();
          mTrackFitter->cleanupAndStartFit();
          mTrackFitter->setMaxChi2Cl(maxChi2TrClu);
          std::vector<const ClusterType*> clusters;
          clusters.reserve(nClus);
          for (int jClu = 0; jClu < nClus; jClu++) {
            int cluID = cell1.cluIDs[jClu];
            const auto& clu = cluArr[cluID];
            int nLay = clu.getLayer() - mLayerStart;
            clusters.push_back(&clu); // store the pointer for later use
            mTrackFitter->addCluster(nLay, clu);
          }

          std::unique_ptr<NA6PTrack> fitTrackOutw(mTrackFitter->fitTrackPointsOutward());
          if (!fitTrackOutw) {
            LOGP(info, "Outward track fit failed");
            continue;
          }
          double chi2ndf = fitTrackOutw->getNormChi2();
          int jClu2 = cell2.cluIDs[2];
          const auto& clu2 = cluArr[jClu2];
          double cluchi2a = computeTrackToClusterChi2(*fitTrackOutw, clu2);
          if (cluchi2a > maxChi2TrClu)
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
                                                       double maxChi2TrClu,
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
        if (mChi2TrackCluOpt != kOnlyInward) {
          const auto& fitNext = ccNext.trackFitFast;
          int cluRefIndex = (dir == ExtendDirection::kInward) ? 2 : 0;
          const auto& cluRef = cluArr[refCell.cluIDs[cluRefIndex]];
          double chi2a = computeTrackToClusterChi2(fitNext, cluRef);
          if (chi2a > maxChi2TrClu)
            continue;
          const auto& fitRef = refCell.trackFitFast;
          int cluNextIndex = (dir == ExtendDirection::kInward) ? 0 : 2;
          const auto& cluNext = cluArr[ccNext.cluIDs[cluNextIndex]];
          double chi2b = computeTrackToClusterChi2(fitRef, cluNext);
          if (chi2b > maxChi2TrClu)
            continue;
        } else {
          // Compute chi2 only for the innermost cluster of the cells to be connected
          const auto& fitCurr = (dir == ExtendDirection::kInward) ? refCell.trackFitFast : ccNext.trackFitFast;
          const auto& cluToAdd = (dir == ExtendDirection::kInward) ? cluArr[ccNext.cluIDs[0]] : cluArr[refCell.cluIDs[0]];
          double chi2 = computeTrackToClusterChi2(fitCurr, cluToAdd);
          if (chi2 > maxChi2TrClu)
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
                              double maxChi2TrClu)
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
                                       std::vector<TrackFitted>& tracks,
                                       double maxChi2TrClu,
                                       int minNClu,
                                       double maxChi2NDF)
{

  std::vector<TrackFitted> fittedTracks;
  fittedTracks.reserve(trackCands.size());

  // fit all track candidates
  for (const auto& cand : trackCands) {
    int nClus = 0;
    std::vector<int> cluIDsForfit;
    for (int cluID : cand.cluIDs) {
      if (cluID >= 0) { // cluID == -1 correspond to no clu on that layer
        nClus++;
        cluIDsForfit.push_back(cluID);
      }
    }
    // fit the track
    //    genfit::Track fitTrack;
    NA6PTrack fitTrackFast;
    //    bool fitSuccess = fitTrackPoints(cluIDsForfit,cluArr,fitter,fitTrack,9999999.,maxChi2NDF);
    bool fitSuccess = fitTrackPointsFast(cluIDsForfit, cluArr, fitTrackFast, maxChi2TrClu, maxChi2NDF);
    if (!fitSuccess)
      continue;                      // skip if fit failed
    nClus = fitTrackFast.getNHits(); // some clusters could have been rejected in the fit
    if (nClus < minNClu)
      continue;
    // const genfit::AbsTrackRep* rep = fitTrack.getTrackRep(0);
    // double chi2 = fitTrack.getFitStatus(rep)->getChi2();
    // double ndf = fitTrack.getFitStatus(rep)->getNdf();
    // double chi2ndf = (ndf > 0) ? chi2 / ndf : 9999.0;
    double chi2ndf = fitTrackFast.getNormChi2();
    fittedTracks.emplace_back(cand.innerLayer, cand.outerLayer, nClus, cand.cluIDs, std::move(fitTrackFast), chi2ndf);
  }
  // sort by chi2 in view of selection
  std::sort(fittedTracks.begin(), fittedTracks.end(), [](const TrackFitted& a, const TrackFitted& b) {
    return a.chi2ndf < b.chi2ndf;
  });
  // treat tracks with shared clus: in case of sharing keep the track with lower chi2
  for (auto& track : fittedTracks) {
    int nShared = 0;
    for (int cluID : track.cluIDs) {
      if (cluID >= 0 && mIsClusterUsed[cluID])
        nShared++;
    }
    if (nShared > mMaxSharedClusters) {
      // tracks should not be stored in the list of selected tracks
      continue;
    }
    // mark the clusters of the selected track as used
    for (int cluID : track.cluIDs) {
      if (cluID >= 0)
        mIsClusterUsed[cluID] = true;
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
    if (counts.size() > 1 && idPartTrack > 0)
      idPartTrack *= -1;
    track.trackFitFast.setParticleID(idPartTrack);
    tracks.push_back(std::move(track));
  }
}

//______________________________________________________________________

template <typename ClusterType>
void NA6PTrackerCA::findTracks(std::vector<ClusterType>& cluArr, TVector3 primVert)
{

  mFinalTracks.clear();
  uint nClus = cluArr.size();
  mPrimVertPos[0] = primVert.X();
  mPrimVertPos[1] = primVert.Y();
  mPrimVertPos[2] = primVert.Z();
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
    if (mVerbose)
      LOGP(info, " -> Iteration {} <-", jIteration);
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
    fitAndSelectTracks(trackCandidates, cluArr, iterationTracks, mMaxChi2TrClCellsCA[jIteration], mMinNClusTracksCA[jIteration], mMaxChi2ndfTracksCA[jIteration]);
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
  double aveClus = 0;
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
           label.c_str(), requiredClus, nGood, nSelected, (double)nGood / (double)nSelected, aveClus / (double)nSelected);
    } else {
      LOGP(info, "No {} having {} clus", label.c_str(), requiredClus);
    }
  } else {
    if (nFound > 0)
      LOGP(info, "Fraction of good {} = {} / {} = {}  --- average clus = {}",
           label.c_str(), nGood, nFound, (double)nGood / (double)nFound, aveClus / (double)nFound);
  }
}

template void NA6PTrackerCA::findTracks<NA6PMuonSpecCluster>(std::vector<NA6PMuonSpecCluster>&, TVector3);
template void NA6PTrackerCA::findTracks<NA6PVerTelCluster>(std::vector<NA6PVerTelCluster>&, TVector3);
