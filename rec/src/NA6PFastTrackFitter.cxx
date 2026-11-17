// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TSystem.h>
#include "NA6PVertex.h"
#include "MagneticField.h"
#include "NA6PPhysConst.h"
#include "NA6PFastTrackFitter.h"
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PLine.h"
#include "Seeder.h"

NA6PFastTrackFitter::NA6PFastTrackFitter()
{
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    Propagator::loadField();
  }
}

void NA6PFastTrackFitter::printClusters() const
{
  for (int jLay = 0; jLay < getNLayers(); ++jLay) {
    if (mClusters[jLay])
      LOGP(info, "Layer {} Cluster position = {} {} {}", jLay, mClusters[jLay]->getX(), mClusters[jLay]->getY(), mClusters[jLay]->getZ());
    else
      LOGP(info, "Layer {} No cluster", jLay);
  }
}

void NA6PFastTrackFitter::setSeed(const float* pos, const float* mom, int charge) // RSTODO
{
  if (!pos || !mom) {
    LOGP(error, "Null pointers passer seed");
    return;
  }
  mSeed.initParam(pos, mom, charge);
}

bool NA6PFastTrackFitter::constrainTrackToVertex(NA6PTrack& trc, const NA6PVertex& pv)
{
  NA6PTrackParCov t = trc;
  if (!Propagator::Instance()->propagateToZ(t, pv.getZ(), mPropOpt)) {
    return false;
  }
  // create a pseudocluster for the vertex, use -1 for layer to distinguish it from detector hits
  NA6PBaseCluster vClu(pv.getX(), pv.getY(), pv.getZ(), 1, -1);
  auto [tx, ty] = t.getSlopesXY(); // project covariance matrix elements using track direction
  auto sxx = pv.getSigmaX2() + (tx * tx * pv.getSigmaZ2()) - (2.0f * tx * pv.getSigmaXZ());
  auto syy = pv.getSigmaY2() + (ty * ty * pv.getSigmaZ2()) - (2.0f * ty * pv.getSigmaYZ());
  auto sxy = pv.getSigmaXY() + (tx * ty * pv.getSigmaZ2()) - (tx * pv.getSigmaYZ()) - (ty * pv.getSigmaXZ());
  vClu.setErr(sxx, sxy, syy);
  auto chi2 = t.getPredictedChi2(vClu);
  mChi2Buffer.push_back(chi2);
  if ((mMaxChi2VertConstr > 0.f && chi2 > mMaxChi2VertConstr) || !t.update(vClu)) {
    trc.getVertexConstrainedParam().invalidate();
    return false;
  }
  trc.setVertexConstrainedParam(t);
  return true;
}

int NA6PFastTrackFitter::getLayersForSeed(int dir, std::array<int, 3>& layForSeed)
{
  // define layers to be used in seed calculation depending on the options, in the order of track fit direction (dir>0 : forward)
  int nLayWithClus = countLayerWithClusters();
  if (nLayWithClus < 2) {
    LOGP(error, "Only {} layers with clusters, at least 2 needed for seeding", nLayWithClus);
    return 0;
  }
  layForSeed[2] = -1; // in case only 2 layers available
  if (nLayWithClus <= 3) {
    for (int i = 0; i < mNClusters; i++) {
      layForSeed[i] = mLayersWithClusters[i];
    }
    if (dir < 0) {
      std::swap(layForSeed[0], layForSeed[nLayWithClus - 1]);
    }
    return nLayWithClus;
  }
  // select layers for seed determination
  if (mSeedOption == kEdgeClusters) {
    if (dir > 0) {
      for (int i = 0; i < 3 && i < nLayWithClus; ++i) {
        layForSeed[i] = mLayersWithClusters[i];
      }
    } else {
      for (int i = 0; i < 3 && i < nLayWithClus; ++i) {
        layForSeed[i] = mLayersWithClusters[nLayWithClus - 1 - i];
      }
    }
  } else { // max lever arm
    layForSeed[0] = mLayersWithClusters[0];
    layForSeed[2] = mLayersWithClusters[nLayWithClus - 1];
    int mid = (layForSeed[0] + layForSeed[2]) / 2, best = -1, bestDist = 999; // Find closest layer to midpoint
    for (int i = 1; i < nLayWithClus - 1; ++i) {
      int dist = std::abs(mLayersWithClusters[i] - mid);
      if (dist < bestDist) {
        bestDist = dist;
        best = mLayersWithClusters[i];
      }
    }
    layForSeed[1] = best;
    if (dir < 0) {
      std::swap(layForSeed[0], layForSeed[2]);
    }
  }
  return 3;
}

bool NA6PFastTrackFitter::computeSeed(int dir, std::array<int, 3>& layForSeed, NA6PTrackPar* seed)
{
  // layers are provided in increasing order if dir>0 (forward fit) and decreasing if dir<0 (backward fit), seed is defined at 1st provided cluster position
  if (!seed) {
    seed = &mSeed;
  }
  int nSeedLayers = layForSeed[2] >= 0 ? 3 : 2;
  const auto &xyz0 = mClusters[layForSeed[0]]->getXYZ(), xyz1 = mClusters[layForSeed[1]]->getXYZ();
  if (nSeedLayers == 2) {                              // straight line seed only
    const auto dpos01 = NA6PLine::getDiff(xyz0, xyz1); // c1 - c0
    seed->setXYZ(xyz0);
    auto pxzI = (dir > 0 ? 1.f : -1.f) / std::hypot(dpos01[0], dpos01[2]);
    seed->setTx(dpos01[0] * pxzI);
    seed->setTy(dpos01[1] * pxzI);
    seed->setQ2Pxz(1.f / mMostProbableP);
    return true;
  }
  const auto& xyz2 = mClusters[layForSeed[2]]->getXYZ();

#ifdef _LOG_MC_TRUTH_MODE_
  auto mcLab = getMCTruthStatus(); // RSREM
  LOGP(info, "seed Dir{} L:{}/{}/{} MCComposedLabel:{}", dir, layForSeed[0], layForSeed[1], layForSeed[2], mcLab.asString());
#endif //_LOG_MC_TRUTH_MODE_
  return Seeder::create(*seed, dir > 0, xyz0, xyz1, xyz2, mSeedImprovePrec, mPID);
}

bool NA6PFastTrackFitter::computeSeed(int dir, NA6PTrackPar* seed)
{
  // compute track seed from 3 (or 2) clusters
  // dir = 1 -> forward dir = -1 -> backward
  if (!seed) {
    seed = &mSeed;
  }
  std::array<int, 3> layForSeed = {-1, -1, -1};
  int nSeedLayers = getLayersForSeed(dir, layForSeed);
  if (nSeedLayers < 2) {
    seed->invalidate();
    LOGP(error, "Cannot compute seed with {} seeding layers", nSeedLayers);
    return false;
  }
  return computeSeed(dir, layForSeed, seed);
}

void NA6PFastTrackFitter::printSeed() const
{
  if (mSeed.isValid()) {
    LOGP(info, "Seed for tracking: {}", mSeed.asString());
  } else
    LOGP(info, "Seed not set");
}

float NA6PFastTrackFitter::fitSeed(NA6PTrackParCov& seed, bool resetCovMat, int dir, bool useLinRef)
{
  // fit track (already seeded but not necessarily at the position of the 1st cluster to fit.
  // If the covariance matrix reset it requested, then the track position is imposed at the 1st cluster to fit (w/o propagating other track params)
  // and the covariance is reset. Otherwise (assuming that the fit is a continuation of the existing fit), the track is propagated to the 1st point
  // (including the cov. matrix propagation).
  // If useLinRef is true, the clone of the seed is used for the KF linearization, otherwise the linearization is done wrt the track being updated.
  // The propagation is done with the PID of the provided seed (or that of linRef, if provided).
  //
  // In case of the sucessful fit return total chi2, otherwise negative number
  if (mNClusters < 2) {
    LOGP(error, "Cannot fit track with only {} clusters\n", mNClusters);
    return -1.f;
  }
  float chi2Tot = 0.f;
  mChi2Buffer.clear();
  auto prop = Propagator::Instance();
  int step = 1, startL = mMinLayerWithCl, stopL = mMaxLayerWithCl + step;
  if (dir < 0) {
    startL = mMaxLayerWithCl;
    step = -1;
    stopL = mMinLayerWithCl - 1;
  }
  if (resetCovMat) {
    auto cl = getCluster(startL);
    seed.setXYZ(cl->getXYZ());
    seed.resetCovariance(NA6PTrackParCov::MaxErrSelRescale);
  }
  NA6PTrackPar linRefLoc(seed);
  mPropOpt.linRef = useLinRef ? &linRefLoc : nullptr;

#ifdef _LOG_MC_TRUTH_MODE_
  auto mcLab = getMCTruthStatus(); // RSREM
#endif                             //_LOG_MC_TRUTH_MODE_

  int cntCl = 0;
  for (int il = startL; il != stopL; il += step) {
    if (!mClusters[il]) {
      continue;
    }
    const auto& cl = *mClusters[il];
    float chi2 = 0;
    if (!prop->propagateToZ(seed, cl.getZ(), mPropOpt)) {
#ifdef _LOG_MC_TRUTH_MODE_
      LOGP(info, "DBG fitSeed fail: propagate il={} z={} MCComposedLabel={} state={}", il, cl.getZ(), mcLab.asString(), seed.asString());
#endif //_LOG_MC_TRUTH_MODE_
      chi2Tot = -1;
      break;
    }
#ifdef _LOG_MC_TRUTH_MODE_
    if (mcLab.isSet()) {
      LOGP(info, "DBGR Before update {}/{}: {}", cntCl, mNClusters, seed.asString());
    }
#endif //_LOG_MC_TRUTH_MODE_
    chi2 = seed.getPredictedChi2(cl);
    mChi2Buffer.push_back(chi2);
#ifdef _LOG_MC_TRUTH_MODE_
    if (mcLab.isSet()) {
      LOGP(info, "DBGR chi2={} for cluster {}", chi2, cl.asString());
    }
#endif //_LOG_MC_TRUTH_MODE_
    if (mMaxChi2Cl > 0.f && chi2 > mMaxChi2Cl) {
#ifdef _LOG_MC_TRUTH_MODE_
      const auto dx = seed.getX() - cl.getX();
      const auto dy = seed.getY() - cl.getY();
      const auto exx = seed.getSigmaX2() + cl.getSigXX();
      const auto eyx = seed.getSigmaYX() + cl.getSigYX();
      const auto eyy = seed.getSigmaY2() + cl.getSigYY();
      const auto det = exx * eyy - eyx * eyx;
      LOGP(info, "DBG fitSeed fail chi2 at {}/{} il={} chi2={:.3f} max={} TrackMCLabel={} cluLabs={} clXYZ=({:.3e},{:.3e},{:.3e}) state={}", cntCl, mNClusters,
           il, chi2, mMaxChi2Cl, mcLab.asString(), clusterLabelsAsString(cl), cl.getX(), cl.getY(), cl.getZ(), seed.asString());
      LOGP(info, "DBG predChi2 new: frame=labXY il={} TrackMCLabel={} cluLabs={} dx={} dy={} trXY=({},{}) clXY=({},{}) trCov=({:.3e},{:.3e},{:.3e}) clCov=({:.3e},{:.3e},{:.3e}) sumcov=({:.3e},{:.3e},{:.3e}) det={} chi2={}",
           il, mcLab.asString(), clusterLabelsAsString(cl), dx, dy, seed.getX(), seed.getY(), cl.getX(), cl.getY(), seed.getSigmaX2(), seed.getSigmaYX(), seed.getSigmaY2(),
           cl.getSigXX(), cl.getSigYX(), cl.getSigYY(), exx, eyx, eyy, det, chi2);
#endif //_LOG_MC_TRUTH_MODE_
      chi2Tot = -1;
      break;
    } else if (!seed.update(cl)) {
#ifdef _LOG_MC_TRUTH_MODE_
      LOGP(info, "DBG fitSeed fail: update il={} chi2={:.3f} TrackMCLabel={} cluLabs={} clXYZ=({:.3e},{:.3e},{:.3e}) state={}",
           il, chi2, mcLab.asString(), clusterLabelsAsString(cl), cl.getX(), cl.getY(), cl.getZ(), seed.asString());
#endif //_LOG_MC_TRUTH_MODE_
      chi2Tot = -1;
      break;
#ifdef _LOG_MC_TRUTH_MODE_
    } else if (mcLab.isSet()) {
      LOGP(info, "DBG fitSeed OK chi2 at {}/{} il={} chi2={:.3f}/{:.3f} max={} TrackMCLabel={} cluLabs={} clXYZ=({:.3e},{:.3e},{:.3e}) state={}", cntCl, mNClusters,
           il, chi2, chi2 + chi2Tot, mMaxChi2Cl, mcLab.asString(), clusterLabelsAsString(cl), cl.getX(), cl.getY(), cl.getZ(), seed.asString());
#endif //_LOG_MC_TRUTH_MODE_
    }
#ifdef _LOG_MC_TRUTH_MODE_
    if (mcLab.isSet()) {
      LOGP(debug, "DBGR After update {}/{}: {}", cntCl, mNClusters, seed.asString());
    }
#endif //_LOG_MC_TRUTH_MODE_
    chi2Tot += chi2;
    cntCl++;
  }

  if (dir < 0 && mPropagateToPrimVert && mIsPrimVertSet && !prop->propagateToZ(seed, mPrimVertZ, mPropOpt)) {
    return -1.f;
  }

  mPropOpt.linRef = nullptr;
  return chi2Tot;
}

void NA6PFastTrackFitter::addClustersToTrack(NA6PTrack& tr)
{
  for (int i = mMinLayerWithCl; i <= mMaxLayerWithCl; i++) {
    const auto* cl = getCluster(i);
    if (cl) {
      tr.addCluster(cl);
    }
  }
}

NA6PMCComposedLabel NA6PFastTrackFitter::getMCTruthStatus()
{
  NA6PMCComposedLabel lbl;
  if (!mCluMCLabels) {
    return lbl;
  }
  std::vector<std::pair<NA6PMCComposedLabel, int>> occurrences;
  for (int jl = mMinLayerWithCl; jl <= mMaxLayerWithCl; ++jl) {
    if (!mClusters[jl]) {
      continue;
    }
    const int cluID = mClusters[jl]->getClusterIndex();
    std::span labels = mCluMCLabels->getLabels(cluID);
    for (const auto currLbl : labels) {
      bool found{false};
      for (auto& occ : occurrences) {
        if (currLbl == occ.first) {
          occ.second++;
          found = true;
          break;
        }
      }
      if (!found) {
        occurrences.emplace_back(currLbl, 1);
      }
    }
  }
  if (occurrences.empty()) {
    return lbl;
  }
  auto occIter = occurrences.begin();
  auto occWin = occIter;
  while (++occIter != occurrences.end()) {
    if (occWin->second < occIter->second) {
      occWin = occIter;
    }
  }
  lbl = occWin->first;
  lbl.setFakeFlag(occWin->second < mNClusters);
  return lbl;
}

std::string NA6PFastTrackFitter::clusterLabelsAsString(const NA6PBaseCluster& cl) const
{
  std::string allLabels = "";
  if (mCluMCLabels) {
    const int cluID = cl.getClusterIndex();
    std::span labels = mCluMCLabels->getLabels(cluID);
    int nLabels = labels.size();
    for (int jLab = 0; jLab < nLabels; jLab++) {
      if (nLabels > 1)
        allLabels.append(fmt::format("({})", jLab));
      NA6PMCComposedLabel lbl = labels[jLab];
      allLabels.append(lbl.asString());
      if (nLabels > 1)
        allLabels.append(",");
    }
  }
  return allLabels;
}
