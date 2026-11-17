// NA6PCCopyright

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TRandom3.h>
#include <fairlogger/Logger.h>
#include "NA6PMuonSpecModularHit.h"
#include "NA6PTrackerCA.h"
#include "NA6PVertexerTracklets.h"
#include "NA6PMuonSpecReconstruction.h"
#include "NA6PLayoutParam.h"

NA6PMuonSpecReconstruction::NA6PMuonSpecReconstruction() : NA6PReconstruction("MuonSpec")
{
  initTracker();
}

NA6PMuonSpecReconstruction::~NA6PMuonSpecReconstruction()
{
  // move here since the NA6PTrackerCA was fwd-declared
}

bool NA6PMuonSpecReconstruction::initTracker()
{
  mMSTracker = std::make_unique<NA6PTrackerCA>();
  mMSTracker->configureFromRecoParamMS();
  const auto& param = NA6PRecoParam::Instance();
  if (param.msDoTrackMSTrackletMIDMatching) {
    mMSTracker->setNLayers(4);
    mMSTracker->setStartLayer(param.vtNLayers);
    for (int j = 0; j < mMSTracker->getNIterations(); ++j) {
      if (mMSTracker->getMinimumNumberOfClusters(j) > 4) {
        mMSTracker->setMinimumNumberOfClusters(j, 4);
      }
    }
    mMSTracker->setDoOutwardPropagation(true);
    const auto& layout = NA6PLayoutParam::Instance();
    const int firstMIDPlane = layout.nMSPlanes - 2;
    mMSTracker->setZForOutwardPropagation(layout.shiftMS[2] + layout.posMSPlaneZ[firstMIDPlane] - layout.thicknessMSPlane[firstMIDPlane] / 2.f);
  }
  createTracksOutput();
  return true;
}

void NA6PMuonSpecReconstruction::createClustersOutput()
{
  useOwnedClusterStorage(); // to write into mClusters
  auto nm = fmt::format("Clusters{}.root", getName());
  mClusFile = TFile::Open(nm.c_str(), "recreate");
  mClusTree = new TTree(fmt::format("clusters{}", getName()).c_str(), fmt::format("{} Clusters", getName()).c_str());
  mClusTree->Branch(getName().c_str(), &hClusPtr);
  mClusTree->Branch(fmt::format("{}MCTruth", getName()).c_str(), &hCluMCLabelsPtr);
  LOGP(info, "Will store {} clusters in {}", getName(), nm);
}

void NA6PMuonSpecReconstruction::writeClusters()
{
  if (mClusTree) {
    mClusTree->Fill();
  }
  LOGP(info, "Saved {} clusters in tree with {} entries", mClusters.size(), mClusTree->GetEntries());
}

void NA6PMuonSpecReconstruction::closeClustersOutput()
{
  if (mClusTree && mClusFile) {
    mClusFile->cd();
    mClusTree->Write();
    delete mClusTree;
    mClusTree = nullptr;
    mClusFile->Close();
    delete mClusFile;
    mClusFile = nullptr;
  }
}

void NA6PMuonSpecReconstruction::setClusters(std::vector<NA6PMuonSpecCluster>& clusters)
{
  hClusPtr = &clusters;
  // Assign a transient index based on position in the vector
  // to be used for storing the cluster indices in the track
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PMuonSpecReconstruction::hitsToRecPoints(const std::vector<NA6PMuonSpecModularHit>& hits, int evID, const std::vector<TParticle>* mcArray)
{
  int nHits = hits.size();
  const auto& layout = NA6PLayoutParam::Instance();

  for (int jHit = 0; jHit < nHits; ++jHit) {
    const auto& hit = hits[jHit];
    double x = hit.getX();
    double y = hit.getY();
    double z = hit.getZ();
    double ex2clu = 5.e-4;
    double ey2clu = 5.e-4;
    if (mCluResX > 0) {
      x = gRandom->Gaus(hit.getX(), mCluResX);
      ex2clu = mCluResX * mCluResX;
    }
    if (mCluResY > 0) {
      y = gRandom->Gaus(hit.getY(), mCluResY);
      ey2clu = mCluResY * mCluResY;
    }
    // very rough cluster size settings
    int clusiz = 1;
    int nDet = hit.getDetectorID();
    // Set layer based on Z position - find closest plane

    float minDist = std::abs(z - layout.posMSPlaneZ[0]);
    int layer = layout.nVerTelPlanes;

    for (int i = 1; i < layout.nMSPlanes; ++i) {
      float dist = std::abs(z - layout.posMSPlaneZ[i]);
      if (dist < minDist) {
        minDist = dist;
        layer = i + layout.nVerTelPlanes;
      }
    }
    // reject electrons produced in the material which are nearby the muon
    bool reject = false;
    if (mDiscardNearbyElectrons && mcArray != nullptr) {
      int mcLabel = hit.getTrackID();
      if (mcLabel >= 0 && mcLabel < static_cast<int>(mcArray->size())) {
        TParticle part = mcArray->at(mcLabel);
        int pdg = part.GetPdgCode();
        if (pdg == 11) {
          // hit from electron
          int mcLabelMoth = part.GetFirstMother();
          if (mcLabelMoth >= 0 && mcLabelMoth < static_cast<int>(mcArray->size())) {
            TParticle mothPart = mcArray->at(mcLabelMoth);
            int pdgMoth = mothPart.GetPdgCode();
            if (pdgMoth != 11 && pdgMoth != 22) {
              // search for hit of the mother particle (unless the mother is a photon or an electron)
              for (int kHit = 0; kHit < nHits; ++kHit) {
                const auto& hitMoth = hits[kHit];
                int candLabel = hitMoth.getTrackID();
                if (candLabel == mcLabelMoth) {
                  double xm = hitMoth.getX();
                  double ym = hitMoth.getY();
                  double zm = hitMoth.getZ();
                  double dist2 = (xm - x) * (xm - x) + (ym - y) * (ym - y);
                  if (dist2 < mMinElectronHitDistance * mMinElectronHitDistance) {
                    reject = true;
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (!reject) {
      int cluID = mClusters.size();
      mClusters.emplace_back(x, y, z, clusiz, layer);
      auto& clu = mClusters.back();
      clu.setErr(ex2clu, 0., ey2clu);
      clu.setDetectorID(nDet);
      clu.setHitID(jHit);
      clu.setClusterIndex(cluID);
      NA6PMCComposedLabel lbl(hit.getTrackID(), evID, 0);
      mCluMCLabels.addElement(cluID, lbl);
    }
  }
}

//____________________________________________________________________________________

void NA6PMuonSpecReconstruction::createTracksOutput()
{
  auto nm = fmt::format("Tracks{}.root", getName());
  mTrackFile = TFile::Open(nm.c_str(), "recreate");
  mTrackTree = new TTree(fmt::format("tracks{}", getName()).c_str(), fmt::format("{} Tracks", getName()).c_str());
  mTrackTree->Branch(getName().c_str(), &hTrackPtr);
  if (mReadMCTruth)
    mTrackTree->Branch(fmt::format("{}MCTruth", getName()).c_str(), &hTrkMCLabelsPtr);
  LOGP(info, "Will store {} tracks in {}", getName(), nm);
}

void NA6PMuonSpecReconstruction::writeTracks()
{
  if (mTrackTree) {
    mTrackTree->Fill();
  }
  LOGP(info, "Saved {} tracks in tree with {} entries", mTracks.size(), mTrackTree->GetEntries());
}

void NA6PMuonSpecReconstruction::closeTracksOutput()
{
  if (mTrackTree && mTrackFile) {
    mTrackFile->cd();
    mTrackTree->Write();
    delete mTrackTree;
    mTrackTree = nullptr;
    mTrackFile->Close();
    delete mTrackFile;
    mTrackFile = nullptr;
  }
}

void NA6PMuonSpecReconstruction::runTracking()
{
  clearTracks();
  if (mReadMCTruth)
    mMSTracker->setClusterMCTruth(hCluMCLabelsPtr);
  const auto& param = NA6PRecoParam::Instance();
  mMSTracker->findTracks(getClusters(), mPrimaryVertex);
  if (param.msDoTrackMSTrackletMIDMatching == false) {
    mTracks = mMSTracker->getTracks();
    for (auto& t : mTracks) {
      t.setStatusMS(NA6PTrack::kFullMSMID);
    }
  } else {
    runMSTrackMIDTrackletMatching();
  }
  if (mReadMCTruth)
    assignMCLabels(mTracks, *hTrkMCLabelsPtr, *hCluMCLabelsPtr);
  writeTracks();
}

void NA6PMuonSpecReconstruction::runMSTrackMIDTrackletMatching()
{
  const auto& param = NA6PRecoParam::Instance();
  const auto& layout = NA6PLayoutParam::Instance();
  auto& clusters = getClusters();
  const auto& trks = mMSTracker->getFinalTracks();

  auto saveTrack = [&trks, this](int j, int flag) -> NA6PTrack& {
    this->mTracks.push_back(trks[j].trackFitFast);
    this->mTracks.back().setStatusMS(flag);
    return this->mTracks.back();
  };

  int nTrks = trks.size();
  int firstMID = layout.nVerTelPlanes + layout.nMSPlanes - 2;
  std::vector<std::pair<NA6PMuonSpecCluster, NA6PMuonSpecCluster>> trkltsMID = mMSTracker->findTracklets(firstMID, firstMID + 1, clusters, mPrimaryVertex);
  int nTrklets = trkltsMID.size();
  NA6PFastTrackFitter* fitter = mMSTracker->getTrackFitter();
  fitter->setMaxChi2Cl(param.msMaxChi2TrClMSTrackletMID);
  fitter->setMaxChi2VertConstr(param.msMaxChi2VertexConstrain);

  // create lookup table of cluster indices
  std::vector<int> clusterLookup(clusters.size(), -1);
  for (size_t i = 0; i < clusters.size(); ++i) {
    int originalID = clusters[i].getClusterIndex();
    if (originalID >= 0 && (size_t)originalID < clusters.size())
      clusterLookup[originalID] = static_cast<int>(i);
  }

  // matching loop
  std::vector<MatchCandidate> candidates;
  std::vector<NA6PTrackParCov> outTrProp; // to reuse propagated outer tracks
  if (nTrklets) {
    for (int jT = 0; jT < nTrks; jT++) {
      const auto& tr = trks[jT].trackFitFast;
      auto& trOutCopy = outTrProp.emplace_back(tr.getOuterParam());
      if (!trOutCopy.isValid()) { // outward propagation is done in the ca tracking
        LOGP(debug, "Track {} has invalid outward propagated parameters, skipping matching", jT);
        continue;
      }
      for (int jS = 0; jS < nTrklets; jS++) {
        const auto& tracklet = trkltsMID[jS];
        NA6PLine trkLine(tracklet.first.getXYZ(), tracklet.second.getXYZ());
        NA6PTrackPar trOutCopyT{trOutCopy};
        if (!Propagator::Instance()->propagateToZ(trOutCopyT, trkLine.mOriginPoint[2], fitter->getPropOpt())) {
          continue;
        }
        auto trackDirCos = trOutCopyT.getDirCos();
        float distXY = std::hypot(trkLine.mOriginPoint[0] - trOutCopyT.getX(), trkLine.mOriginPoint[1] - trOutCopyT.getY());
        float costh = NA6PLine::getScalar(trkLine.mCosinesDirector, trackDirCos);
        if (distXY < param.msMaxDistTrackMSTrackletMID && costh > param.msMinCosThetaTrackMSTrackletMID) {
          costh = std::clamp(costh, -1.f, 1.f);
          float score = distXY / 2. + (1 - costh) / 1.e-4;
          candidates.push_back({jT, jS, score});
        }
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [](auto& a, auto& b) { return a.score < b.score; });

  std::vector<bool> trackUsed(nTrks, false);
  std::vector<bool> segmentUsed(nTrklets, false);
  std::vector<int> trackMatch(nTrks, -1);

  auto saveUnrefittedMatch = [&saveTrack](int trackIndex,
                                          const std::pair<NA6PMuonSpecCluster, NA6PMuonSpecCluster>& tracklet, NA6PTrack::MuonSpecTrackCase caseType) {
    auto& savedTrack = saveTrack(trackIndex, caseType);
    savedTrack.addCluster(&tracklet.first);
    savedTrack.addCluster(&tracklet.second);
  };

  for (auto& c : candidates) {
    if (trackUsed[c.track])
      continue;
    if (segmentUsed[c.segment])
      continue;
    trackMatch[c.track] = c.segment;
    trackUsed[c.track] = true;
    segmentUsed[c.segment] = true;
  }

  for (int jT = 0; jT < nTrks; jT++) {
    int jS = trackMatch[jT];
    if (jS < 0) {
      // track not matched to MID: save with specific status flag
      saveTrack(jT, NA6PTrack::kMSNotMatchedToMID);
      continue;
    }

    const auto& tracklet = trkltsMID[jS];
    if (!param.msDoRefitMSTrackletMID) {
      // Store the MID cluster IDs without including them in the fit.
      saveUnrefittedMatch(jT, tracklet, NA6PTrack::kMSMatchedToMIDnotRefitted);
      continue;
    }

    auto& trOut = outTrProp[jT]; // propagated tracks are guaranteed to exist if there was at least 1 match (or even tracklet)
    fitter->cleanupAndStartFit();
    fitter->addCluster(tracklet.first);
    fitter->addCluster(tracklet.second);
    if (fitter->fitSeedOutward(trOut, false) < 0.) {
      // fit failed: save with specific status flag
      saveUnrefittedMatch(jT, tracklet, NA6PTrack::kMSMatchedToMIDRefitFailed);
      continue;
    }

    // add MS clusters for refit inward
    const auto& trMS = trks[jT].trackFitFast;
    for (int jl = 0; jl < 4; ++jl) {
      int originalID = trMS.getClusterIndex(jl + param.vtNLayers);
      if (originalID >= 0 && (size_t)originalID < clusterLookup.size()) {
        int jNewPos = clusterLookup[originalID];
        if (jNewPos < 0) {
          LOGP(error, "Cluster originalID={} not found for track {} layer {}", originalID, jT, jl);
          continue;
        }
        const auto& cl = clusters[jNewPos];
        fitter->addCluster(cl);
      }
    }
    auto refitInw = trOut;
    float chi2Refit = fitter->fitSeedInward(refitInw, true);
    if (chi2Refit < 0.f) {
      // fit failed: save with specific status flag
      saveUnrefittedMatch(jT, tracklet, NA6PTrack::kMSMatchedToMIDRefitFailed);
      continue;
    }
    saveTrack(jT, NA6PTrack::kMSMatchedToMIDRefitted);
    auto& trMSMID = mTracks.back();
    trMSMID.setInwardParam(refitInw);
    trMSMID.setChi2(chi2Refit);
    trMSMID.setOuterParam(trOut);
    trMSMID.resetClusters();
    fitter->addClustersToTrack(trMSMID);
    if (param.msDoConstrainedTrack) {
      fitter->constrainTrackToVertex(trMSMID, *mPrimaryVertex);
    } else {
      trMSMID.getVertexConstrainedParam().invalidate();
    }
  }
}
