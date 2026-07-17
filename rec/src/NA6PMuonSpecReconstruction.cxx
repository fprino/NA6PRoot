// NA6PCCopyright

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <fairlogger/Logger.h>
#include "NA6PMuonSpecModularHit.h"
#include "NA6PTrackerCA.h"
#include "NA6PVertexerTracklets.h"
#include "NA6PMuonSpecReconstruction.h"
#include "NA6PLayoutParam.h"

ClassImp(NA6PMuonSpecReconstruction)

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
  createTracksOutput();
  return true;
}

void NA6PMuonSpecReconstruction::createClustersOutput()
{
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

void NA6PMuonSpecReconstruction::hitsToRecPoints(const std::vector<NA6PMuonSpecModularHit>& hits)
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
    int idPart = hit.getTrackID();
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
    int cluID = mClusters.size();
    mClusters.emplace_back(x, y, z, clusiz, layer);
    auto& clu = mClusters.back();
    clu.setErr(ex2clu, 0., ey2clu);
    clu.setDetectorID(nDet);
    clu.setParticleID(idPart);
    clu.setHitID(jHit);
    NA6PMCComposedLabel lbl(hit.getTrackID(), 0, 0);
    mCluMCLabels.addElement(cluID, lbl);
  }
}

//____________________________________________________________________________________

void NA6PMuonSpecReconstruction::createTracksOutput()
{
  auto nm = fmt::format("Tracks{}.root", getName());
  mTrackFile = TFile::Open(nm.c_str(), "recreate");
  mTrackTree = new TTree(fmt::format("tracks{}", getName()).c_str(), fmt::format("{} Tracks", getName()).c_str());
  mTrackTree->Branch(getName().c_str(), &hTrackPtr);
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
  mMSTracker->findTracks(*hClusPtr, *hCluMCLabelsPtr, mPrimaryVertex);
  mTracks = mMSTracker->getTracks();
  writeTracks();
}
