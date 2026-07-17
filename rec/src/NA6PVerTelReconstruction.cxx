// NA6PCCopyright

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <fairlogger/Logger.h>
#include "NA6PVerTelHit.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PVerTelDigit.h"
#include "NA6PVerTelSegmentation.h"
#include "NA6PVerTelDigitizer.h"
#include "NA6PTrackerCA.h"
#include "NA6PVertexerTracklets.h"
#include "NA6PVerTelReconstruction.h"

NA6PVerTelReconstruction::NA6PVerTelReconstruction() : NA6PReconstruction("VerTel")
{
  initAll();
  mClusterizer.initGeometry();
}

NA6PVerTelReconstruction::~NA6PVerTelReconstruction()
{
  // move to cxx since NA6PVertexerTracklets and NA6PTrackerCA were fwd-declared
}

bool NA6PVerTelReconstruction::initAll()
{
  bool retV = initVertexer();
  bool retC = initTracker();
  return (retV & retC);
}

bool NA6PVerTelReconstruction::initVertexer()
{
  if (!mVTTrackletVertexer) {
    mVTTrackletVertexer = std::make_unique<NA6PVertexerTracklets>();
  }
  mVTTrackletVertexer->configureFromRecoParam();
  createVerticesOutput();
  return true;
}

bool NA6PVerTelReconstruction::initTracker()
{
  if (!mVTTracker) {
    mVTTracker = std::make_unique<NA6PTrackerCA>();
  }
  mVTTracker->setNLayers(5);
  mVTTracker->setStartLayer(0);
  mVTTracker->configureFromRecoParamVT();
  createTracksOutput();
  return true;
}

void NA6PVerTelReconstruction::createClustersOutput()
{
  auto nm = fmt::format("Clusters{}.root", getName());
  mClusFile = TFile::Open(nm.c_str(), "recreate");
  mClusTree = new TTree(fmt::format("clusters{}", getName()).c_str(), fmt::format("{} Clusters", getName()).c_str());
  mClusTree->Branch(getName().c_str(), &hClusPtr);
  mClusTree->Branch(fmt::format("{}MCTruth", getName()).c_str(), &hCluMCLabelsPtr);
  LOGP(info, "Will store {} clusters in {}", getName(), nm);
}

void NA6PVerTelReconstruction::writeClusters()
{
  if (mClusTree) {
    mClusTree->Fill();
  }
  LOGP(info, "Saved {} clusters in tree with {} entries", mClusters.size(), mClusTree->GetEntries());
}

void NA6PVerTelReconstruction::closeClustersOutput()
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

void NA6PVerTelReconstruction::setClusters(std::vector<NA6PVerTelCluster>& clusters)
{
  hClusPtr = &clusters;
  // Assign a transient index based on position in the vector
  // to be used for storing the cluster indices in the track
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PVerTelReconstruction::hitsToRecPoints(const std::vector<NA6PVerTelHit>& hits)
{
  int nHits = hits.size();
  NA6PVerTelSegmentation vt;
  if (mEmulateNoGaps)
    vt.setStaggered(true);
  NA6PVerTelDigitizer dig;
  dig.initGeometry();

  for (int jHit = 0; jHit < nHits; ++jHit) {
    const auto& hit = hits[jHit];
    double x = hit.getX();
    double y = hit.getY();
    double z = hit.getZ();
    double xyzLocS[3], xyzLocE[3];
    dig.getHitLocalCoord(hit, xyzLocS, xyzLocE);
    if (vt.isInAcc(xyzLocS[0], xyzLocS[1]) != 1)
      continue;
    double ex2clu = 5.e-4;
    double ey2clu = 5.e-4;
    if (mCluRes > 0) {
      x = gRandom->Gaus(hit.getX(), mCluRes);
      y = gRandom->Gaus(hit.getY(), mCluRes);
      ex2clu = mCluRes * mCluRes;
      ey2clu = mCluRes * mCluRes;
    }
    double eloss = hit.getHitValue();
    // very rough cluster size settings
    int clusiz = 2;
    if (eloss > 2.e-5 && eloss < 5.e-5)
      clusiz = 3;
    else if (eloss > 5.e-5)
      clusiz = 4;
    int nDet = hit.getDetectorID();
    int idPart = hit.getTrackID();
    int layer = nDet / 4;
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
void NA6PVerTelReconstruction::digitsToRecPoints(const std::vector<NA6PVerTelDigit>& vtDigits, const NA6PMCTruthContainer& digMCLabels)
{
  int nDigits = vtDigits.size();
  int nMClabels = digMCLabels.getNElements();
  mClusterizer.process(vtDigits, digMCLabels, mClusters, mCluMCLabels);
  int nCluMClabels = mCluMCLabels.getNElements();
  int nClusters = mClusters.size();
  LOGP(info, " --> nClusters = {} nCluMCLabels = {}", nClusters, nCluMClabels);
}

//____________________________________________________________________________________

void NA6PVerTelReconstruction::createVerticesOutput()
{
  auto nm = fmt::format("Vertices{}.root", getName());
  mVertexFile = TFile::Open(nm.c_str(), "recreate");
  mVertexTree = new TTree(fmt::format("vertices{}", getName()).c_str(), fmt::format("{} Vertices", getName()).c_str());
  mVertexTree->Branch(getName().c_str(), &hVerticesPtr);
  LOGP(info, "Will store {} vertices in {}", getName(), nm);
}

void NA6PVerTelReconstruction::writeVertices()
{
  if (mVertexTree) {
    mVertexTree->Fill();
  }
  LOGP(info, "Saved {} vertices in tree with {} entries", mVertices.size(), mVertexTree->GetEntries());
}

void NA6PVerTelReconstruction::closeVerticesOutput()
{
  if (mVertexTree && mVertexFile) {
    mVertexFile->cd();
    mVertexTree->Write();
    delete mVertexTree;
    mVertexTree = nullptr;
    mVertexFile->Close();
    delete mVertexFile;
    mVertexFile = nullptr;
  }
}

void NA6PVerTelReconstruction::runVertexerTracklets()
{
  if (!mVTTrackletVertexer) {
    LOGP(info, "Tracklet vertexer not initialized, will call initialization");
    initVertexer();
  }
  clearVertices();
  mVTTrackletVertexer->findVertices(*hClusPtr, mVertices);
  writeVertices();
}
//____________________________________________________________________________________

void NA6PVerTelReconstruction::createTracksOutput()
{
  auto nm = fmt::format("Tracks{}.root", getName());
  mTrackFile = TFile::Open(nm.c_str(), "recreate");
  mTrackTree = new TTree(fmt::format("tracks{}", getName()).c_str(), fmt::format("{} Tracks", getName()).c_str());
  mTrackTree->Branch(getName().c_str(), &hTrackPtr);
  LOGP(info, "Will store {} tracks in {}", getName(), nm);
}

void NA6PVerTelReconstruction::writeTracks()
{
  if (mTrackTree) {
    mTrackTree->Fill();
  }
  LOGP(info, "Saved {} tracks in tree with {} entries", mTracks.size(), mTrackTree->GetEntries());
}

void NA6PVerTelReconstruction::closeTracksOutput()
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

void NA6PVerTelReconstruction::runTracking()
{
  clearTracks();
  mVTTracker->findTracks(*hClusPtr, *hCluMCLabelsPtr, mPrimaryVertex);
  mTracks = mVTTracker->getTracks();
  writeTracks();
}
