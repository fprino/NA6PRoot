#include "NA6PMCTruthContainer.h"
#include "NA6PVerTelHit.h"
#include "NA6PVerTelDigit.h"
#include "NA6PVerTelClusterizer.h"
#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>

std::vector<NA6PVerTelCluster> mClusters, *hClusPtr = &mClusters; // vector of clusters
NA6PMCTruthContainer mMCLabels, *hMCLabelsPtr = &mMCLabels;       // MC labels
TFile* mClusFile = nullptr;                                       // file with clusters
TTree* mClusTree = nullptr;                                       // tree of clusters

std::string getName()
{
  return "VerTel";
}

void createClustersOutput()
{
  auto nm = fmt::format("Clusters{}.root", getName());
  mClusFile = TFile::Open(nm.c_str(), "recreate");
  mClusTree = new TTree(fmt::format("clusters{}", getName()).c_str(), fmt::format("{} Clusters", getName()).c_str());
  mClusTree->Branch(getName().c_str(), &hClusPtr);
  mClusTree->Branch(fmt::format("{}MCTruth", getName()).c_str(), &hMCLabelsPtr);
  LOGP(info, "Will store {} clusters in {}", getName(), nm);
}

void clearClusters()
{
  mClusters.clear();
  mMCLabels.clear_andfreememory();
}

void writeClusters()
{
  if (mClusTree) {
    mClusTree->Fill();
  }
  LOGP(info, "Saved {} clusters in tree with {} entries", mClusters.size(), mClusTree->GetEntries());
}

void closeClustersOutput()
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

void digitsToClusters(const char* dirSimu = ".")
{

  TFile* fd = TFile::Open(Form("%s/DigitsVerTel.root", dirSimu));
  TTree* td = (TTree*)fd->Get("digitsVerTel");
  std::vector<NA6PVerTelDigit> vtDigits, *vtDigitsPtr = &vtDigits;
  NA6PMCTruthContainer digMCLabels, *digMCLabelsPtr = &digMCLabels;
  td->SetBranchAddress("VerTel", &vtDigitsPtr);
  td->SetBranchAddress("VerTelMCTruth", &digMCLabelsPtr);
  int nEv = td->GetEntries();
  NA6PVerTelClusterizer cl;
  cl.initGeometry();
  createClustersOutput();

  for (int jEv = 0; jEv < nEv; jEv++) {
    td->GetEvent(jEv);
    int nDigits = vtDigits.size();
    int nMClabels = digMCLabels.getNElements();
    LOGP(info, "Event {} nDigits = {} nDigMClabels = {}", jEv, nDigits, nMClabels);
    clearClusters();
    cl.process(vtDigits, digMCLabels, mClusters, mMCLabels);
    int nCluMClabels = mMCLabels.getNElements();
    int nClusters = mClusters.size();
    LOGP(info, " --> nClusters = {} nCluMCLabels = {}", nClusters, nCluMClabels);
    writeClusters();
  }
  closeClustersOutput();
  fd->Close();
}
