#ifndef NA6P_EVENT_READER_H
#define NA6P_EVENT_READER_H

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"

#include "NA6PVertex.h"
#include "NA6PTrack.h"
#include "NA6PMatching.h"
#include "NA6PMCEventHeader.h"

#include <cstdint>
#include <memory>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <algorithm>

class NA6PEventReader {
public:
  NA6PEventReader(
    const char* fileNameVerTel         = "VerticesVerTel.root",
    const char* fileNameTracksVerTel   = "TracksVerTel.root",
    const char* fileNameTracksMuonSpec = "TracksMuonSpec.root",
    const char* fileNameTracksMatching = "TracksMatching.root",
    const char* fileNameMC             = "MCKine.root",
    bool readVerTel = true,
    bool readTracksVerTel = true,
    bool readTracksMuonSpec = true,
    bool readTracksMatching = true,
    bool readMC = true
  )
  {
    if (readVerTel) {
      openFile(mFileVerTel, fileNameVerTel);
      setupTree(mFileVerTel.get(), "verticesVerTel", mTreeVerTel);
    }

    if (readTracksVerTel) {
      openFile(mFileTracksVerTel, fileNameTracksVerTel);
      setupTree(mFileTracksVerTel.get(), "tracksVerTel", mTreeTracksVerTel);
    }

    if (readTracksMuonSpec) {
      openFile(mFileTracksMuonSpec, fileNameTracksMuonSpec);
      setupTree(mFileTracksMuonSpec.get(), "tracksMuonSpec", mTreeTracksMuonSpec);
    }

    if (readTracksMatching) {
      openFile(mFileTracksMatching, fileNameTracksMatching);
      setupTree(mFileTracksMatching.get(), "tracksMatching", mTreeTracksMatching);
    }

    if (readMC) {
      openFile(mFileMC, fileNameMC);
      setupTree(mFileMC.get(), "mckine", mTreeMC);
    }

    if (mTreeVerTel) {
      mTreeVerTel->SetBranchAddress("VerTel", &mVerticesVerTel);
    }

    if (mTreeTracksVerTel) {
      mTreeTracksVerTel->SetBranchAddress("VerTel", &mTracksVerTel);
    }

    if (mTreeTracksMuonSpec) {
      mTreeTracksMuonSpec->SetBranchAddress("MuonSpec", &mTracksMuonSpec);
    }

    if (mTreeTracksMatching) {
      mTreeTracksMatching->SetBranchAddress("Matching", &mMatches);
    }

    if (mTreeMC) {
      mTreeMC->SetBranchAddress("header", &mMCHeader);
      mTreeMC->SetBranchAddress("tracks", &mMCParticles);
    }

    determineEntries();
  }

  std::int64_t entries() const { return mEntries; }

  bool loadEvent(std::int64_t i)
  {
    if (i < 0 || i >= mEntries) {
      return false;
    }

    if (mTreeMC) {
      mTreeMC->GetEntry(i);
    }

    if (mTreeVerTel) {
      mTreeVerTel->GetEntry(i);
    }

    if (mTreeTracksVerTel) {
      mTreeTracksVerTel->GetEntry(i);
    }

    if (mTreeTracksMuonSpec) {
      mTreeTracksMuonSpec->GetEntry(i);
    }

    if (mTreeTracksMatching) {
      mTreeTracksMatching->GetEntry(i);
    }

    mCurrentEntry = i;
    return true;
  }

  std::int64_t currentEntry() const { return mCurrentEntry; }

  NA6PMCEventHeader* mcHeader() const { return mMCHeader; }

  const std::vector<NA6PVertex>& verticesVerTel() const
  {
    static const std::vector<NA6PVertex> empty;
    return mVerticesVerTel ? *mVerticesVerTel : empty;
  }

  const std::vector<NA6PTrack>& tracksVerTel() const
  {
    static const std::vector<NA6PTrack> empty;
    return mTracksVerTel ? *mTracksVerTel : empty;
  }

  const std::vector<NA6PTrack>& tracksMuonSpec() const
  {
    static const std::vector<NA6PTrack> empty;
    return mTracksMuonSpec ? *mTracksMuonSpec : empty;
  }

  const std::vector<NA6PMatch>& matches() const
  {
    static const std::vector<NA6PMatch> empty;
    return mMatches ? *mMatches : empty;
  }

  const std::vector<TParticle>& mcParticles() const
  {
    static const std::vector<TParticle> empty;
    return mMCParticles ? *mMCParticles : empty;
  }

  const NA6PTrack* getVTTrack(size_t i) const
  {
    const auto& tr = tracksVerTel();
    return i < tr.size() ? &tr[i] : nullptr;
  }

  const NA6PTrack* getMSTrack(size_t i) const
  {
    const auto& tr = tracksMuonSpec();
    return i < tr.size() ? &tr[i] : nullptr;
  }

  const TParticle* getMCParticle(int i) const
  {
    const auto& mc = mcParticles();
    return (i >= 0 && static_cast<size_t>(i) < mc.size()) ? &mc[i] : nullptr;
  }

private:
  static void openFile(std::unique_ptr<TFile>& file, const char* name)
  {
    file.reset(TFile::Open(name, "READ"));

    if (!file || file->IsZombie()) {
      throw std::runtime_error(std::string("Cannot open file: ") + name);
    }
  }

  static void setupTree(TFile* file, const char* treeName, TTree*& tree)
  {
    tree = nullptr;

    if (!file) {
      return;
    }

    file->GetObject(treeName, tree);

    if (!tree) {
      throw std::runtime_error(std::string("Cannot find tree: ") + treeName);
    }
  }

  void determineEntries()
  {
    std::vector<std::int64_t> n;

    if (mTreeMC) {
      n.push_back(mTreeMC->GetEntries());
    }

    if (mTreeVerTel) {
      n.push_back(mTreeVerTel->GetEntries());
    }

    if (mTreeTracksVerTel) {
      n.push_back(mTreeTracksVerTel->GetEntries());
    }

    if (mTreeTracksMuonSpec) {
      n.push_back(mTreeTracksMuonSpec->GetEntries());
    }

    if (mTreeTracksMatching) {
      n.push_back(mTreeTracksMatching->GetEntries());
    }

    if (n.empty()) {
      throw std::runtime_error("No input trees were loaded");
    }

    const auto [minIt, maxIt] = std::minmax_element(n.begin(), n.end());

    if (*minIt != *maxIt) {
      std::cerr << "WARNING: input trees have different number of entries. "
                << "Using min entries = " << *minIt
                << ", max entries = " << *maxIt
                << std::endl;
    }

    mEntries = *minIt;
  }

  bool hasVerticesVerTel() const { return mTreeVerTel != nullptr; }
  bool hasTracksVerTel() const { return mTreeTracksVerTel != nullptr; }
  bool hasTracksMuonSpec() const { return mTreeTracksMuonSpec != nullptr; }
  bool hasMatching() const { return mTreeTracksMatching != nullptr; }
  bool hasMC() const { return mTreeMC != nullptr; }

private:
  std::unique_ptr<TFile> mFileVerTel;
  std::unique_ptr<TFile> mFileTracksVerTel;
  std::unique_ptr<TFile> mFileTracksMuonSpec;
  std::unique_ptr<TFile> mFileTracksMatching;
  std::unique_ptr<TFile> mFileMC;

  TTree* mTreeVerTel = nullptr;
  TTree* mTreeTracksVerTel = nullptr;
  TTree* mTreeTracksMuonSpec = nullptr;
  TTree* mTreeTracksMatching = nullptr;
  TTree* mTreeMC = nullptr;

  std::vector<NA6PVertex>* mVerticesVerTel = nullptr;
  std::vector<NA6PTrack>* mTracksVerTel = nullptr;
  std::vector<NA6PTrack>* mTracksMuonSpec = nullptr;
  std::vector<NA6PMatch>* mMatches = nullptr;
  std::vector<TParticle>* mMCParticles = nullptr;
  NA6PMCEventHeader* mMCHeader = nullptr;

  std::int64_t mEntries = 0;
  std::int64_t mCurrentEntry = -1;
};

#endif
