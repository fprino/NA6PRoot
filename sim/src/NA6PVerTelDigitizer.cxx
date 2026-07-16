// NA6PCCopyright

#include <ranges>
#include <numeric>
#include <fairlogger/Logger.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TTree.h>

#include "NA6PLayoutParam.h"
#include "NA6PVerTelDigitizer.h"

void NA6PVerTelDigitizer::init(const char* filename, const char* geoname)
{
  const auto& param = NA6PLayoutParam::Instance();
  mNumberOfModules = param.nVerTelPlanes * NA6PGeometryManager::kNVTModulesPerLayer;
  mModules.resize(mNumberOfModules);
  mThresholds.assign(mNumberOfModules * NA6PVerTelSegmentation::NTilesPerModule, kDefaultThresholdEl);
  initGeometry(filename, geoname);
  createDigitsOutput();
}

void NA6PVerTelDigitizer::initGeometry(const char* filename, const char* geoname)
{
  if (!mGeoManager.loadGeometry(filename, geoname)) {
    LOGP(fatal, "Load of geometry not successful");
  }
}

void NA6PVerTelDigitizer::createDigitsOutput()
{
  auto nm = fmt::format("Digits{}.root", getName());
  mDigitsFile = TFile::Open(nm.c_str(), "recreate");
  mDigitsTree = new TTree(fmt::format("digits{}", getName()).c_str(), fmt::format("{} Digits", getName()).c_str());
  mDigitsTree->Branch(getName().c_str(), &hDigitsPtr);
  mDigitsTree->Branch(fmt::format("{}MCTruth", getName()).c_str(), &hMCLabelsPtr);
  LOGP(info, "Will store {} digits in {}", getName(), nm);
}

void NA6PVerTelDigitizer::closeDigitsOutput()
{
  if (mDigitsTree && mDigitsFile) {
    mDigitsFile->cd();
    mDigitsTree->Write();
    delete mDigitsTree;
    mDigitsTree = nullptr;
    mDigitsFile->Close();
    delete mDigitsFile;
    mDigitsFile = nullptr;
  }
}

void NA6PVerTelDigitizer::writeDigits()
{
  if (mDigitsTree) {
    mDigitsTree->Fill();
  }
  LOGP(info, "Saved {} digits in tree with {} entries", mDigits.size(), mDigitsTree->GetEntries());
}

void NA6PVerTelDigitizer::process(const std::vector<NA6PVerTelHit>& hits, int layer)
{
  clearDigits();
  int nHits = hits.size();
  std::vector<int> hitIdx(nHits);
  std::iota(std::begin(hitIdx), std::end(hitIdx), 0);
  // sort hits to improve memory access
  std::sort(hitIdx.begin(), hitIdx.end(), [&hits](auto lhs, auto rhs) {
    return hits[lhs].getDetectorID() < hits[rhs].getDetectorID();
  });
  for (int i : hitIdx | std::views::filter([&](int idx) {
                 if (layer < 0)
                   return true;
                 return detID2Layer(hits[idx].getDetectorID()) == layer;
               })) {
    processHit(hits[i]);
  }
  finalizeDigits();
  writeDigits();
}

void NA6PVerTelDigitizer::getHitLocalCoord(const NA6PVerTelHit& hit, double xyzLocS[3], double xyzLocE[3])
{
  auto modID = hit.getDetectorID();
  auto& matrix = mGeoManager.getMatrix(modID);

  double xyzGloS[3] = {hit.getXIn(), hit.getYIn(), hit.getZIn()};
  double xyzGloE[3] = {hit.getXOut(), hit.getYOut(), hit.getZOut()};

  matrix.MasterToLocal(xyzGloS, xyzLocS);
  xyzLocS[0] += mGeoManager.getModuleHalfX(modID);
  xyzLocS[1] += mGeoManager.getModuleHalfY(modID);
  matrix.MasterToLocal(xyzGloE, xyzLocE);
  xyzLocE[0] += mGeoManager.getModuleHalfX(modID);
  xyzLocE[1] += mGeoManager.getModuleHalfY(modID);
}

void NA6PVerTelDigitizer::processHit(const NA6PVerTelHit& hit)
{
  auto modID = hit.getDetectorID();
  auto& mod = mModules[modID];
  mSegmentation.setDetectorID(modID);
  if (mod.isDisabled()) {
    LOGP(info, "Skipping disabled module {}", modID);
    return;
  }
  double xyzLocS[3], xyzLocE[3];
  getHitLocalCoord(hit, xyzLocS, xyzLocE);
  // we assume for the time being that the charge is deposited uniformly along all the 50 um Si thickness
  // -> to be improved once we have a better knowledge of charge collection in MOSAIX
  double deltaX = xyzLocE[0] - xyzLocS[0];
  double deltaY = xyzLocE[1] - xyzLocS[1];
  float pitchX = NA6PVerTelSegmentation::ActiveDX / NA6PVerTelSegmentation::NColsPerTile;
  float pitchY = NA6PVerTelSegmentation::ActiveDYTile / NA6PVerTelSegmentation::NRowsPerTile;
  int nSteps = std::max(std::abs(deltaX) / pitchX, std::abs(deltaY) / pitchY);
  if (nSteps < 1)
    nSteps = 1;
  if (nSteps > 100)
    nSteps = 100;
  float chargePerStep = hit.getHitValue() * kGeVToEl / nSteps;
  float stepX = deltaX / nSteps;
  float stepY = deltaY / nSteps;
  float x = static_cast<float>(xyzLocS[0]) + 0.5f * stepX;
  float y = static_cast<float>(xyzLocS[1]) + 0.5f * stepY;
  for (int iStep = 0; iStep < nSteps; ++iStep) {
    UShort_t rsu, tile, row, col;
    bool digOk = mSegmentation.localToIndices(x, y, rsu, tile, row, col);
    if (digOk) {
      auto key = mod.getOrderingKey(rsu, tile, row, col);
      NA6PMCComposedLabel lbl(hit.getTrackID(), 0, 0);
      PreDigit* pd = mod.findDigit(key);
      if (!pd) {
        mod.addDigit(key, rsu, tile, row, col, chargePerStep, lbl);
      } else { // there is already a digit at this slot, account as PreDigitExtra contribution
        pd->addContribution(chargePerStep, lbl);
      }
    }
    x += stepX;
    y += stepY;
  }
}

void NA6PVerTelDigitizer::finalizeDigits()
{
  for (int jMod = 0; jMod < mNumberOfModules; ++jMod) {
    auto& mod = mModules[jMod];
    if (mod.isDisabled()) {
      LOGP(info, "Skipping disabled module {}", jMod);
      continue;
    }
    auto& buffer = mod.getPreDigits();
    if (buffer.empty()) {
      continue;
    }
    auto itBeg = buffer.begin();
    auto iter = itBeg;
    for (; iter != buffer.end(); ++iter) {
      auto& preDig = iter->second;
      int jTile = NA6PVerTelSegmentation::getTileId(jMod, preDig.pixID.getRsu(), preDig.pixID.getTile());
      if (jTile >= 0 && preDig.charge >= mThresholds[jTile]) {
        preDig.sortLabelsByEnergy();
        int digID = mDigits.size();
        mDigits.emplace_back(static_cast<uint16_t>(jMod), preDig.pixID);
        for (auto& plab : preDig.labels) {
          // here we can add selections on the labels to be stored in case of multiple labels per digit
          // e.g. remove delta electrons, remove particles contributing with much smaller energy than the others, etc.
          mMCLabels.addElement(digID, plab.second);
        }
      }
    }
    buffer.erase(itBeg, iter); // erase processed entries; iter == end() for now
  }
}
