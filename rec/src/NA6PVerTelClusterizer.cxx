// NA6PCCopyright

#include "NA6PLayoutParam.h"
#include "NA6PGeometryManager.h"
#include "NA6PVerTelClusterizer.h"

NA6PVerTelClusterizer::NA6PVerTelClusterizer()
{
  mChipDigits.reserve(10000);
  int nModules = NA6PLayoutParam::Instance().nVerTelPlanes * NA6PGeometryManager::kNVTModulesPerLayer;
  int totChips = nModules * nTilesPerModule;
  mChipSlices.reserve(totChips);
}

void NA6PVerTelClusterizer::initGeometry(const char* filename, const char* geoname)
{
  if (!mGeoManager.loadGeometry(filename, geoname)) {
    LOGP(fatal, "Load of geometry not successful");
  }
}

void NA6PVerTelClusterizer::process(const std::vector<NA6PVerTelDigit>& vtDigits,
                                    const NA6PMCTruthContainer& mcDigLabels,
                                    std::vector<NA6PVerTelCluster>& outClusters,
                                    NA6PMCTruthContainer& mcCluLabels)
{
  // read digits and store them in the vector mChipDigits in the compact format of DigitProxy
  mDigitsPtr = &vtDigits;
  mDigMCLabelsPtr = &mcDigLabels;
  mClustersPtr = &outClusters;
  mCluMCLabelsPtr = &mcCluLabels;
  int nDigits = vtDigits.size();
  mChipDigits.clear();
  mChipSlices.clear();
  for (int jDig = 0; jDig < nDigits; ++jDig) {
    const auto& dig = vtDigits.at(jDig);
    uint16_t mod = dig.getDetectorID();
    uint32_t rsu = dig.getRSU();
    uint32_t tile = dig.getTile();
    uint32_t row = dig.getRow();
    uint32_t col = dig.getCol();
    int chipID = tile + nTilesPerRSU * rsu + nTilesPerModule * mod;
    mChipDigits.push_back({jDig, chipID, col, row});
  }
  // sort mChipDigits to allow grouping by chip
  std::sort(mChipDigits.begin(), mChipDigits.end(), [](const DigitProxy& a, const DigitProxy& b) {
    if (a.chipID != b.chipID)
      return a.chipID < b.chipID;
    if (a.col != b.col)
      return a.col < b.col;
    return a.row < b.row;
  });
  // fill mChipSlices with the first and last index of each chip
  if (!mChipDigits.empty()) {
    size_t start = 0;
    int currentChip = mChipDigits[0].chipID;
    for (size_t i = 1; i < mChipDigits.size(); ++i) {
      if (mChipDigits[i].chipID != currentChip) {
        mChipSlices.push_back({currentChip, start, i - 1});
        start = i;
        currentChip = mChipDigits[i].chipID;
      }
    }
    mChipSlices.push_back({currentChip, start, mChipDigits.size() - 1});
  }
  // clusterize chip-by-chip
  for (size_t i = 0; i < mChipSlices.size(); ++i) {
    const auto& slice = mChipSlices[i];
    processChip(slice.startIdx, slice.endIdx);
  }
}

void NA6PVerTelClusterizer::processChip(size_t first,
                                        size_t last)
{
  initChip(first);
  for (size_t j = first + 1; j <= last; ++j) {
    updateChip(j);
  }
  finishChip();
}

void NA6PVerTelClusterizer::initChip(uint32_t jPix)
{
  // init chip with the 1st unmasked pixel (entry "from" in the mChipData)
  mColumn1.fill(-1);
  mColumn2.fill(-1);
  mPrevColDataPtr = mColumn1.data() + 1;
  mCurrColDataPtr = mColumn2.data() + 1;
  mNoLeftCol = true;
  const auto& pix = mChipDigits[jPix];
  mCurrCol = pix.col;
  // this is the first cluster on this chip, so it takes index 0
  mCurrColDataPtr[pix.row] = 0;
  mPixels.clear();
  mPixels.emplace_back(-1, jPix); // id of current pixel
  mPreClusters.clear();
  mPreClusters.emplace_back();
}

void NA6PVerTelClusterizer::updateChip(uint32_t jPix)
{
  const auto& pix = mChipDigits[jPix];
  int row = pix.row;         // should be int because then we compute row-1, which can be negative
  if (mCurrCol != pix.col) { // switch the buffers
    swapColumnBuffers();
    resetColumn(mCurrColDataPtr);
    mNoLeftCol = false;
    if (pix.col > mCurrCol + 1) {
      // no connection with previous column, this pixel cannot belong to any of the
      // existing preclusters, create a new precluster and flag to check only the row above for next pixels of this column
      mCurrCol = pix.col;
      addNewPrecluster(jPix, row);
      mNoLeftCol = true;
      return;
    }
    mCurrCol = pix.col;
  }

  if (mNoLeftCol) { // check only the row above
    if (mCurrColDataPtr[row - 1] >= 0) {
      expandPreCluster(jPix, row, mCurrColDataPtr[row - 1]); // attach to the precluster of the previous row
    } else {
      addNewPrecluster(jPix, row); // start new precluster
    }
  } else {
    // row above should be always checked
    int nNeighb = 0, lowestIndex = mCurrColDataPtr[row - 1];
    int *nbrCol[4], nbrRow[4];
    if (lowestIndex >= 0) {
      nbrCol[nNeighb] = mCurrColDataPtr;
      nbrRow[nNeighb++] = row - 1;
    } else {
      lowestIndex = 0x7ffff;
    }
    // loop to allow for diagonal clusters
    for (int i : {-1, 0, 1}) {
      auto v = mPrevColDataPtr[row + i];
      if (v >= 0) {
        nbrCol[nNeighb] = mPrevColDataPtr;
        nbrRow[nNeighb] = row + i;
        if (v < lowestIndex) {
          lowestIndex = v;
        }
        nNeighb++;
      }
    }
    if (!nNeighb) {
      // no neighbours, create new precluster
      addNewPrecluster(jPix, row);
    } else {
      // attach to the adjacent precluster with smallest index
      expandPreCluster(jPix, row, lowestIndex);
      if (nNeighb > 1) {
        for (int jNeighb = 0; jNeighb < nNeighb; jNeighb++) {
          // reassign precluster index to smallest one, replicating updated values to columns caches
          int neighborClusterID = (nbrCol[jNeighb])[nbrRow[jNeighb]];
          mPreClusters[neighborClusterID].index = lowestIndex;
          (nbrCol[jNeighb])[nbrRow[jNeighb]] = lowestIndex;
        }
      }
    }
  }
}

void NA6PVerTelClusterizer::finishChip()
{
  int nPreclusters = mPreClusters.size();
  // account for the eventual reindexing of preClusters: Id2 might have been reindexed to Id1, which later was reindexed to Id0
  for (int i = 1; i < nPreclusters; i++) {
    if (mPreClusters[i].index != i) { // reindexing is always done towards smallest index
      mPreClusters[i].index = mPreClusters[mPreClusters[i].index].index;
    }
  }
  for (int i1 = 0; i1 < nPreclusters; ++i1) {
    auto& preCluster = mPreClusters[i1];
    auto ci = preCluster.index;
    if (ci < 0) {
      continue;
    }
    BBox bbox;
    int nlab = 0;
    int next = preCluster.head;
    mPixArrBuff.clear();
    while (next >= 0) {
      const auto& pixEntry = mPixels[next];
      const auto pix = mChipDigits[pixEntry.second];
      mPixArrBuff.push_back(pix); // needed for cluster topology
      bbox.adjust(pix.row, pix.col);
      next = pixEntry.first;
    }
    preCluster.index = -1;
    for (int i2 = i1 + 1; i2 < nPreclusters; ++i2) {
      auto& preCluster2 = mPreClusters[i2];
      if (preCluster2.index != ci) {
        continue;
      }
      next = preCluster2.head;
      while (next >= 0) {
        const auto& pixEntry = mPixels[next];
        const auto pix = mChipDigits[pixEntry.second];
        mPixArrBuff.push_back(pix); // needed for cluster topology
        bbox.adjust(pix.row, pix.col);
        next = pixEntry.first;
      }
      preCluster2.index = -1;
    }
    if (bbox.isAcceptableSize()) {
      if (pixelsToRecPoint()) {
        // save cluster
        // not doable now because of parallel open Pull Request
      }
    } else {
      // special treatment of large clusters if needed
    }
  }
}

bool NA6PVerTelClusterizer::pixelsToRecPoint()
{
  float aveCol = 0.f, aveRow = 0.f;
  size_t clusiz = mPixArrBuff.size();
  int refMod = -999;
  int refRsu = -999;
  int refTile = -999;
  for (const auto& pix : mPixArrBuff) {
    int jDig = pix.originalIndex;
    const auto& dig = mDigitsPtr->at(jDig);
    uint16_t mod = dig.getDetectorID();
    uint32_t rsu = dig.getRSU();
    uint32_t tile = dig.getTile();
    uint32_t row = dig.getRow();
    uint32_t col = dig.getCol();
    aveCol += col;
    aveRow += row;
    if (refMod < 0)
      refMod = mod;
    if (mod != refMod) {
      LOGP(error, "Mismatch in modID {} vs. {}", mod, refMod);
      return false;
    }
    if (refRsu < 0)
      refRsu = rsu;
    if (rsu != refRsu) {
      LOGP(error, "Mismatch in RSU {} vs. {}", rsu, refRsu);
      return false;
    }
    if (refTile < 0)
      refTile = tile;
    if (tile != refTile) {
      LOGP(error, "Mismatch in Tile {} vs. {}", tile, refTile);
      return false;
    }
  }
  aveCol /= clusiz;
  aveRow /= clusiz;
  float xLocClu, yLocClu;
  UShort_t row0 = static_cast<UShort_t>(std::floor(aveRow));
  UShort_t col0 = static_cast<UShort_t>(std::floor(aveCol));
  bool cluOk = mSegmentation.indicesToLocal(refRsu, refTile, row0, col0, xLocClu, yLocClu);
  xLocClu += (aveCol - col0) * pitchX;
  yLocClu += (aveRow - row0) * pitchY;
  double xyzGloClu[3];
  getClusterGlobalCoord(refMod, xLocClu, yLocClu, xyzGloClu);
  int layer = refMod / 4;
  int cluID = mClustersPtr->size();

  mClustersPtr->emplace_back(xyzGloClu[0], xyzGloClu[1], xyzGloClu[2], clusiz, layer);
  auto& clu = mClustersPtr->back();
  clu.setDetectorID(refMod);
  // Very preliminary definition of cluster uncertainties
  // -> use pitch / sqrt(12) independently of cluster shape / size
  //    to avoid penalizing large clusters
  //    better tuning as a function of cluster shape should be done based
  //    on MOSAIX characterization measurements
  float sigX = pitchX / std::sqrt(12);
  float sigY = pitchY / std::sqrt(12);
  clu.setErr(sigX * sigX, 0., sigY * sigY);
  // add MC labels
  int nCluLabels = 0;
  for (const auto& pix : mPixArrBuff) {
    if (nCluLabels >= MaxLabels)
      break;
    int jDig = pix.originalIndex;
    std::span labels = mDigMCLabelsPtr->getLabels(jDig);
    int nLabels = labels.size();
    for (int jLab = 0; jLab < nLabels; jLab++) {
      NA6PMCComposedLabel lbl = labels[jLab];
      // Check if this label has already been added to our cluster accumulator
      bool isDuplicate = false;
      for (int jLabC = 0; jLabC < nCluLabels; ++jLabC) {
        if (lbl == mLabelsBuff[jLabC]) {
          isDuplicate = true;
          break;
        }
      }
      if (!isDuplicate) {
        mLabelsBuff[nCluLabels++] = lbl;
        if (nCluLabels >= MaxLabels)
          break;
      }
    }
  }
  for (int jLabC = 0; jLabC < nCluLabels; ++jLabC) {
    mCluMCLabelsPtr->addElement(cluID, mLabelsBuff[jLabC]);
  }
  if (nCluLabels >= 0) {
    clu.setParticleID(mLabelsBuff[0].getTrackID());
  }
  return true;
}

void NA6PVerTelClusterizer::getClusterGlobalCoord(int modID, float xloc, float yloc, double xyzGlo[3]) const
{
  auto& matrix = mGeoManager.getMatrix(modID);
  double xyzLoc[3];
  xyzLoc[0] = xloc - mGeoManager.getModuleHalfX(modID);
  xyzLoc[1] = yloc - mGeoManager.getModuleHalfY(modID);
  xyzLoc[2] = 0.;
  matrix.LocalToMaster(xyzLoc, xyzGlo);
}
