// NA6PCCopyright

// Based on:
// Copyright 2019-2030 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_VERTEL_CLUSTERIZER_H
#define NA6P_VERTEL_CLUSTERIZER_H

#include "NA6PVerTelSegmentation.h"
#include "NA6PGeometryManager.h"
#include "NA6PVerTelDigit.h"
#include "NA6PVerTelCluster.h"
#include "NA6PMCTruthContainer.h"
#include <Rtypes.h>

struct DigitProxy {
  int originalIndex;
  int chipID;
  uint32_t col;
  uint32_t row;
};

struct ChipSlice {
  int chipID;
  size_t startIdx;
  size_t endIdx;
};

struct PreCluster {
  int head = 0; // index of precluster head in the pixels
  int index = 0;
};

struct BBox {
  uint32_t rowMin = 0xffffffff;
  uint32_t rowMax = 0;
  uint32_t colMin = 0xffffffff;
  uint32_t colMax = 0;

  bool isInside(uint32_t row, uint32_t col) const { return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax; }
  auto rowSpan() const { return rowMax - rowMin + 1; }
  auto colSpan() const { return colMax - colMin + 1; }
  bool isAcceptableSize() const
  {
    // Always returns true for now since we are skipping the huge cluster treatment
    return true;
  }
  void clear()
  {
    rowMin = colMin = 0xffffffff;
    rowMax = colMax = 0;
  }
  void adjust(uint32_t row, uint32_t col)
  {
    if (row < rowMin)
      rowMin = row;
    if (row > rowMax)
      rowMax = row;
    if (col < colMin)
      colMin = col;
    if (col > colMax)
      colMax = col;
  }
};

class NA6PVerTelClusterizer
{
 public:
  NA6PVerTelClusterizer();
  ~NA6PVerTelClusterizer() = default;

  void initGeometry(const char* filename = "geometry.root", const char* geoname = "NA6P");

  void resetColumn(int* buff)
  {
    std::fill_n(buff, NA6PVerTelSegmentation::NRowsPerTile, -1);
  }

  ///< swap current and previous column buffers
  void swapColumnBuffers() { std::swap(mPrevColDataPtr, mCurrColDataPtr); }

  ///< add new precluster at given row of current column for the fired pixel with index ip in chipdigits
  void addNewPrecluster(uint32_t jPix, uint32_t row)
  {
    int lastIndex = mPreClusters.size();
    mPreClusters.emplace_back(mPixels.size(), lastIndex);
    // new head does not point yet (-1) on other pixels, store just the entry of the pixel in chipdigits
    mPixels.emplace_back(-1, jPix);
    mCurrColDataPtr[row] = lastIndex; // store index of the new precluster in the current column buffer
  }

  ///< add pixel at row (entry jPix in chipdigits) to the precluster with given index
  void expandPreCluster(uint32_t jPix, uint32_t row, int preClusIndex)
  {
    int masterIndex = mPreClusters[preClusIndex].index;
    int firstIndex = mPreClusters[masterIndex].head;
    mPixels.emplace_back(firstIndex, jPix);
    mPreClusters[masterIndex].head = mPixels.size() - 1;
    mCurrColDataPtr[row] = preClusIndex;
  }

  void process(const std::vector<NA6PVerTelDigit>& vtDigits, const NA6PMCTruthContainer& mcDigLabels, std::vector<NA6PVerTelCluster>& outClusters, NA6PMCTruthContainer& mcCluLabels);
  void processChip(size_t first, size_t last);
  void initChip(uint32_t jPix);
  void updateChip(uint32_t jPix);
  void finishChip();
  bool pixelsToRecPoint();
  void getClusterGlobalCoord(int modID, float xloc, float yloc, double xyzGlo[3]) const;

 protected:
  static constexpr int nRSU = NA6PVerTelSegmentation::NYSensors * NA6PVerTelSegmentation::NXSegments / NA6PVerTelSegmentation::NSegmentsPerRSU;
  static constexpr int nTilesPerRSU = NA6PVerTelSegmentation::NTilesPerRSU;
  static constexpr int nTilesPerModule = NA6PVerTelSegmentation::NTilesPerModule;
  static constexpr float pitchX = NA6PVerTelSegmentation::ActiveDX / NA6PVerTelSegmentation::NColsPerTile;
  static constexpr float pitchY = NA6PVerTelSegmentation::ActiveDYTile / NA6PVerTelSegmentation::NRowsPerTile;
  static constexpr int MaxLabels = 10;

  const std::vector<NA6PVerTelDigit>* mDigitsPtr = nullptr; // transient pointer to digits
  const NA6PMCTruthContainer* mDigMCLabelsPtr = nullptr;    // transient pointer to MC truth of digits
  std::vector<NA6PVerTelCluster>* mClustersPtr = nullptr;   // transient pointer to output clusters
  NA6PMCTruthContainer* mCluMCLabelsPtr = nullptr;          // transient pointer to MC truth of clusters
  NA6PVerTelSegmentation mSegmentation;                     ///< segmentation class
  NA6PGeometryManager mGeoManager;                          ///< geometry manager
  std::vector<DigitProxy> mChipDigits;                      // vector to store compacted digits during clusterization
  std::vector<ChipSlice> mChipSlices;                       // vector to store limits of chips inside mChipDigits
  // buffers for entries in preClusterIndices in 2 columns, to avoid boundary checks, we reserve
  // extra elements in the beginning and the end
  std::array<int, NA6PVerTelSegmentation::NRowsPerTile + 2> mColumn1;
  std::array<int, NA6PVerTelSegmentation::NRowsPerTile + 2> mColumn2;
  int* mCurrColDataPtr = nullptr; // pointer on the 1st row of currently processed columnsX
  int* mPrevColDataPtr = nullptr; // pointer on the 1st row of previously processed columnsX
  bool mNoLeftCol = true;         ///< flag that there is no column on the left to check
  uint32_t mCurrCol = 0xffff;     ///< Column being processed
  // mPixels[].first is the index of the next pixel of the same precluster in the pixels
  // mPixels[].second is the index of the referred pixel in the chipdigits vector
  std::vector<std::pair<int, uint32_t>> mPixels;
  std::vector<PreCluster> mPreClusters;                   //! preclusters info
  std::vector<DigitProxy> mPixArrBuff;                    //! temporary buffer for pattern calc.
  std::array<NA6PMCComposedLabel, MaxLabels> mLabelsBuff; //! temporary buffer for building cluster labels

  ClassDefNV(NA6PVerTelClusterizer, 1);
};

#endif
