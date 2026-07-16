// NA6PCCopyright

#include "NA6PVerTelSegmentation.h"
#include "NA6PGeometryManager.h"
#include "NA6PLayoutParam.h"

NA6PVerTelSegmentation::NA6PVerTelSegmentation()
{
  const auto& layout = NA6PLayoutParam::Instance();
  mOffsXFront = layout.pixChipOffsXFront;
  mOffsXBack = layout.pixChipOffsXBack;
  mOffsX = mOffsXFront;
  mOffsY = -layout.pixChipOffsY;
}

void NA6PVerTelSegmentation::setDetectorID(int detID)
{
  mOffsX = isBackChip(detID) ? mOffsXBack : mOffsXFront;
}

bool NA6PVerTelSegmentation::isBackChip(int detID)
{
  const int chipInPlane = detID % NA6PGeometryManager::kNVTModulesPerLayer;
  return chipInPlane == 0 || chipInPlane == 2;
}

void NA6PVerTelSegmentation::setStaggered(bool val)
{
  // emulate the configuration with sensor acceptance overlaps
  mStaggered = val;
  if (val) {
    mDeadXLongEff = DeadXLong + DeadXShort;
    mDeadXShortEff = 0.f;
    mDeadYTopEff = DeadYTop + DeadYBottom;
    mDeadYBottomEff = 0.f;
    mInterChipGap = 0.f;
  } else {
    mDeadXLongEff = DeadXLong;
    mDeadXShortEff = DeadXShort;
    mDeadYTopEff = DeadYTop;
    mDeadYBottomEff = DeadYBottom;
  }
}

bool NA6PVerTelSegmentation::computePixelIndices(float xloc, float yloc,
                                                 float deadYBottom, float deadYTop,
                                                 UShort_t& rsu, UShort_t& tile,
                                                 UShort_t& row, UShort_t& col) const
{
  int sens = int(yloc / DYSens);
  yloc -= sens * DYSens;
  if (yloc < deadYBottom || yloc > DYSens - deadYTop)
    return false;
  yloc -= deadYBottom;

  int topbot = int(yloc / ActiveDYHalf);
  yloc -= topbot * ActiveDYHalf;
  if (yloc < DeadTopBotHalves / 2)
    return false;
  yloc -= DeadTopBotHalves / 2;
  int segm = int(xloc / DXSegment);
  xloc -= segm * DXSegment;
  if (xloc > DXSegment - DeadXDataBackbone)
    return false;
  int tileIdx = int(xloc / DXTile);
  xloc -= tileIdx * DXTile;
  if (xloc < DeadXTile)
    return false;
  xloc -= DeadXTile;

  // map segm [0-11] and topbot [0-1] and tileIdx [0-2] to rsu [0-41] and tile [0-11]
  int rsuX = segm / NSegmentsPerRSU;
  int segInRsu = segm % NSegmentsPerRSU;
  rsu = static_cast<UShort_t>(sens * (NXSegments / NSegmentsPerRSU) + rsuX);
  tile = static_cast<UShort_t>(topbot * (NSegmentsPerRSU * NTilesPerSegment) + segInRsu * NTilesPerSegment + tileIdx);
  col = static_cast<UShort_t>(xloc / ActiveDX * NColsPerTile);
  row = static_cast<UShort_t>(yloc / ActiveDYTile * NRowsPerTile);
  return true;
}

int NA6PVerTelSegmentation::isInAcc(float xloc, float yloc) const
{
  // method used in hitsToRecPoints to account for inactive regions of MOSAIX

  xloc -= mInterChipGap / 2;
  yloc -= mInterChipGap / 2;

  if (xloc < 0 || xloc > XSizeTot || yloc < 0 || yloc > YSizeTot)
    return 0;

  if (xloc < mDeadXShortEff || xloc > XSizeTot - mDeadXLongEff || yloc < mDeadYBottomEff || yloc > YSizeTot - mDeadYTopEff)
    return -1;
  xloc -= mDeadXShortEff;

  UShort_t rsu, tile, row, col;
  return computePixelIndices(xloc, yloc, mDeadYBottomEff, mDeadYTopEff, rsu, tile, row, col) ? 1 : -1;
}

bool NA6PVerTelSegmentation::localToIndices(float xloc, float yloc, UShort_t& rsu, UShort_t& tile, UShort_t& row, UShort_t& col) const
{
  rsu = tile = row = col = -1;
  if (xloc < 0 || xloc > XSizeTot || yloc < 0 || yloc > YSizeTot)
    return false;

  if (xloc < DeadXShort || xloc > XSizeTot - DeadXLong || yloc < DeadYBottom || yloc > YSizeTot - DeadYTop)
    return false;

  xloc -= DeadXShort;

  bool isOk = computePixelIndices(xloc, yloc, DeadYBottom, DeadYTop, rsu, tile, row, col);
  if (!isOk || col >= NColsPerTile || row >= NRowsPerTile || rsu >= NYSensors * (NXSegments / NSegmentsPerRSU) || tile >= NTilesPerRSU) {
    return false;
  }
  return true;
}

bool NA6PVerTelSegmentation::indicesToLocal(UShort_t rsu, UShort_t tile, UShort_t row, UShort_t col, float& xloc, float& yloc) const
{
  xloc = yloc = -9999.f;
  if (rsu >= NYSensors * (NXSegments / NSegmentsPerRSU))
    return false;
  if (tile >= NTilesPerSegment * NSegmentsPerRSU * 2)
    return false;
  if (row >= NRowsPerTile)
    return false;
  if (col >= NColsPerTile)
    return false;

  float xInTile = (col + 0.5f) * ActiveDX / NColsPerTile;
  float yInTile = (row + 0.5f) * ActiveDYTile / NRowsPerTile;
  int tileIdx = tile % NTilesPerSegment;
  int segInRsu = (tile / NTilesPerSegment) % NSegmentsPerRSU;
  int topbot = tile / (NTilesPerSegment * NSegmentsPerRSU);
  float xInSegm = xInTile + DeadXTile + tileIdx * DXTile;
  int rsuX = rsu % (NXSegments / NSegmentsPerRSU);
  int segm = NSegmentsPerRSU * rsuX + segInRsu;
  xloc = xInSegm + segm * DXSegment + DeadXShort;
  float yinRsu = yInTile + DeadTopBotHalves / 2 + topbot * ActiveDYHalf;
  int rsuY = rsu / (NXSegments / NSegmentsPerRSU);
  yloc = yinRsu + rsuY * DYSens + DeadYBottom;
  return true;
}
