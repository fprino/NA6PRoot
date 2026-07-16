// NA6PCCopyright

#ifndef NA6P_COOLING_PLATE_H_
#define NA6P_COOLING_PLATE_H_

#include "NA6PPixelStation.h"

#include <TGeoMatrix.h>
#include <RtypesCore.h>

class TGeoMedium;
class TGeoVolume;

struct NA6PStationGeometryParams {
  Double_t halfFrameX = NA6PPixelStation::FrameHX;
  Double_t halfFrameY = NA6PPixelStation::FrameHY;
  Double_t halfFrameZ = NA6PPixelStation::FrameHZ;
  Double_t halfAlCutZ = 0.5;
  Double_t outerCarbonFoamDz = 0.406;
  Double_t outerCarbonFiberDz = NA6PPixelStation::CarbonPlateDz;
  Double_t pipeZ = 0.0;
  Double_t grooveFitClear = NA6PPixelStation::kGeomEps;
  Double_t innerTubeOffsetX = 10.9 / 2.0;
  Double_t outerTubeOffsetX = 11.3 / 2.0;
  Double_t tubeOffsetY = 14.085;
  Double_t outerOutletAngleDeg = 65.0;
  Double_t fiberCenterHoleR = 0.3;

  Double_t outerCarbonFoamHalfZ = outerCarbonFoamDz / 2.0;
  Double_t outerCarbonFiberHalfZ = outerCarbonFiberDz / 2.0;
  Double_t outerFrontCarbonFiberSurfaceZ = -outerCarbonFoamHalfZ - outerCarbonFiberDz;
  Double_t outerBackCarbonFiberSurfaceZ = outerCarbonFoamHalfZ + outerCarbonFiberDz;
  Double_t outerSensorPocketDz = halfFrameZ - outerBackCarbonFiberSurfaceZ;
  Double_t outerFrontSensorZ = outerFrontCarbonFiberSurfaceZ - NA6PPixelStation::FrameGlueGap - NA6PPixelStation::PixChipDz / 2.0;
  Double_t outerBackSensorZ = outerBackCarbonFiberSurfaceZ + NA6PPixelStation::FrameGlueGap + NA6PPixelStation::PixChipDz / 2.0;
};

struct NA6PCoolingPipeParams {
  Double_t innerRadius = 0.125;
  Double_t outerRadius = 0.150;
  TGeoMedium* medSteel = nullptr;
  TGeoMedium* medWater = nullptr;
};

struct NA6PCoolingPlateMaterials {
  TGeoMedium* plate = nullptr;
  TGeoMedium* carbonFoam = nullptr;
  TGeoMedium* carbonFiber = nullptr;
  TGeoMedium* air = nullptr;
};

struct NA6PPipeSegment {
  enum class Kind {
    Straight,
    Curved
  };

  Kind kind = Kind::Straight;
  Double_t halfLength = 0.0;
  Double_t rBend = 0.0;
  Double_t phiStart = 0.0;
  Double_t phiRange = 0.0;
  Double_t tx = 0.0;
  Double_t ty = 0.0;
  Double_t phi = 90.0;
  Double_t theta = 90.0;
};

TGeoVolume* BuildNA6PCoolingPlate(NA6PPixelStation::CoolingPlateType type,
                                  const char* tag,
                                  const NA6PCoolingPipeParams& pipeP,
                                  const NA6PCoolingPlateMaterials& materials,
                                  const NA6PStationGeometryParams& geom);

#endif
