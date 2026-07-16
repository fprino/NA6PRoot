// NA6PCCopyright

#ifndef NA6P_PIXEL_STATION_H_
#define NA6P_PIXEL_STATION_H_

#include <string>

class NA6PModule;
class TGeoVolume;

class NA6PPixelStation
{
 public:
  enum class CoolingPlateType {
    FirstStations,
    OuterStations
  };

  struct Materials {
    std::string silicon;
    std::string carbonFiber;
    std::string carbonFoam;
    std::string air;
    std::string steel;
    std::string water;
    std::string coolingPlate;
  };

  struct Placement {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    float frameZ = 0.f;
    float downstreamCarbonZ = 0.f;
    float backSensorZ = 0.f;
  };

  static constexpr float PixChipContainerDX = 30.0f;
  static constexpr float PixChipContainerDY = 30.0f;
  static constexpr float PixChipContainerDZ = 1.f;
  static constexpr float PixChipDz = 50e-4f;
  static constexpr float CarbonPlateDz = 400e-4f;
  static constexpr double FrameHX = 19.0;
  static constexpr double FrameHY = 16.0;
  static constexpr double FrameHZ = 0.4;
  static constexpr float FrameGlueGap = 0.01f;
  static constexpr double kGeomEps = 0.01; // linear overshoot/clearance, cm
  static constexpr double kAngEps = 0.01;  // angular overshoot for TGeoTubeSeg phi bounds, deg

  NA6PPixelStation(const NA6PModule& module, const Materials& materials);

  void addTo(TGeoVolume* vtContainer, int layer, const Placement& placement, int sensorsPerPlane) const;
  static CoolingPlateType coolingPlateTypeForLayer(int layer);

 private:
  const NA6PModule& mModule;
  Materials mMaterials;
};

#endif
