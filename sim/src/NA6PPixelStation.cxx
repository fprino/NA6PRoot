// NA6PCCopyright

#include "NA6PPixelStation.h"

#include "NA6PCoolingPlate.h"
#include "NA6PLayoutParam.h"
#include "NA6PModule.h"
#include "NA6PTGeoHelper.h"
#include "NA6PVerTelSegmentation.h"

#include <TColor.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <TString.h>

#include <array>

namespace
{
constexpr int NSensorsPerStation = 4;
} // namespace

NA6PPixelStation::NA6PPixelStation(const NA6PModule& module, const Materials& materials)
  : mModule(module), mMaterials(materials)
{
}

NA6PPixelStation::CoolingPlateType NA6PPixelStation::coolingPlateTypeForLayer(int layer)
{
  return layer < 2 ? CoolingPlateType::FirstStations : CoolingPlateType::OuterStations;
}

void NA6PPixelStation::addTo(TGeoVolume* vtContainer, int layer, const Placement& placement, int sensorsPerPlane) const
{
  const auto& layout = NA6PLayoutParam::Instance();
  const NA6PStationGeometryParams geom;

  const float pixChipFullYLayer = sensorsPerPlane * NA6PVerTelSegmentation::DYSens;
  const float pixChipHalfX = NA6PVerTelSegmentation::XSizeTot / 2.f;
  const float pixChipHalfYLayer = pixChipFullYLayer / 2.f;
  const float pixChipHalfZ = PixChipDz / 2.f;

  TString sensorShName = TString::Format("SensorShape_Pl%d", layer);
  TString stationVName = TString::Format("PixelStationVol_Pl%d", layer);
  TString sensorVName = TString::Format("PixelSensor_Pl%d", layer);

  const float frontSensorZ = layer < 2 ? placement.z : placement.frameZ + geom.outerFrontSensorZ;
  const float backSensorZ = layer < 2 ? placement.backSensorZ : placement.frameZ + geom.outerBackSensorZ;
  const float stationCenterZ = 0.5f * (frontSensorZ + backSensorZ);
  const float frontSensorLocalZ = frontSensorZ - stationCenterZ;
  const float backSensorLocalZ = backSensorZ - stationCenterZ;
  auto* sensorShape = new TGeoBBox(sensorShName, pixChipHalfX, pixChipHalfYLayer, pixChipHalfZ);

  // Contract with NA6PVerTel::setAlignableEntries(): PixelStationVol_PlX must
  // contain PixelSensor_PlX_Y as direct daughters for the fixed alignable path.
  TGeoVolume* pixelStationVol = new TGeoVolumeAssembly(stationVName);
  auto* pixelSensor = new TGeoVolume(sensorVName, sensorShape, NA6PTGeoHelper::instance().getMedium(mMaterials.silicon));
  pixelSensor->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(mMaterials.silicon));

  const std::array<float, NSensorsPerStation> alpdx = {pixChipHalfX + layout.pixChipOffsXBack, -pixChipHalfX + layout.pixChipOffsXFront,
                                                       -pixChipHalfX - layout.pixChipOffsXBack, pixChipHalfX - layout.pixChipOffsXFront};
  const std::array<float, NSensorsPerStation> alpdy = {pixChipHalfYLayer - layout.pixChipOffsY, pixChipHalfYLayer + layout.pixChipOffsY,
                                                       -pixChipHalfYLayer + layout.pixChipOffsY, -pixChipHalfYLayer - layout.pixChipOffsY};

  std::vector<float> phi = {0., 90., 180., 0.};
  std::vector<float> theta = {0., 180., 0., 180.};
  std::vector<float> psi = {0., -90., 0., 0.};

  for (size_t ii = 0; ii < alpdx.size(); ++ii) {
    const float sensorLocalZ = ((ii == 0 || ii == 2)) ? backSensorLocalZ : frontSensorLocalZ;
    pixelStationVol->AddNode(pixelSensor, mModule.composeSensorVolID(ii),
                             new TGeoCombiTrans(alpdx[ii], alpdy[ii], sensorLocalZ, new TGeoRotation(Form("rotQ%d", ii + 1), phi[ii], theta[ii], psi[ii])));
  }
  vtContainer->AddNode(pixelStationVol, mModule.composeNonSensorVolID(layer),
                       new TGeoCombiTrans(placement.x, placement.y, stationCenterZ,
                                          NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0)));

  NA6PCoolingPipeParams pipeP;
  pipeP.medSteel = NA6PTGeoHelper::instance().getMedium(mMaterials.steel);
  pipeP.medWater = NA6PTGeoHelper::instance().getMedium(mMaterials.water);

  NA6PCoolingPlateMaterials plateMaterials;
  plateMaterials.plate = NA6PTGeoHelper::instance().getMedium(mMaterials.coolingPlate);
  plateMaterials.carbonFoam = NA6PTGeoHelper::instance().getMedium(mMaterials.carbonFoam);
  plateMaterials.carbonFiber = NA6PTGeoHelper::instance().getMedium(mMaterials.carbonFiber);
  plateMaterials.air = NA6PTGeoHelper::instance().getMedium(mMaterials.air);

  TString planeTag = TString::Format("Pl%d", layer);
  TGeoVolume* alAsm = BuildNA6PCoolingPlate(coolingPlateTypeForLayer(layer),
                                            planeTag.Data(),
                                            pipeP,
                                            plateMaterials,
                                            geom);

  vtContainer->AddNode(alAsm, mModule.composeNonSensorVolID(layer + 40),
                       new TGeoCombiTrans(placement.x, placement.y, placement.frameZ,
                                          NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0)));
}
