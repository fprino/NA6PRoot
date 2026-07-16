// NA6PCCopyright

#include "NA6PVerTel.h"
#include "NA6PDetector.h"
#include "NA6PPixelStation.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"

#include <TVirtualMC.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TColor.h>
#include <fairlogger/Logger.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

void NA6PVerTel::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;
  nameM = addName("Silicon");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kCyan + 1);
  }
  nameM = addName("CarbonFoam");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 1, 0.5);
    mixt->AddElement(12.01, 6, 1.0);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue - 6);
  }
  nameM = addName("CarbonFiber");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 1, 1.91);
    mixt->AddElement(12.01, 6, 1.0);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 1);
  }
  nameM = addName("Air");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM);
  }
  nameM = addName("Aluminium");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 26.98, 13, 2.7);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kAzure - 9);
  }
  nameM = addName("Steel");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 55.84, 26, 7.87);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 2);
  }
  nameM = addName("Water");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 1.0);
    mixt->AddElement(new TGeoElement("H", "Hydrogen", 1, 1.008), 2);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 1);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue);
  }
}

void NA6PVerTel::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();

  TGeoMedium* medAir = NA6PTGeoHelper::instance().getMedium(addName("Air"));

  // ── VT container ─────────────────────────────────────────────────
  const float vtZMargin = NA6PPixelStation::PixChipContainerDZ;
  const float downstreamCarbonFullZ = NA6PPixelStation::CarbonPlateDz + NA6PPixelStation::FrameGlueGap;
  const float vtFullZ = param.posVerTelPlaneZ[param.nVerTelPlanes - 1] - param.posVerTelPlaneZ[0] + 2 * vtZMargin + 2 * NA6PPixelStation::FrameHZ + downstreamCarbonFullZ;
  const float vtHalfZ = vtFullZ / 2.f;

  const float maxCarbonPlateDX = 2.f * NA6PPixelStation::FrameHX;
  const float maxCarbonPlateDY = 2.f * NA6PPixelStation::FrameHY;
  const float vtHalfX = TMath::Max(static_cast<float>(NA6PPixelStation::FrameHX + 4.5), maxCarbonPlateDX / 2.f + 1.f);
  const float vtHalfY = TMath::Max(static_cast<float>(NA6PPixelStation::FrameHY + 3.f), maxCarbonPlateDY / 2.f + 1.f);

  auto* vtShape = new TGeoBBox("VTContainer",
                               vtHalfX,
                               vtHalfY,
                               vtHalfZ);
  TGeoVolume* vtContainer = new TGeoVolume("VTContainer", vtShape, medAir);

  const float vtZCenterGlobal = param.posVerTelPlaneZ[0] + vtHalfZ - vtZMargin;
  auto* vtTransform = new TGeoCombiTrans(
    param.shiftVerTel[0],
    param.shiftVerTel[1],
    param.shiftVerTel[2] + vtZCenterGlobal,
    NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
  world->AddNode(vtContainer, composeNonSensorVolID(0), vtTransform);

  // ── Per-plane loop ────────────────────────────────────────────────
  const float vtZOffset = vtZCenterGlobal;
  const float frameRelZ = static_cast<float>(NA6PPixelStation::PixChipDz / 2 + NA6PPixelStation::CarbonPlateDz + NA6PPixelStation::FrameHZ + NA6PPixelStation::FrameGlueGap);
  const float downstreamCarbonRelZ = frameRelZ + NA6PPixelStation::FrameHZ + NA6PPixelStation::FrameGlueGap + NA6PPixelStation::CarbonPlateDz / 2.f;
  const float backSensorRelZ = downstreamCarbonRelZ + NA6PPixelStation::CarbonPlateDz / 2.f + NA6PPixelStation::PixChipDz / 2.f;

  const NA6PPixelStation::Materials stationMaterials = {
    addName("Silicon"),
    addName("CarbonFiber"),
    addName("CarbonFoam"),
    addName("Air"),
    addName("Steel"),
    addName("Water"),
    addName("Aluminium")};
  const NA6PPixelStation pixelStation(*this, stationMaterials);

  for (int ll = 0; ll < param.nVerTelPlanes; ++ll) {
    const NA6PPixelStation::Placement placement = {
      param.posVerTelPlaneX[ll],
      param.posVerTelPlaneY[ll],
      param.posVerTelPlaneZ[ll] - vtZOffset,
      param.posVerTelPlaneZ[ll] - vtZOffset + frameRelZ,
      param.posVerTelPlaneZ[ll] - vtZOffset + downstreamCarbonRelZ,
      param.posVerTelPlaneZ[ll] - vtZOffset + backSensorRelZ};
    pixelStation.addTo(vtContainer, ll, placement, param.sensorsPerPlane[ll]);
  }
}

void NA6PVerTel::setAlignableEntries()
{
  const auto& param = NA6PLayoutParam::Instance();
  std::string topNodeName = gGeoManager->GetTopNode()->GetName();
  int svolCnt = 0;
  for (int ll = 0; ll < param.nVerTelPlanes; ++ll) {
    for (int ii = 0; ii < 4; ++ii) {
      int id = getActiveID() * 100 + svolCnt;
      std::string nm = fmt::format("VT_Lr{}_Sens{}", ll, ii);
      std::string path = fmt::format("/{}/VTContainer_{}/PixelStationVol_Pl{}_{}/PixelSensor_Pl{}_{}",
                                     topNodeName,
                                     composeNonSensorVolID(0),
                                     ll, composeNonSensorVolID(ll),
                                     ll, composeSensorVolID(ii));
      TGeoPNEntry* entry = gGeoManager->SetAlignableEntry(nm.c_str(), path.c_str(), id);
      if (entry) {
        LOGP(info, "Successfully added {} {} as alignable sensor {}", nm, path, id);
      } else {
        LOGP(error, "FAILED to add alignable entry {} {}", nm, path);
      }
      svolCnt++;
    }
  }
}

bool NA6PVerTel::stepManager(int volID)
{
  int sensID = NA6PModule::volID2SensID(volID);
  if (sensID < 0) {
    LOGP(fatal, "Non-sensor volID={} was provided to stepManager of {}", volID, getName());
  }
  auto mc = TVirtualMC::GetMC();
  bool startHit = false, stopHit = false;
  unsigned char status = 0;
  if (mc->IsTrackEntering()) {
    status |= NA6PBaseHit::kTrackEntering;
  }
  if (mc->IsTrackInside()) {
    status |= NA6PBaseHit::kTrackInside;
  }
  if (mc->IsTrackExiting()) {
    status |= NA6PBaseHit::kTrackExiting;
  }
  if (mc->IsTrackOut()) {
    status |= NA6PBaseHit::kTrackOut;
  }
  if (mc->IsTrackStop()) {
    status |= NA6PBaseHit::kTrackStopped;
  }
  if (mc->IsTrackAlive()) {
    status |= NA6PBaseHit::kTrackAlive;
  }

  if ((status & NA6PBaseHit::kTrackEntering) || (status & NA6PBaseHit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((status & (NA6PBaseHit::kTrackExiting | NA6PBaseHit::kTrackOut | NA6PBaseHit::kTrackStopped))) {
    stopHit = true;
  }
  if (!startHit) {
    mTrackData.mEnergyLoss += mc->Edep();
  }
  if (!(startHit | stopHit)) {
    return false; // do nothing
  }
  auto stack = (NA6PMCStack*)mc->GetStack();
  if (startHit) {
    mTrackData.mEnergyLoss = 0.;
    mc->TrackMomentum(mTrackData.mMomentumStart);
    mc->TrackPosition(mTrackData.mPositionStart);
    mTrackData.mTrkStatusStart = status;
    mTrackData.mHitStarted = true;
  }
  if (stopHit) {
    TLorentzVector positionStop, momentumStop;
    mc->TrackMomentum(momentumStop);
    mc->TrackPosition(positionStop);
    int stationID(-1);
    mc->CurrentVolOffID(1, stationID);
    stationID = volID2NonSensID(stationID);

    int chipindex = NChipsPerStation * stationID + sensID;
    auto* p = addHit(stack->GetCurrentTrackNumber(), chipindex, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                     mTrackData.mMomentumStart.Vect(), momentumStop.Vect(), positionStop.T(),
                     mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart, status);
    if (mVerbosity > 0) {
      LOGP(info, "{} Tr{} {}", getName(), stack->GetCurrentTrackNumber(), p->asString());
    }
    // register the points in TParticle
    stack->addHit(getActiveIDBit());
    return true;
  }
  return false;
}

NA6PVerTelHit* NA6PVerTel::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom, const TVector3& endMom,
                                  float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endMom, endTime, eLoss, startStatus, endStatus);
  return &(mHits.back());
}

void NA6PVerTel::createHitsOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}Hits{}.root", outDir, getName());
  mHitsFile = TFile::Open(nm.c_str(), "recreate");
  mHitsTree = new TTree(fmt::format("hits{}", getName()).c_str(), fmt::format("{} Hits", getName()).c_str());
  mHitsTree->Branch(getName().c_str(), &hHitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), nm);
}

void NA6PVerTel::closeHitsOutput()
{
  if (mHitsTree && mHitsFile) {
    mHitsFile->cd();
    mHitsTree->Write();
    delete mHitsTree;
    mHitsTree = 0;
    mHitsFile->Close();
    delete mHitsFile;
    mHitsFile = 0;
  }
}

void NA6PVerTel::writeHits(const std::vector<int>& remapping)
{
  int nh = mHits.size();
  for (int i = 0; i < nh; i++) {
    auto& h = mHits[i];
    if (remapping[h.getTrackID()] < 0) {
      LOGP(error, "Track {} hit {} in {} was not remapped!", h.getTrackID(), i, getName());
    }
    h.setTrackID(remapping[h.getTrackID()]);
  }
  if (mHitsTree) {
    mHitsTree->Fill();
  }
}
