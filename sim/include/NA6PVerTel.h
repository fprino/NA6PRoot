// NA6PCCopyright

#ifndef NA6P_VERTEL_H_
#define NA6P_VERTEL_H_

#include "NA6PModule.h"
#include "NA6PVerTelHit.h"
#include <vector>

class TFile;
class TTree;
class TGeoVolume;

class NA6PVerTel : public NA6PModule
{
 public:
  NA6PVerTel() : NA6PModule("VerTel") { setActiveID(0); }
  ~NA6PVerTel() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume* world) override;
  bool stepManager(int volID) override;
  size_t getNHits() const override { return mHits.size(); }
  void createHitsOutput(const std::string& outDir) override;
  void closeHitsOutput() override;
  void writeHits(const std::vector<int>& remapping) override;
  void setAlignableEntries() override;
  void clearHits() override { mHits.clear(); }

  const auto& getHits() const { return mHits; }

  NA6PVerTelHit* addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom, const TVector3& endMom,
                        float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus);

 private:
  static constexpr int NChipsPerStation = 4;

  std::vector<NA6PVerTelHit> mHits, *hHitsPtr = &mHits;
  TFile* mHitsFile = nullptr;
  TTree* mHitsTree = nullptr;
};

#endif