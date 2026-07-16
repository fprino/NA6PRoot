// NA6PCCopyright

#include "NA6PCoolingPlate.h"
#include "NA6PLayoutParam.h"

#include <TColor.h>
#include <TGeoBBox.h>
#include <TGeoBoolNode.h>
#include <TGeoCompositeShape.h>
#include <TGeoTube.h>
#include <TGeoTorus.h>
#include <TGeoVolume.h>
#include <TGeoXtru.h>
#include <TMath.h>
#include <TString.h>

#include <cmath>
#include <vector>

namespace
{
Int_t gPipeNodeIdx = 0;
std::vector<TString> gGrooveTerms;
std::vector<TString> gGrooveTermsOuter;
std::vector<TString>* gActiveGrooveTerms = &gGrooveTerms;
TString gGrooveNamePrefix = "Inner";

void AddStraightTube(TGeoVolume* assembly,
                     const NA6PCoolingPipeParams& p,
                     Double_t halfLength,
                     Double_t tx, Double_t ty,
                     Double_t phi = 90,
                     Double_t theta = 90,
                     Double_t pipeZ = 0.0,
                     Double_t halfFrameZ = NA6PPixelStation::FrameHZ,
                     Double_t fitClear = NA6PPixelStation::kGeomEps)
{
  Int_t idx = ++gPipeNodeIdx;
  TString tag = TString::Format("%s_%04d", gGrooveNamePrefix.Data(), idx);

  auto* shSteel = new TGeoTube("shTubeS_" + tag, p.innerRadius, p.outerRadius, halfLength);
  auto* volSteel = new TGeoVolume("VolTubeS_" + tag, shSteel, p.medSteel);
  volSteel->SetLineColor(kGray + 2);

  auto* shWater = new TGeoTube("shTubeW_" + tag, 0, p.innerRadius, halfLength);
  auto* volWater = new TGeoVolume("VolTubeW_" + tag, shWater, p.medWater);
  volWater->SetLineColor(kBlue);

  auto* rot = new TGeoRotation("rotT_" + tag, phi, theta, 0);
  assembly->AddNode(volSteel, idx, new TGeoCombiTrans(tx, ty, pipeZ, rot));
  assembly->AddNode(volWater, idx, new TGeoCombiTrans(tx, ty, pipeZ, rot));

  // Fit-clearance cylinder around the pipe, used to cut its housing groove
  // out of the frame/foam/fiber stack.
  auto* shFit = new TGeoTube("GrFit_" + tag, 0, p.outerRadius + fitClear, halfLength);
  auto* ctCutTube = new TGeoCombiTrans("GrFitT_" + tag, tx, ty, pipeZ, rot);
  ctCutTube->RegisterYourself();
  gActiveGrooveTerms->push_back(TString::Format("%s:%s", shFit->GetName(), ctCutTube->GetName()));

  // Channel box, extended by kGeomEps past the tube end to avoid a
  // coincident-surface artifact at the groove termination.
  auto* shChan = new TGeoBBox("GrChan_" + tag,
                              p.outerRadius + fitClear,
                              halfFrameZ / 2. + NA6PPixelStation::kGeomEps,
                              halfLength + NA6PPixelStation::kGeomEps);
  auto* ctCutBox = new TGeoCombiTrans("GrChanT_" + tag,
                                      tx, ty, pipeZ + halfFrameZ / 2., rot);
  ctCutBox->RegisterYourself();
  gActiveGrooveTerms->push_back(TString::Format("%s:%s", shChan->GetName(), ctCutBox->GetName()));
}

void AddCurvedTube(TGeoVolume* assembly,
                   const NA6PCoolingPipeParams& p,
                   Double_t rBend,
                   Double_t phiStart, Double_t phiRange,
                   Double_t tx, Double_t ty,
                   Double_t pipeZ = 0.0,
                   Double_t halfFrameZ = NA6PPixelStation::FrameHZ,
                   Double_t chanWidth = 0.125,
                   Double_t fitClear = NA6PPixelStation::kGeomEps)
{
  Int_t idx = ++gPipeNodeIdx;
  TString tag = TString::Format("%s_%04d", gGrooveNamePrefix.Data(), idx);

  auto* shSteel = new TGeoTorus("shTorS_" + tag, rBend, p.innerRadius, p.outerRadius, phiStart, phiRange);
  auto* volSteel = new TGeoVolume("VolTorS_" + tag, shSteel, p.medSteel);
  volSteel->SetLineColor(kGray + 2);

  auto* shWater = new TGeoTorus("shTorW_" + tag, rBend, 0, p.innerRadius, phiStart, phiRange);
  auto* volWater = new TGeoVolume("VolTorW_" + tag, shWater, p.medWater);
  volWater->SetLineColor(kBlue);

  assembly->AddNode(volSteel, idx, new TGeoTranslation(tx, ty, pipeZ));
  assembly->AddNode(volWater, idx, new TGeoTranslation(tx, ty, pipeZ));

  auto* shFit = new TGeoTorus("GrFit_" + tag, rBend, 0, p.outerRadius + fitClear, phiStart, phiRange);
  auto* trFit = new TGeoTranslation("GrFitT_" + tag, tx, ty, pipeZ);
  trFit->RegisterYourself();
  gActiveGrooveTerms->push_back(TString::Format("%s:%s", shFit->GetName(), trFit->GetName()));

  // Phi boundaries are widened by kAngEps to avoid a
  // coincident-surface artifact at the arc start/end.
  auto* shChan = new TGeoTubeSeg("GrChan_" + tag,
                                 rBend - chanWidth - fitClear,
                                 rBend + chanWidth + fitClear,
                                 halfFrameZ / 2. + NA6PPixelStation::kGeomEps,
                                 phiStart - NA6PPixelStation::kAngEps,
                                 phiStart + phiRange + NA6PPixelStation::kAngEps);
  auto* trChan = new TGeoTranslation("GrChanT_" + tag, tx, ty, pipeZ + halfFrameZ / 2.);
  trChan->RegisterYourself();
  gActiveGrooveTerms->push_back(TString::Format("%s:%s", shChan->GetName(), trChan->GetName()));
}

TGeoXtru* CreateRoundedRect(const char* name, double halfX, double halfY, double cornerRadius, double halfZ)
{
  const int nSegments = 20;
  const int nPoints = 4 * nSegments;
  double x[nPoints];
  double y[nPoints];
  for (int i = 0; i < nSegments; i++) {
    double phi = i * (TMath::Pi() / 2) / (nSegments - 1);
    x[i] = (halfX - cornerRadius) + cornerRadius * std::cos(phi);
    y[i] = (halfY - cornerRadius) + cornerRadius * std::sin(phi);
    x[i + nSegments] = -(halfX - cornerRadius) + cornerRadius * std::cos(phi + TMath::Pi() / 2);
    y[i + nSegments] = (halfY - cornerRadius) + cornerRadius * std::sin(phi + TMath::Pi() / 2);
    x[i + 2 * nSegments] = -(halfX - cornerRadius) + cornerRadius * std::cos(phi + TMath::Pi());
    y[i + 2 * nSegments] = -(halfY - cornerRadius) + cornerRadius * std::sin(phi + TMath::Pi());
    x[i + 3 * nSegments] = (halfX - cornerRadius) + cornerRadius * std::cos(phi + 3 * TMath::Pi() / 2);
    y[i + 3 * nSegments] = -(halfY - cornerRadius) + cornerRadius * std::sin(phi + 3 * TMath::Pi() / 2);
  }
  auto* xtru = new TGeoXtru(2);
  xtru->SetName(name);
  xtru->DefinePolygon(nPoints, x, y);
  xtru->DefineSection(0, -halfZ);
  xtru->DefineSection(1, halfZ);
  return xtru;
}

using PipeRoute = std::vector<NA6PPipeSegment>;

NA6PPipeSegment StraightSegment(Double_t halfLength, Double_t tx, Double_t ty, Double_t phi = 90.0, Double_t theta = 90.0)
{
  NA6PPipeSegment segment;
  segment.kind = NA6PPipeSegment::Kind::Straight;
  segment.halfLength = halfLength;
  segment.tx = tx;
  segment.ty = ty;
  segment.phi = phi;
  segment.theta = theta;
  return segment;
}

NA6PPipeSegment CurvedSegment(Double_t rBend, Double_t phiStart, Double_t phiRange, Double_t tx, Double_t ty)
{
  NA6PPipeSegment segment;
  segment.kind = NA6PPipeSegment::Kind::Curved;
  segment.rBend = rBend;
  segment.phiStart = phiStart;
  segment.phiRange = phiRange;
  segment.tx = tx;
  segment.ty = ty;
  return segment;
}

void AddPipeRoute(TGeoVolume* assembly,
                  const NA6PCoolingPipeParams& pipeP,
                  const PipeRoute& route,
                  const NA6PStationGeometryParams& geom,
                  Double_t chanWidth)
{
  for (const auto& segment : route) {
    if (segment.kind == NA6PPipeSegment::Kind::Straight) {
      AddStraightTube(assembly, pipeP, segment.halfLength, segment.tx, segment.ty,
                      segment.phi, segment.theta, geom.pipeZ, geom.halfFrameZ, geom.grooveFitClear);
    } else {
      AddCurvedTube(assembly, pipeP, segment.rBend, segment.phiStart, segment.phiRange,
                    segment.tx, segment.ty, geom.pipeZ, geom.halfFrameZ, chanWidth, geom.grooveFitClear);
    }
  }
}

TGeoVolume* BuildFirstStationsCoolingPlate(const char* tag,
                                           const NA6PCoolingPipeParams& pipeP,
                                           TGeoMedium* medPlate,
                                           TGeoMedium* medCarbonFiber,
                                           const NA6PStationGeometryParams& geom)
{
  gActiveGrooveTerms = &gGrooveTerms;
  gGrooveTerms.clear();
  gGrooveNamePrefix = TString::Format("Inner_%s", tag);

  const Double_t halfFrameX = geom.halfFrameX;
  const Double_t halfFrameY = geom.halfFrameY;
  const Double_t halfFrameZ = geom.halfFrameZ;
  const Double_t halfAlCutZ = geom.halfAlCutZ;
  const Double_t rBend = 1.985;
  const Double_t chanWidth = pipeP.outerRadius;
  const Double_t fiberCenterHoleR = geom.fiberCenterHoleR;

  const Double_t halfHor1 = 35.1 / 2.0;
  const Double_t halfVer1 = 24.2 / 2.0;
  const Double_t halfHor2 = 24.2 / 2.0;
  const Double_t halfVer2 = 15.7 / 2.0;
  const Double_t halfHor3 = 3.0 / 2.0;
  const Double_t halfVer3 = 8.7 / 2.0;
  const Double_t halfHor4 = 10.2 / 2.0;
  const Double_t halfVer4 = 10.2 / 2.0;
  const Double_t halfHor5 = 17.2 / 2.0;
  const Double_t halfVer5 = 1.8 / 2.0;
  const Double_t halfHor6 = 6.9 / 2.0;

  auto TX = [&](Double_t x0) { return x0 + geom.innerTubeOffsetX; };
  auto TY = [&](Double_t y0) { return y0 + geom.tubeOffsetY; };

  TString asmName = TString::Format("FrameAsm_%s", tag);
  TGeoVolume* assembly = new TGeoVolumeAssembly(asmName);
  assembly->SetVisibility(kFALSE);

  const PipeRoute route = {
    StraightSegment(halfHor1, TX(0), TY(0)),
    CurvedSegment(rBend, 90, 90, TX(-halfHor1), TY(-rBend)),
    StraightSegment(halfVer1, TX(-halfHor1 - rBend), TY(-rBend - halfVer1), 0),
    CurvedSegment(rBend, 180, 90, TX(-halfHor1), TY(-rBend - 2 * halfVer1)),
    StraightSegment(halfHor2, TX(-(halfHor1 - halfHor2)), TY(-rBend - 2 * halfVer1 - rBend)),
    CurvedSegment(rBend, 270, 90, TX(2 * halfHor2 - halfHor1), TY(-rBend - 2 * halfVer1)),
    StraightSegment(halfVer2, TX(2 * halfHor2 - halfHor1 + rBend), TY(-rBend - 2 * halfVer1 + halfVer2), 0),
    CurvedSegment(rBend, 0, 90, TX(2 * halfHor2 - halfHor1), TY(-rBend - 2 * halfVer1 + 2 * halfVer2)),
    StraightSegment(halfHor3, TX(2 * halfHor2 - halfHor1 - halfHor3), TY(-2 * halfVer1 + 2 * halfVer2)),
    CurvedSegment(rBend, 90, 90, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3), TY(-rBend - 2 * halfVer1 + 2 * halfVer2)),
    StraightSegment(halfVer3, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - rBend), TY(-rBend - 2 * halfVer1 + 2 * halfVer2 - halfVer3), 0),
    CurvedSegment(rBend, 270, 90, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend), TY(-rBend - 2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3)),
    StraightSegment(halfHor4, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend - halfHor4), TY(-2 * rBend - 2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3)),
    CurvedSegment(rBend, 180, 90, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend - 2 * halfHor4), TY(-rBend - 2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3)),
    StraightSegment(halfVer4, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 3 * rBend - 2 * halfHor4), TY(-rBend - 2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3 + halfVer4), 0),
    CurvedSegment(rBend, 90, 90, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend - 2 * halfHor4), TY(-rBend - 2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3 + 2 * halfVer4)),
    StraightSegment(halfHor5, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend - 2 * halfHor4 + halfHor5), TY(-2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3 + 2 * halfVer4)),
    CurvedSegment(rBend, 270, 90, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend - 2 * halfHor4 + 2 * halfHor5), TY(-2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3 + 2 * halfVer4 + rBend)),
    StraightSegment(halfVer5, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend - 2 * halfHor4 + 2 * halfHor5 + rBend), TY(-2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3 + 2 * halfVer4 + rBend + halfVer5), 0),
    CurvedSegment(rBend, 90, 90, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend - 2 * halfHor4 + 2 * halfHor5 + 2 * rBend), TY(-2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3 + 2 * halfVer4 + rBend + 2 * halfVer5)),
    StraightSegment(halfHor6, TX(2 * halfHor2 - halfHor1 - 2 * halfHor3 - 2 * rBend - 2 * halfHor4 + 2 * halfHor5 + 2 * rBend + halfHor6), TY(-2 * halfVer1 + 2 * halfVer2 - 2 * halfVer3 + 2 * halfVer4 + 2 * rBend + 2 * halfVer5))};
  AddPipeRoute(assembly, pipeP, route, geom, chanWidth);

  TString base = TString::Format("FrameBase_%s", tag);
  TString hCn = TString::Format("HoleC_%s", tag);
  TString hHn = TString::Format("HoleH_%s", tag);
  TString hVn = TString::Format("HoleV_%s", tag);
  TString fiberFront = TString::Format("CarbonFiberFront_%s", tag);
  TString fiberBack = TString::Format("CarbonFiberBack_%s", tag);
  TString fiberCenterHole = TString::Format("CarbonFiberCenterHole_%s", tag);
  TString fiberWindowHole = TString::Format("CarbonFiberWindowHole_%s", tag);
  TString fiberFrontShapeName = TString::Format("CarbonFiberFrontShape_%s", tag);
  TString fiberBackShapeName = TString::Format("CarbonFiberBackShape_%s", tag);
  TString fiberFrontVolName = TString::Format("CarbonFiberFrontVol_%s", tag);
  TString fiberBackVolName = TString::Format("CarbonFiberBackVol_%s", tag);
  TString tTopN = TString::Format("TTop_%s", tag);
  TString tBotN = TString::Format("TBot_%s", tag);
  TString tLN = TString::Format("TL_%s", tag);
  TString tRN = TString::Format("TR_%s", tag);
  TString tFiberFrontN = TString::Format("TFiberFront_%s", tag);
  TString tFiberBackN = TString::Format("TFiberBack_%s", tag);
  TString tFiberFrontTLN = TString::Format("TFiberFrontTL_%s", tag);
  TString tFiberFrontBRN = TString::Format("TFiberFrontBR_%s", tag);
  TString tFiberBackTRN = TString::Format("TFiberBackTR_%s", tag);
  TString tFiberBackBLN = TString::Format("TFiberBackBL_%s", tag);

  new TGeoBBox(base, halfFrameX, halfFrameY, halfFrameZ);
  const Double_t cHolePlateR = 1.5;
  const Double_t cHolePlateSideX = 9.8;
  const Double_t cHolePlateSideY = 10.;
  const Double_t halfCHolePlateX = cHolePlateSideX / 2 + cHolePlateR;
  const Double_t halfCHolePlateY = cHolePlateSideY / 2 + cHolePlateR;
  CreateRoundedRect(hCn, halfCHolePlateX, halfCHolePlateY, cHolePlateR, halfAlCutZ + 2 * NA6PPixelStation::kGeomEps);

  constexpr double alHSlotR = 2.5;
  constexpr double alHSlotL = 20.6;
  const Double_t halfAlHSlotX = alHSlotL / 2 + alHSlotR;
  const Double_t halfAlHSlotY = alHSlotR;
  CreateRoundedRect(hHn, halfAlHSlotX, halfAlHSlotY, alHSlotR, halfAlCutZ + 3 * NA6PPixelStation::kGeomEps);
  constexpr double alCornerR = 1.0;
  constexpr double alVSlotY = 9.;
  constexpr double alVSlotX = 3.;
  const Double_t halfAlVSlotX = alVSlotX / 2 + alCornerR;
  const Double_t halfAlVSlotY = alVSlotY / 2 + alCornerR;
  CreateRoundedRect(hVn, halfAlVSlotX, halfAlVSlotY, alCornerR, halfAlCutZ + 2 * NA6PPixelStation::kGeomEps);

  const Double_t halfFiberZ = 0.04 / 2;
  const Double_t fiberR = 0.;
  const Double_t fiberSideX = 36.;
  const Double_t fiberSideY = 17.;
  const Double_t halfFiberX = fiberSideX / 2 + fiberR;
  const Double_t halfFiberY = fiberSideY / 2 + fiberR;
  CreateRoundedRect(fiberFront, halfFiberX, halfFiberY, fiberR, halfFiberZ);
  CreateRoundedRect(fiberBack, halfFiberX, halfFiberY, fiberR, halfFiberZ);
  const Double_t fiberWindowR = 0.1;
  const Double_t fiberWindowSideX = 5.8;
  const Double_t fiberWindowSideY = 5.6;
  const Double_t halfFiberWindowX = fiberWindowSideX / 2 + fiberWindowR;
  const Double_t halfFiberWindowY = fiberWindowSideY / 2 + fiberWindowR;
  CreateRoundedRect(fiberWindowHole, halfFiberWindowX, halfFiberWindowY, fiberWindowR, halfFiberZ + NA6PPixelStation::kGeomEps);
  new TGeoTube(fiberCenterHole, 0., fiberCenterHoleR, halfFiberZ + NA6PPixelStation::kGeomEps);

  auto* tTop = new TGeoTranslation(tTopN, 0, 11., 0);
  tTop->RegisterYourself();
  auto* tBot = new TGeoTranslation(tBotN, 0, -11., 0);
  tBot->RegisterYourself();
  auto* tL = new TGeoTranslation(tLN, -10.3, -0.5, 0);
  tL->RegisterYourself();
  auto* tR = new TGeoTranslation(tRN, 10.3, -0.5, 0);
  tR->RegisterYourself();
  auto* tFiberFront = new TGeoTranslation(tFiberFrontN, 0, 0, -halfFrameZ - halfFiberZ);
  tFiberFront->RegisterYourself();
  auto* tFiberBack = new TGeoTranslation(tFiberBackN, 0, 0, halfFrameZ + halfFiberZ);
  tFiberBack->RegisterYourself();

  const Double_t fiberWindowX = fiberWindowSideX / 2 + fiberWindowR;
  const Double_t shiftFiberWindow = NA6PLayoutParam::Instance().shiftFiberWindow;
  const Double_t fiberWindowY = fiberWindowSideY / 2 + fiberWindowR;

  auto* tFiberFrontTL = new TGeoTranslation(tFiberFrontTLN, +fiberWindowX + shiftFiberWindow, +fiberWindowY - fiberCenterHoleR, 0);
  tFiberFrontTL->RegisterYourself();
  auto* tFiberFrontBR = new TGeoTranslation(tFiberFrontBRN, -fiberWindowX - shiftFiberWindow, -fiberWindowY + fiberCenterHoleR, 0);
  tFiberFrontBR->RegisterYourself();

  auto* tFiberBackTR = new TGeoTranslation(tFiberBackTRN, -fiberWindowX + shiftFiberWindow, +fiberWindowY + fiberCenterHoleR, 0);
  tFiberBackTR->RegisterYourself();
  auto* tFiberBackBL = new TGeoTranslation(tFiberBackBLN, +fiberWindowX - shiftFiberWindow, -fiberWindowY - fiberCenterHoleR, 0);
  tFiberBackBL->RegisterYourself();

  TString grooveExpr = "";
  for (const auto& term : gGrooveTerms) {
    grooveExpr += " + " + term;
  }

  TString fiberFrontExpr = fiberFront + " - (" + fiberCenterHole + " + " + fiberWindowHole + ":" + tFiberFrontTLN + " + " + fiberWindowHole + ":" + tFiberFrontBRN + ")";
  TString fiberBackExpr = fiberBack + " - (" + fiberCenterHole + " + " + fiberWindowHole + ":" + tFiberBackTRN + " + " + fiberWindowHole + ":" + tFiberBackBLN + ")";
  auto* carbonFiberFrontShape = new TGeoCompositeShape(fiberFrontShapeName, fiberFrontExpr.Data());
  auto* carbonFiberBackShape = new TGeoCompositeShape(fiberBackShapeName, fiberBackExpr.Data());
  auto* carbonFiberFront = new TGeoVolume(fiberFrontVolName, carbonFiberFrontShape, medCarbonFiber);
  carbonFiberFront->SetLineColor(kGray + 2);
  auto* carbonFiberBack = new TGeoVolume(fiberBackVolName, carbonFiberBackShape, medCarbonFiber);
  carbonFiberBack->SetLineColor(kGray + 2);

  TString alExpr = base + " - (" + hCn + " + " + hHn + ":" + tTopN + " + " + hHn + ":" + tBotN + " + " + hVn + ":" + tLN + " + " + hVn + ":" + tRN + grooveExpr + ")";

  TString shapeName = TString::Format("FrameShape_%s", tag);
  TString volName = TString::Format("Frame_%s", tag);

  auto* frameShape = new TGeoCompositeShape(shapeName, alExpr.Data());
  auto* frame = new TGeoVolume(volName, frameShape, medPlate);
  frame->SetLineColor(kAzure - 9);

  assembly->AddNode(frame, 1, new TGeoTranslation(0, 0, 0));
  assembly->AddNode(carbonFiberFront, 1, new TGeoTranslation(0, 0, -halfFrameZ - halfFiberZ));
  assembly->AddNode(carbonFiberBack, 1, new TGeoTranslation(0, 0, +halfFrameZ + halfFiberZ));
  return assembly;
}

TGeoVolume* BuildOuterStationsCoolingPlate(const char* tag,
                                           const NA6PCoolingPipeParams& pipeP,
                                           TGeoMedium* medPlate,
                                           TGeoMedium* medCarbonFoam,
                                           TGeoMedium* medCarbonFiber,
                                           TGeoMedium* medAir,
                                           const NA6PStationGeometryParams& geom)
{
  gActiveGrooveTerms = &gGrooveTermsOuter;
  gGrooveTermsOuter.clear();
  gGrooveNamePrefix = TString::Format("Outer_%s", tag);

  // Placeholder for the more complex plate used by the last three stations.
  // Keep it routed separately so its geometry can diverge without touching
  // the station assembly logic.
  const Double_t halfFrameX = geom.halfFrameX;
  const Double_t halfFrameY = geom.halfFrameY;
  const Double_t halfFrameZ = geom.halfFrameZ;
  const Double_t halfAlCutZ = geom.halfAlCutZ;
  const Double_t rBend = 1.985;
  const Double_t rBend2 = 1.485;
  const Double_t chanWidth = pipeP.outerRadius;
  const Double_t fiberCenterHoleR = geom.fiberCenterHoleR;

  const Double_t halfHor1 = 34.6 / 2.0;
  const Double_t halfVer1 = 24.2 / 2.0;
  const Double_t halfHor2 = 23.2 / 2.0;
  const Double_t halfVer2 = 19.7 / 2.0;
  const Double_t halfVer3 = 16.2 / 2.0;
  const Double_t halfHor3 = 17.2 / 2.0;
  const Double_t halfVer4 = 17.2 / 2.0;
  const Double_t halfHor4 = 21.236 / 2.0;
  const Double_t halfHor5 = 6.9 / 2.0;

  auto TX = [&](Double_t x0) { return x0 + geom.outerTubeOffsetX; };
  auto TY = [&](Double_t y0) { return y0 + geom.tubeOffsetY; };

  TString asmName = TString::Format("FrameAsm_%s", tag);
  TGeoVolume* assembly = new TGeoVolumeAssembly(asmName);
  assembly->SetVisibility(kFALSE);

  PipeRoute route = {
    StraightSegment(halfHor1, TX(0), TY(0)),
    CurvedSegment(rBend, 90, 90, TX(-halfHor1), TY(-rBend)),
    StraightSegment(halfVer1, TX(-halfHor1 - rBend), TY(-rBend - halfVer1), 0),
    CurvedSegment(rBend, 180, 90, TX(-halfHor1), TY(-rBend - 2 * halfVer1)),
    StraightSegment(halfHor2, TX(-(halfHor1 - halfHor2)), TY(-rBend - 2 * halfVer1 - rBend)),
    CurvedSegment(rBend, 270, 90, TX(2 * halfHor2 - halfHor1), TY(-rBend - 2 * halfVer1)),
    StraightSegment(halfVer2, TX(2 * halfHor2 - halfHor1 + rBend), TY(-rBend - 2 * halfVer1 + halfVer2), 0)};
  const Double_t bend2CenterX = 2 * halfHor2 - halfHor1 + rBend - rBend2;
  const Double_t bend2CenterY = -rBend - 2 * halfVer1 + 2 * halfVer2;
  const Double_t bend2SecondCenterX = bend2CenterX;
  const Double_t bend2SecondCenterY = bend2CenterY;
  const Double_t afterBend2X = bend2SecondCenterX - rBend2;
  const Double_t afterBend2Y = bend2SecondCenterY;
  route.push_back(CurvedSegment(rBend2, 0, 90, TX(bend2CenterX), TY(bend2CenterY)));
  route.push_back(CurvedSegment(rBend2, 90, 90, TX(bend2SecondCenterX), TY(bend2SecondCenterY)));
  route.push_back(StraightSegment(halfVer3, TX(afterBend2X), TY(afterBend2Y - halfVer3), 0));
  route.push_back(CurvedSegment(rBend, 270, 90, TX(afterBend2X - rBend), TY(afterBend2Y - 2 * halfVer3)));
  route.push_back(StraightSegment(halfHor3, TX(afterBend2X - rBend - halfHor3), TY(afterBend2Y - 2 * halfVer3 - rBend)));
  route.push_back(CurvedSegment(rBend, 180, 90, TX(afterBend2X - rBend - 2 * halfHor3), TY(afterBend2Y - 2 * halfVer3)));
  route.push_back(StraightSegment(halfVer4, TX(afterBend2X - 2 * rBend - 2 * halfHor3), TY(afterBend2Y - 2 * halfVer3 + halfVer4), 0));
  route.push_back(CurvedSegment(rBend, 90, 90, TX(afterBend2X - rBend - 2 * halfHor3), TY(afterBend2Y - 2 * halfVer3 + 2 * halfVer4)));
  route.push_back(StraightSegment(halfHor4, TX(afterBend2X - rBend - 2 * halfHor3 + halfHor4), TY(afterBend2Y - 2 * halfVer3 + 2 * halfVer4 + rBend)));

  const Double_t angle = geom.outerOutletAngleDeg;
  const Double_t arcStartX = afterBend2X - rBend - 2 * halfHor3 + 2 * halfHor4;
  const Double_t arcStartY = afterBend2Y - 2 * halfVer3 + 2 * halfVer4 + rBend;
  const Double_t arcCenterX = arcStartX;
  const Double_t arcCenterY = arcStartY + rBend;
  const Double_t arcCenterStepX = 2 * rBend * TMath::Sin(angle * TMath::DegToRad());
  const Double_t arcCenterStepY = 2 * rBend * TMath::Cos(angle * TMath::DegToRad());
  route.push_back(CurvedSegment(rBend, 270, angle, TX(arcCenterX), TY(arcCenterY)));
  route.push_back(CurvedSegment(rBend, 90, angle, TX(arcCenterX + arcCenterStepX), TY(arcCenterY - arcCenterStepY)));
  const Double_t outletStartX = arcCenterX + arcCenterStepX;
  const Double_t outletHalfLength = 0.5 * (halfHor1 - outletStartX);
  route.push_back(StraightSegment(outletHalfLength,
                                  TX(outletStartX + outletHalfLength),
                                  TY(arcCenterY - arcCenterStepY + rBend)));
  AddPipeRoute(assembly, pipeP, route, geom, chanWidth);

  TString base = TString::Format("FrameBase_%s", tag);
  TString hCn = TString::Format("HoleCDeep_%s", tag);
  TString hCnF = TString::Format("HoleCShallowFront_%s", tag);
  TString hCnB = TString::Format("HoleCShallowBack_%s", tag);
  TString foamBase = TString::Format("CarbonFoamBase_%s", tag);
  TString foamDiskHole = TString::Format("CarbonFoamDiskHole_%s", tag);
  TString foamShapeName = TString::Format("CarbonFoamNoGrooveShape_%s", tag);
  TString foamVolName = TString::Format("CarbonFoamNoGroove_%s", tag);
  TString fiberCenterHole = TString::Format("CarbonFiberCenterHole_%s", tag);
  TString fiberWindowHole = TString::Format("CarbonFiberWindowHole_%s", tag);
  TString fiberFrontShapeName = TString::Format("CarbonFiberFrontWithHoleShape_%s", tag);
  TString fiberBackShapeName = TString::Format("CarbonFiberBackWithHoleShape_%s", tag);
  TString fiberFrontVolName = TString::Format("CarbonFiberFrontVol_%s", tag);
  TString fiberBackVolName = TString::Format("CarbonFiberBackVol_%s", tag);
  TString hHn = TString::Format("HoleH_%s", tag);
  TString hVn = TString::Format("HoleV_%s", tag);
  TString tFrontN = TString::Format("TFront_%s", tag);
  TString tBackN = TString::Format("TBack_%s", tag);
  TString tFoamN = TString::Format("TFoam_%s", tag);

  TString fiberFront = TString::Format("CarbonFiberFront_%s", tag);
  TString fiberBack = TString::Format("CarbonFiberBack_%s", tag);
  TString tFiberFrontN = TString::Format("TFiberFront_%s", tag);
  TString tFiberBackN = TString::Format("TFiberBack_%s", tag);
  TString tFiberFrontTLN = TString::Format("TFiberFrontTL_%s", tag);
  TString tFiberFrontTRN = TString::Format("TFiberFrontTR_%s", tag);
  TString tFiberFrontBLN = TString::Format("TFiberFrontBL_%s", tag);
  TString tFiberFrontBRN = TString::Format("TFiberFrontBR_%s", tag);
  TString tFiberBackTLN = TString::Format("TFiberBackTL_%s", tag);
  TString tFiberBackTRN = TString::Format("TFiberBackTR_%s", tag);
  TString tFiberBackBLN = TString::Format("TFiberBackBL_%s", tag);
  TString tFiberBackBRN = TString::Format("TFiberBackBR_%s", tag);

  new TGeoBBox(base, halfFrameX, halfFrameY, halfFrameZ);
  // Create the holes in the aluminum frame
  const Double_t cHoleR = 1.5;
  const Double_t cHoleSideX = 23;
  const Double_t cHoleSideY = 24.2;
  const Double_t halfCHoleX = cHoleSideX / 2 + cHoleR;
  const Double_t halfCHoleY = cHoleSideY / 2 + cHoleR;
  const Double_t halfCHoleCutZ = halfFrameZ + 2 * NA6PPixelStation::kGeomEps;
  CreateRoundedRect(hCn, halfCHoleX, halfCHoleY, cHoleR, halfCHoleCutZ);
  const Double_t fbHoleR = 0.5;
  const Double_t fbHoleSideX = 31.2;
  const Double_t fbHoleSideY = 28.2;
  const Double_t halfFbHoleX = fbHoleSideX / 2 + fbHoleR;
  const Double_t halfFbHoleY = fbHoleSideY / 2 + fbHoleR;
  const Double_t halfFrontHoleZ = geom.outerSensorPocketDz / 2;
  const Double_t halfBackHoleZ = geom.outerSensorPocketDz / 2;
  const Double_t halfFrontHoleCutZ = halfFrontHoleZ + 4 * NA6PPixelStation::kGeomEps;
  const Double_t halfBackHoleCutZ = halfBackHoleZ + 4 * NA6PPixelStation::kGeomEps;
  CreateRoundedRect(hCnF, halfFbHoleX, halfFbHoleY, fbHoleR, halfFrontHoleCutZ);
  CreateRoundedRect(hCnB, halfFbHoleX, halfFbHoleY, fbHoleR, halfBackHoleCutZ);

  // Create the carbon foam
  const Double_t halfFoamZ = geom.outerCarbonFoamHalfZ;
  const Double_t foamR = 2;
  const Double_t foamSideX = 21.8;
  const Double_t foamSideY = 23.;
  const Double_t halfFoamX = foamSideX / 2 + foamR;
  const Double_t halfFoamY = foamSideY / 2 + foamR;
  CreateRoundedRect(foamBase, halfFoamX, halfFoamY, foamR, halfFoamZ);

  // Create the carbon fiber and the holes in it
  const Double_t halfFiberZ = geom.outerCarbonFiberHalfZ;
  const Double_t fiberR = 0.5;
  const Double_t fiberSideX = 31.2;
  const Double_t fiberSideY = 28.2;
  const Double_t halfFiberX = fiberSideX / 2 + fiberR;
  const Double_t halfFiberY = fiberSideY / 2 + fiberR;

  CreateRoundedRect(fiberFront, halfFiberX, halfFiberY, fiberR, halfFiberZ);
  CreateRoundedRect(fiberBack, halfFiberX, halfFiberY, fiberR, halfFiberZ);

  const Double_t fiberWindowR = 0.5;
  const Double_t fiberWindowSide = 10;
  const Double_t halfFiberWindow = fiberWindowSide / 2 + fiberWindowR;
  CreateRoundedRect(fiberWindowHole, halfFiberWindow, halfFiberWindow, fiberWindowR, halfFiberZ + NA6PPixelStation::kGeomEps);

  const Double_t foamHoleR = 10.;
  new TGeoTube(foamDiskHole, 0., foamHoleR, halfFrameZ + 2 * NA6PPixelStation::kGeomEps);
  new TGeoTube(fiberCenterHole, 0., fiberCenterHoleR, halfFiberZ + NA6PPixelStation::kGeomEps);

  auto* tFront = new TGeoTranslation(tFrontN, 0, 0, -halfFrameZ + halfFrontHoleZ);
  tFront->RegisterYourself();
  auto* tBack = new TGeoTranslation(tBackN, 0, 0, halfFrameZ - halfBackHoleZ);
  tBack->RegisterYourself();
  auto* tFoam = new TGeoTranslation(tFoamN, 0, 0, 0);
  tFoam->RegisterYourself();
  auto* tFiberFront = new TGeoTranslation(tFiberFrontN, 0, 0, -halfFoamZ - halfFiberZ);
  tFiberFront->RegisterYourself();
  auto* tFiberBack = new TGeoTranslation(tFiberBackN, 0, 0, halfFoamZ + halfFiberZ);
  tFiberBack->RegisterYourself();

  const Double_t fiberWindowX = fiberWindowSide / 2 + fiberWindowR;
  const Double_t shiftFiberWindow = NA6PLayoutParam::Instance().shiftFiberWindow;
  const Double_t fiberWindowY = fiberWindowSide / 2 + fiberWindowR;

  auto* tFiberFrontTL = new TGeoTranslation(tFiberFrontTLN, +fiberWindowX + shiftFiberWindow, +fiberWindowY - fiberCenterHoleR, 0);
  tFiberFrontTL->RegisterYourself();
  auto* tFiberFrontBR = new TGeoTranslation(tFiberFrontBRN, -fiberWindowX - shiftFiberWindow, -fiberWindowY + fiberCenterHoleR, 0);
  tFiberFrontBR->RegisterYourself();

  auto* tFiberBackTR = new TGeoTranslation(tFiberBackTRN, -fiberWindowX + shiftFiberWindow, +fiberWindowY + fiberCenterHoleR, 0);
  tFiberBackTR->RegisterYourself();
  auto* tFiberBackBL = new TGeoTranslation(tFiberBackBLN, +fiberWindowX - shiftFiberWindow, -fiberWindowY - fiberCenterHoleR, 0);
  tFiberBackBL->RegisterYourself();

  TString grooveExpr = "";
  TString alGrooveExpr = "";
  for (const auto& term : gGrooveTermsOuter) {
    grooveExpr += " + " + term;
    alGrooveExpr += " - " + term;
  }

  TString foamExpr = foamBase + " - (" + foamDiskHole + grooveExpr + ")";
  auto* carbonFoamShape = new TGeoCompositeShape(foamShapeName, foamExpr.Data());
  auto* carbonFoam = new TGeoVolume(foamVolName, carbonFoamShape, medCarbonFoam);
  carbonFoam->SetLineColor(kBlue - 6);

  TString fiberFrontExpr = fiberFront + " - (" + fiberCenterHole + " + " + fiberWindowHole + ":" + tFiberFrontTLN + " + " + fiberWindowHole + ":" + tFiberFrontBRN + ")";
  TString fiberBackExpr = fiberBack + " - (" + fiberCenterHole + " + " + fiberWindowHole + ":" + tFiberBackTRN + " + " + fiberWindowHole + ":" + tFiberBackBLN + ")";
  auto* carbonFiberBackShape = new TGeoCompositeShape(fiberBackShapeName, fiberBackExpr.Data());
  auto* carbonFiberFrontShape = new TGeoCompositeShape(fiberFrontShapeName, fiberFrontExpr.Data());
  auto* carbonFiberBack = new TGeoVolume(fiberBackVolName, carbonFiberBackShape, medCarbonFiber);
  carbonFiberBack->SetLineColor(kGray + 2);
  auto* carbonFiberFront = new TGeoVolume(fiberFrontVolName, carbonFiberFrontShape, medCarbonFiber);
  carbonFiberFront->SetLineColor(kGray + 2);

  TString alExpr = base + " - " + hCn + " - " + hCnF + ":" + tFrontN + " - " + hCnB + ":" + tBackN +
                   " - " + foamBase + ":" + tFoamN +
                   " - " + fiberFront + ":" + tFiberFrontN +
                   " - " + fiberBack + ":" + tFiberBackN +
                   alGrooveExpr;

  TString shapeName = TString::Format("FrameShape_%s", tag);
  TString volName = TString::Format("Frame_%s", tag);

  auto* frameShape = new TGeoCompositeShape(shapeName, alExpr.Data());
  auto* frame = new TGeoVolume(volName, frameShape, medPlate);
  frame->SetLineColor(kAzure - 9);

  assembly->AddNode(frame, 1, new TGeoTranslation(0, 0, 0));
  assembly->AddNode(carbonFoam, 1, new TGeoTranslation(0, 0, 0));
  assembly->AddNode(carbonFiberFront, 1, new TGeoTranslation(0, 0, -halfFoamZ - halfFiberZ));
  assembly->AddNode(carbonFiberBack, 1, new TGeoTranslation(0, 0, halfFoamZ + halfFiberZ));

  return assembly;
}

TGeoVolume* BuildCoolingPlateImpl(NA6PPixelStation::CoolingPlateType type,
                                  const char* tag,
                                  const NA6PCoolingPipeParams& pipeP,
                                  const NA6PCoolingPlateMaterials& materials,
                                  const NA6PStationGeometryParams& geom)
{
  switch (type) {
    case NA6PPixelStation::CoolingPlateType::FirstStations:
      return BuildFirstStationsCoolingPlate(tag,
                                            pipeP,
                                            materials.plate,
                                            materials.carbonFiber,
                                            geom);
    case NA6PPixelStation::CoolingPlateType::OuterStations:
      return BuildOuterStationsCoolingPlate(tag,
                                            pipeP,
                                            materials.plate,
                                            materials.carbonFoam,
                                            materials.carbonFiber,
                                            materials.air,
                                            geom);
  }
  return nullptr;
}
} // namespace

TGeoVolume* BuildNA6PCoolingPlate(NA6PPixelStation::CoolingPlateType type,
                                  const char* tag,
                                  const NA6PCoolingPipeParams& pipeP,
                                  const NA6PCoolingPlateMaterials& materials,
                                  const NA6PStationGeometryParams& geom)
{
  return BuildCoolingPlateImpl(type, tag, pipeP, materials, geom);
}
