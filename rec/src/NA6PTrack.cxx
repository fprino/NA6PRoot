// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include <TGeoGlobalMagField.h>
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PTrack.h"
#include "NA6PLayoutParam.h"

ClassImp(NA6PTrack)

  //_______________________________________________________________________
  NA6PTrack::NA6PTrack() : mMass{0.140},
                           mMatchChi2{0.0},
                           mChi2VT{0.0},
                           mChi2MS{0.0},
                           mClusterMap{0},
                           mNClusters{0},
                           mNClustersVT{0},
                           mNClustersMS{0},
                           mNClustersTR{0},
                           mClusterIndices{},
                           mClusterPartID{},
                           mParticleID{-2},
                           mExtTrack{}
{
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
}

//_______________________________________________________________________
NA6PTrack::NA6PTrack(const double* xyz, const double* pxyz, int sign, double errLoose) : mMass{0.140},
                                                                                         mMatchChi2{0.0},
                                                                                         mChi2VT{0.0},
                                                                                         mChi2MS{0.0},
                                                                                         mClusterMap{0},
                                                                                         mNClusters{0},
                                                                                         mNClustersVT{0},
                                                                                         mNClustersMS{0},
                                                                                         mNClustersTR{0},
                                                                                         mClusterIndices{},
                                                                                         mClusterPartID{},
                                                                                         mParticleID{-2},
                                                                                         mExtTrack{}
{
  // initialize arrays
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
  // initialize the track parameters and covariance
  init(xyz, pxyz, sign, errLoose);
}

//_______________________________________________________________________
void NA6PTrack::reset()
{
  mMass = 0.14;
  mMatchChi2 = 0;
  mChi2VT = 0;
  mChi2MS = 0;
  mChi2VTOuter = 0;
  mChi2MSOuter = 0;
  mChi2VTRefit = 0;
  mChi2MSRefit = 0;
  mClusterMap = 0;
  mExtTrack.Reset();
  mOuter.Reset();
  mConstrained.Reset();
  mStatusRefitInward = false;
  mStatusConstrained = false;
  resetCovariance();
  mParticleID = -2;
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
  mNClusters = mNClusters = mNClustersMS = mNClustersTR = 0;
}

//_______________________________________________________________________
Bool_t NA6PTrack::init(const double* xyz, const double* pxyz, int sign, double errLoose)
{
  // Init with track position/momentum in usual Lab frame
  // If errLoose>0 then scale initially small errors by this amount
  double xyzL[3], pxyzL[3];
  lab2trk(xyz, xyzL);
  lab2trk(pxyz, pxyzL);
  double cov[21] = {1.e-6, // assign small errors first
                    0., 1.e-6,
                    0., 0., 1.e-6,
                    0., 0., 0., 1.e-4,
                    0., 0., 0., 0., 1.e-4,
                    0., 0., 0., 0., 0., 1.e-3};

  mExtTrack.set(xyzL, pxyzL, cov, sign);
  if (errLoose > 0)
    mExtTrack.ResetCovariance(errLoose);
  mExtTrack.Rotate(pxyz[2] > 0 ? 0. : TMath::Pi());
  return true;
}

std::array<float, 3> NA6PTrack::getCircleParams(float by) const
{
  // get circle params as {R, Xcenter, Zcenter}, for straight line just set to current coordinates
  constexpr double ZeroCurv = 1e-8;
  float crv = mExtTrack.getC(by);
  if (std::abs(crv) > ZeroCurv) {
    float crvi = 1.0 / crv;
    float sn = mExtTrack.getSnp();
    float cs = std::sqrt((1. - sn) * (1. + sn));
    if (mExtTrack.Px() < 0) cs = -cs;
    float xcLab = mExtTrack.getY() + cs * crvi;
    float zcLab = mExtTrack.getX() - sn * crvi;
    return {std::abs(crvi), xcLab, zcLab};
  } else {
    return {0., static_cast<float>(mExtTrack.getY()),static_cast<float>( mExtTrack.getX())};
  }
}
//_______________________________________________________________________
Bool_t NA6PTrack::propagateToZBxByBz(double z, double maxDZ, double xOverX0, double xTimesRho, bool outer, double density, double atomicZ, double zOverA)
{
  // propagate the track to position Z in uniform material with xOverX0 rad lgt and xTimesRho lgt*density
  double zCurr = outer ? getZLabOuter() : getZLab();
  double dz = z - zCurr;

  if (TMath::Abs(dz) < kAlmost0) {
    return true;
  }

  int nz = TMath::Abs(dz) / maxDZ + 1;
  double zstep = dz / nz;
  double xyz[3], bxyz[3], bxyzFwd[3];

  for (int iz = 0; iz < nz; iz++) {
    if (outer) {
      getXYZOuter(xyz); // coordinates in Lab frame
    } else {
      getXYZ(xyz); // coordinates in Lab frame
    }
    TGeoGlobalMagField::Instance()->Field(xyz, bxyz);
    lab2trk(bxyz, bxyzFwd); // field in Fwd frame

    zCurr += zstep;

    Bool_t success;
    if (outer) {
      success = propagateToZBxByBzOuter(zCurr, bxyzFwd);
    } else {
      success = propagateToZBxByBz(zCurr, bxyzFwd);
    }

    if (!success) {
      return false;
    }

    if (TMath::Abs(xTimesRho) > 1e-6 && !correctForMeanMaterial(xOverX0 / nz, xTimesRho / nz, density, atomicZ, zOverA)) {
      return false;
    }
  }
  return true;
}

//_______________________________________________________________________
Bool_t NA6PTrack::propagateToDCA(NA6PTrack* partner)
{
  // propagate the track to position Z of closest approach to partner track
  //
  double xyz[3], bxyz[3], bxyzFwd[3];
  getXYZ(xyz); // coordinates in Lab frame
  TGeoGlobalMagField::Instance()->Field(xyz, bxyz);
  lab2trk(bxyz, bxyzFwd); // field in Fwd frame
  double zthis = 0, zpartner = 0;
  double dca = mExtTrack.getDCA(&partner->mExtTrack, bxyzFwd[2], zthis, zpartner);

  if (!propagateToZBxByBz(zthis) || !partner->propagateToZBxByBz(zpartner))
    return false;
  return true;
  //
}

bool NA6PTrack::propagateToDCABeamAxis(float beamX, float beamY, float maxDCA)
{
  // iterative propagation to DCA with beam axis (beamX, beamY, *) in z_lab steps 
  
  double xyz[3], pxyz[3];
  getXYZ(xyz);
  getPXYZ(pxyz);

  const int maxIter = 10;
  for (int iter = 0; iter < maxIter; ++iter) {
    double dx  = xyz[0] - beamX;
    double dy  = xyz[1] - beamY;
    double px  = pxyz[0];
    double py  = pxyz[1];
    double pz  = pxyz[2];
    double pt2 = px*px + py*py;
    if (pt2 < 1e-12) break;  // track parallel to beam, can't converge

    // derivative of transverse distance w.r.t. z_lab:
    // d/dz [(x-bx)^2 + (y-by)^2] = 2(dx*px/pz + dy*py/pz)
    double ddistdz = dx * (px / pz) + dy * (py / pz);
    // second derivative (approximated as pt^2/pz^2):
    double d2distdz2 = pt2 / (pz*pz);

    double dz = -ddistdz / d2distdz2;  // Newton step in z_lab

    // dampen large steps
    if (std::abs(dz) > 5.0) dz = std::copysign(5.0, dz);
    if (std::abs(dz) < 1e-4) break;  // converged

    double zNew = xyz[2] + dz;
    if (!propagateToZBxByBz(zNew)) return false;

    getXYZ(xyz);
    getPXYZ(pxyz);
  }

  // final DCA check
  double dx = xyz[0] - beamX, dy = xyz[1] - beamY;
  return (dx*dx + dy*dy) < maxDCA * maxDCA;
}

//_______________________________________________________________________
double NA6PTrack::getSigmaP2() const
{
  // error^2 on total momentum, P = sqrt(1+tgl^2)/(1/pt)
  double pinv = mExtTrack.getSigned1Pt();
  if (TMath::Abs(pinv) < kAlmost0)
    return 0;
  double tgl = mExtTrack.getTgl();
  double tglE = mExtTrack.getSigmaTgl2();
  double pinvE = mExtTrack.getSigma1Pt2();
  double pinvtgE = mExtTrack.getSigma1PtTgl();
  //
  double tp12 = TMath::Sqrt(1. + tgl * tgl);
  double dt = tgl / pinv / tp12;   // dP/dtgl
  double dc = -tp12 / pinv / pinv; // dP/dC
  double err2 = dt * dt * tglE + dc * dc * pinvE + dt * dc * pinvtgE;
  return err2;
  //
}

//_______________________________________________________________________
double NA6PTrack::getSigmaPX2() const
{
  // error^2 on Px in Lav frame (Z in the tracking frame, = tgl/(1/pt)
  double pinv = mExtTrack.getSigned1Pt();
  if (TMath::Abs(pinv) < kAlmost0)
    return 0;
  double tgl = mExtTrack.getTgl();
  double tglE = mExtTrack.getSigmaTgl2();
  double pinvE = mExtTrack.getSigma1Pt2();
  double pinvtgE = mExtTrack.getSigma1PtTgl();
  //
  double dt = 1. / pinv;          // dP/dtgl
  double dc = -tgl / pinv / pinv; // dP/dC
  double err2 = dt * dt * tglE + dc * dc * pinvE + dt * dc * pinvtgE;
  return err2;
  //
}

//_______________________________________________________________________
double NA6PTrack::getSigmaPY2() const
{
  // error^2 on Py in Lab frame (Y in the tracking frame, = sinp/(1/pt)
  double pinv = mExtTrack.getSigned1Pt();
  if (TMath::Abs(pinv) < kAlmost0)
    return 0;
  double cosAlp = TMath::Cos(mExtTrack.getAlpha()), sinAlp = TMath::Sin(mExtTrack.getAlpha());
  double snp = mExtTrack.getSnp();
  double csp = TMath::Sqrt((1. - snp) * (1. + snp));
  double snpE = mExtTrack.getSigmaSnp2();
  double pinvE = mExtTrack.getSigma1Pt2();
  double pinvsnpE = mExtTrack.getSigma1PtSnp();
  //
  double ds = (cosAlp - sinAlp * snp / csp) / pinv;         // dP/dsnp
  double dc = -(snp * cosAlp + csp * sinAlp) / pinv / pinv; // dP/dC
  double err2 = ds * ds * snpE + dc * dc * pinvE + ds * dc * pinvsnpE;
  return err2;

  //
}

//_______________________________________________________________________
double NA6PTrack::getSigmaPZ2() const
{
  // error^2 on Pz in Lab frame (X in the tracking frame, = (sqrt(1-snp^2)*cosAlp-snp*sinAlp)/(1/pt)
  double pinv = mExtTrack.getSigned1Pt();
  if (TMath::Abs(pinv) < kAlmost0)
    return 0;
  double cosAlp = TMath::Cos(mExtTrack.getAlpha()), sinAlp = TMath::Sin(mExtTrack.getAlpha());
  double snp = mExtTrack.getSnp();
  double snpE = mExtTrack.getSigmaSnp2();
  double pinvE = mExtTrack.getSigma1Pt2();
  double pinvsnpE = mExtTrack.getSigma1PtSnp();
  //
  double csp = TMath::Sqrt((1. - snp) * (1. + snp));
  double ds = -snp / pinv * (cosAlp / csp + sinAlp);        // dP/dsnp
  double dc = -(csp * cosAlp - snp * sinAlp) / pinv / pinv; // dP/dC
  double err2 = ds * ds * snpE + dc * dc * pinvE + ds * dc * pinvsnpE;
  return err2;
  //
}

//_______________________________________________________________________
void NA6PTrack::resetCovariance(float err)
{
  // reset cov matrix
  double* trCov = mExtTrack.getCovarianceForUpdate();
  const double* trPars = mExtTrack.getParameter();
  const double kLargeErr2Coord = 50 * 50;
  const double kLargeErr2Dir = 0.6 * 0.6;
  const double kLargeErr2PtI = 2.;
  for (int ic = 15; ic--;)
    trCov[ic] = 0.;
  trCov[kY2] = trCov[kZ2] = err < 0 ? kLargeErr2Coord : err * err;
  trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
  trCov[kPtI2] = kLargeErr2PtI * trPars[kPtI] * trPars[kPtI];
  mExtTrack.checkCovariance();
}

//_______________________________________________________________________
std::string NA6PTrack::asString() const
{
  double pxyz[3];
  getPXYZ(pxyz);
  return fmt::format("Track: Nclusters:{} NVTclusters:{} NMSclusters:{} NTRclusters:{} match chi2:{} chi2VT:{} chi2MS:{} pos:{:.4f},{:.4f},{:.4f} mom:{:.3f},{:.3f},{:.3f}",
                     mNClusters, mNClustersVT, mNClustersMS, mNClustersTR, mMatchChi2, mChi2VT, mChi2MS, getXLab(), getYLab(), getZLab(), pxyz[0], pxyz[1], pxyz[2]);
}

//_______________________________________________________________________
void NA6PTrack::print() const
{
  LOGP(info, "{}", asString());
}

//_______________________________________
template <typename ClusterType>
void NA6PTrack::addCluster(const ClusterType* clu, int cluIndex, double chi2)
{

  mNClusters++;
  int trackID = clu->getParticleID();

  int nLay = clu->getLayer();
  mClusterPartID[nLay] = trackID;
  mClusterIndices[nLay] = cluIndex;
  mClusterMap |= (1 << nLay);

  if (nLay < NA6PLayoutParam::Instance().nVerTelPlanes) {
    mNClustersVT++;
    mChi2VT += chi2;
  } else if (nLay >= NA6PLayoutParam::Instance().nVerTelPlanes + NA6PLayoutParam::Instance().nMSPlanes - 2) {
    mNClustersTR++;
    mChi2MS += chi2;
  } else {
    mNClustersMS++;
    mChi2MS += chi2;
  }
}

template void NA6PTrack::addCluster<NA6PBaseCluster>(const NA6PBaseCluster*, int, double);
template void NA6PTrack::addCluster<NA6PVerTelCluster>(const NA6PVerTelCluster*, int, double);
template void NA6PTrack::addCluster<NA6PMuonSpecCluster>(const NA6PMuonSpecCluster*, int, double);

//____________________________________________
void NA6PTrack::imposeKinematics(const double* xyzLab, const double* cosinesLab,
                                 double en, double mass, int charge)
{
  // RS: note: we assume p=e here to avoid problem with e+ e- interpreted as muon
  double p = en * en - mass * mass;
  if (p <= 0) {
    printf("Anomalous kinematics: E:%e M:%e", en, mass);
    exit(1);
  }
  p = TMath::Sqrt(p);
  double pxyz[3] = {p * cosinesLab[0], p * cosinesLab[1], p * cosinesLab[2]};
  init(xyzLab, pxyz, charge, 1.e4);
  setMass(mass);
  resetCovariance(); // reset cov.matrix
}
