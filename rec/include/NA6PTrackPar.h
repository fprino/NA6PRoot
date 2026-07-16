// NA6PCCopyright

#ifndef NA6P_TRACKPAR_XZ_H
#define NA6P_TRACKPAR_XZ_H

#include <Rtypes.h>
#include <array>
#include <cmath>
#include <string>
#include <cfloat>
#include <algorithm>
#include "PID.h"

// Units: cm, kGauss, GeV.
// Optimized for uniform dipole field B = (0, By, 0).
//
// State a = { x, y, tx, ty, q/pxz } evaluated at z,
//   tx = px/pxz = sin(psi)   (psi is momentum angle in x–z plane, pz>0 => cos(psi)>0)
//   ty = py/pxz
//   q/pxz = charge / |p_xz|
//
// Invariants in pure By magnetic field:
//   pxz = const,  py = const,  q/pxz = const,  ty = py/pxz = const
//
// Curvature in x–z plane:
//   kappa = kB2C * By * (q/pxz)
//
// Exact z-geometry gives:
//   sin(psi1) - sin(psi0) = kappa * dz  =>  tx1 = tx0 + kappa*dz  (must keep |tx|<=1)
//   cos(psi) = sqrt(1 - tx^2)  (positive since pz>0)
//
// Exact x update (same circular-arc identity):
//   x1 = x0 + (cos(psi0) - cos(psi1)) / kappa     (if |kappa|>0)
//
// y motion: drift along y with dy/ds_perp = ty, where s_perp is arc length in x–z plane
//   s_perp = (psi1 - psi0) / kappa
//   y1 = y0 + ty * (psi1 - psi0) / kappa          (if |kappa|>0)
// For |kappa|~0 (straight line): dx/dz = tan(psi)=tx/cos, dy/dz = ty/cos.

class NA6PTrackPar
{
 public:
  enum ParLabels : int { kX,
                         kY,
                         kTx,
                         kTy,
                         kQ2Pxz };

  static constexpr float kB2C = -0.299792458e-3f; // GeV/(kG*cm)
  static constexpr float kTinyF = 1e-15f;
  static constexpr float kAlmost1F = 1.f - FLT_EPSILON;
  static constexpr float kSmallBend = 1e-6f;               // threshold on |kappa*dz|
  static constexpr float kSmallKappa = 1e-14f;             // threshold on |kappa|
  static constexpr float ELoss2EKinThreshInv = 1. / 0.025; // do not allow E.Loss correction step with dE/Ekin above the inverse of this value
  static constexpr float InvalidZ = -99999.f;              // tracks with this Z were invalidated
  static constexpr int MaxELossIter = 50;                  // max number of iteration for the ELoss to account for BB dependence on beta*gamma
  static constexpr float kMaxBetaGammaELoss = 1000.;       // limit eloss calculation to this bg.

  NA6PTrackPar() = default;
  NA6PTrackPar(const NA6PTrackPar&) = default;
  NA6PTrackPar(float z, const std::array<float, 5>& par) : mZ{z}, mP{par} {}
  NA6PTrackPar(const std::array<float, 3>& xyz, const std::array<float, 3>& pxyz, int q) { initParam(xyz, pxyz, q); }
  NA6PTrackPar(const float* xyz, const float* pxyz, int q) { initParam(xyz, pxyz, q); }

  ~NA6PTrackPar() = default;
  NA6PTrackPar& operator=(const NA6PTrackPar&) = default;

  void initParam(const std::array<float, 3> xyz, const std::array<float, 3> pxyz, int sign) { initParam(xyz.data(), pxyz.data(), sign); }
  void initParam(const float* xyz, const float* pxyz, int sign);

  void invalidate() { mZ = InvalidZ; }
  bool isValid() const { return mZ != InvalidZ; }

  void setPID(PID pid) { mPID = pid; }
  PID getPID() const { return mPID; }

  void setUserField(uint8_t v) { mUserField = v; }
  uint8_t getUserField() const { return mUserField; }

  float getZ() const { return mZ; }
  float getX() const { return mP[kX]; }
  float getY() const { return mP[kY]; }
  float getTx() const { return mP[kTx]; }       // px/pxz = sin(psi)
  float getTy() const { return mP[kTy]; }       // py/pxz
  float getQ2Pxz() const { return mP[kQ2Pxz]; } // q/pxz
  float getR2() const { return getX() * getX() + getY() * getY(); }

  template <typename T = float>
  std::array<T, 3> getCircleParams(float by) const;

  const auto& getParam() const { return mP; }
  float getParam(int l) const { return mP[l]; }
  void setParam(int l, float v) { mP[l] = v; }
  void incParam(int l, float v) { mP[l] += v; }

  void setZ(float v) { mZ = v; }
  void setX(float v) { mP[kX] = v; }
  void setY(float v) { mP[kY] = v; }
  void setTx(float v) { mP[kTx] = v; }       // px/pxz = sin(psi)
  void setTy(float v) { mP[kTy] = v; }       // py/pxz
  void setQ2Pxz(float v) { mP[kQ2Pxz] = v; } // q/pxz

  template <typename T = float>
  void setXYZ(const T* pos)
  {
    setX(pos[0]);
    setY(pos[1]);
    setZ(pos[2]);
  }
  template <typename T = float>
  void setXYZ(const std::array<T, 3>& pos)
  {
    setXYZ(pos.data());
  }

  float getPxz() const { return mP[kQ2Pxz] != 0.f ? 1.f / std::abs(mP[kQ2Pxz]) : 0.f; }
  int getSign() const { return mP[kQ2Pxz] < 0.f ? -1 : 1; }
  int getCharge() const { return getSign(); }

  // Derived 3D momentum magnitudes
  float getP2Pxz2() const { return 1.f + getTy() * getTy(); }
  float getP2Pxz() const { return std::sqrt(getP2Pxz2()); }
  float getP() const { return getPxz() * getP2Pxz(); }
  float getPz() const { return getPxz() * getCosPsi(); } // pz = pxz*cos(psi)
  float getPx() const { return getPxz() * getTx(); }     // px = pxz*sin(psi)
  float getPy() const { return getPxz() * getTy(); }     // py = pxz*ty
  float getPt() const { return std::hypot(getPx(), getPy()); } // pt = sqrt(px*px + py*py)
  float getQ2P() const;
  float getCharge2Pxz() const { return getQ2Pxz(); } // RSTODO in case of q=0 return 0, while  getQ2Pxz() returns 1/p_xz

  float getSinPsi2() const { return getTx() * getTx(); }
  float getSinPsi() const { return getTx(); }
  float getCosPsi2() const { return getCos2FromSin(getTx()); }
  float getCosPsi() const { return getCosFromSin(getTx()); }
  float getPsi() const { return std::asin(std::max(-1.f, std::min(1.f, getTx()))); }
  float getPhi() const { return std::atan2(getTy(), getTx()); }
  float getTheta() const { return std::acos(getCosPsi() / getP2Pxz()); }
  float getEta() const { return -std::log(std::tan(getTheta() * 0.5f)); }
  float getEnergy() const { return std::hypot(getP(), mPID.getMass()); }
  float getRapidity() const;
  float getCurvature(float by) const { return kB2C * by * getQ2Pxz(); }

  float getR() const { return std::sqrt(getR2()); }
  float getSlopeX() const { return getTx() / getCosPsi(); } // px / pz
  float getSlopeY() const { return getTy() / getCosPsi(); } // py / pz
  std::tuple<float, float> getSlopesXY() const
  {
    auto cosI = 1.f / getCosPsi();
    return {getTx() * cosI, getTy() * cosI};
  }

  template <typename T = float>
  std::array<T, 3> getXYZ() const
  {
    return std::array<T, 3>{getX(), getY(), getZ()};
  }

  template <typename T = float>
  std::array<T, 3> getPXYZ() const
  {
    auto pxz = getPxz();
    return {getTx() * pxz, getTy() * pxz, getCosPsi() * pxz};
  }

  bool getPosDirGlo(std::array<float, 7>& posdirp) const;

  bool propagateParamToZ(float z, float by); // Propagation in dipole field B=(0,By,0)
  bool propagateParamToZ(float z, const float* bxyz);
  bool propagateParamToZ(float z, const std::array<float, 3>& bxyz) { return propagateParamToZ(z, bxyz.data()); }
  bool correctForELoss(float xrho, float density, float atomicZ, float zOverA, bool anglecorr = false);
  bool getZPCAToLine(float XL, float YL, float by, float& zca) const;
  std::string asString() const;

  float getdEdxBB(float betagamma) { return BetheBlochSolid(betagamma); }
  float getdEdxBBOpt(float betagamma, float meanZA = 0.49848) { return BetheBlochSolid(betagamma, meanZA); }
  float getBetheBlochSolidDerivativeApprox(float dedx, float bg, float meanZA) { return BetheBlochSolidDerivative(dedx, bg, meanZA); }

 protected:
  static float getCos2FromSin(float s);
  static float getCosFromSin(float s);

  void clampBetaGamma(float& bg) const { bg = std::min(bg, kMaxBetaGammaELoss); }
  float BetheBlochSolid(float bg, float rho = 2.33, float kp1 = 0.20, float kp2 = 3.00, float meanI = 173e-9, float meanZA = 0.49848);
  float BetheBlochSolidDerivative(float dedx, float bg, float meanZA);
  void g3helx3(float qfield, float step, std::array<float, 7>& vect);

  float mZ = 0.f;
  std::array<float, 5> mP{{0.f, 0.f, 0.f, 0.f, 0.f}};
  PID mPID{PID::Pion};
  uint8_t mUserField = 0;

  ClassDefNV(NA6PTrackPar, 1);
};

inline float NA6PTrackPar::getCos2FromSin(float s)
{
  s = std::clamp(s, -1.f, 1.f);
  return (1.f + s) * (1.f - s);
}

inline float NA6PTrackPar::getCosFromSin(float s)
{
  return std::sqrt(getCos2FromSin(s));
}

inline float NA6PTrackPar::getQ2P() const
{
  const float denom = getP2Pxz();
  return denom > 0.f ? getQ2Pxz() / denom : 0.f; // q/p = (q/pxz)/sqrt(1+ty^2)
}

template <typename T>
std::array<T, 3> NA6PTrackPar::getCircleParams(float by) const
{
  // get circle params as {R, Xcenter, Zcenter}, for straight line just set to current coordinates
  auto crv = getCurvature(by);
  if (std::abs(crv) > kSmallKappa) {
    T crvi = T(1) / crv, sn = getTx(), cs = std::sqrt((T(1) - sn) * (T(1) + sn));
    return {std::abs(crvi), getX() + cs * crvi, getZ() - sn * crvi};
  } else {
    return {T(0), T(getX()), T(getZ())};
  }
}

#endif
