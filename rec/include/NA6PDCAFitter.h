#ifndef _NA6P_DCA_FITTERN_H
#define _NA6P_DCA_FITTERN_H

#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include "NA6PTrack.h"
#include "NA6PHelixHelper.h"

struct FwdTrackCovI {
  float sxx, syy, sxy, szz;

  FwdTrackCovI(const NA6PTrack& trc, float zerrFactor = 1.f) { set(trc, zerrFactor); }
  FwdTrackCovI() = default;
  bool set(const NA6PTrack& trc, float zerrFactor = 1.f)
  {
    float cxx = trc.getSigmaX2(), cyy = trc.getSigmaY2(), cxy = trc.getSigmaYX(), czz = cyy * zerrFactor;
    float detXY = cxx * cyy - cxy * cxy;
    bool res = true;
    if (detXY <= 0.) {
      cxy = std::sqrt(cxx * cyy) * (cxy > 0 ? 0.98f : -0.98f);
      detXY = cxx * cyy - cxy * cxy;
      res = false;
    } 
    auto detXYI = 1. / detXY;
    sxx = cyy * detXYI;
    syy = cxx * detXYI;
    sxy = -cxy * detXYI;
    szz = 1. / czz;
    return res;
  }
};

struct FwdTrackDeriv {
  float dxdz, dydz, d2xdz2, d2ydz2;
  FwdTrackDeriv() = default;
  FwdTrackDeriv(const NA6PTrack& trc, float by) { set(trc, by); }
  void set(const NA6PTrack& trc, float by)
  {
    auto pxyz = trc.getPXYZ();
    if (std::abs(pxyz[2]) < 1e-9) return;
    dxdz = pxyz[0] / pxyz[2];
    dydz = pxyz[1] / pxyz[2];
    // Second derivatives:
    // d^2x/dz^2 = d/dz(px/pz) = 1/pz * (dpx/dz) - px/pz^2 (dpz/dz)
    // Tangents to circumference of radius R in two points separated by dz
    // theta = angle spanned by R
    // dpx/dz = p_xz/R
    // dpz = pxz*(1-cos(theta)) = pxz*2*sin^2(theta/2) ~ pxz*2*(theta/2)^2 ~ (p_xz/R)*dz^2 / 2
    // dpz/dz ~ (p_xz/R)*dz/2 --> negligible
    // d2xdz2 ~ 1/pz * (dpx/dz) = p_xz/pz * 1/R = sqrt((px/pz)^2 + 1) * 1/R
    float crv2c = trc.getCurvature(by);
    d2xdz2 = crv2c * std::sqrt(1.f + dxdz * dxdz);
    //  d^2y/dz^2 = d/dz(py/pz) = 0 (dipole field, py is conserved, and dpz/dz is negligible)
    d2ydz2 = 0.f;
  }
};

template <int N, typename... Args>
class NA6PDCAFitterN
{
  static constexpr double NMin = 2;
  static constexpr double NMax = 4;
  static constexpr double NInv = 1. / N;
  static constexpr int MAXHYP = 2;
  static constexpr float ZerrFactor = 5.; // factor for conversion of track covYY to dummy covXX
  using Vec3D = ROOT::Math::SVector<double, 3>;
  using VecND = ROOT::Math::SVector<double, N>;
  using MatSymND = ROOT::Math::SMatrix<double, N, N, ROOT::Math::MatRepSym<double, N>>;
  using ArrTrack = std::array<NA6PTrack, N>;         // container for prongs (tracks) at single vertex cand.
  using ArrTrackCovI = std::array<FwdTrackCovI, N>;  // container for inv.cov.matrices at single vertex cand.
  using ArrTrDer = std::array<FwdTrackDeriv, N>;     // container of Track 1st and 2nd derivative over their Z param
  using ArrTrPos = std::array<Vec3D, N>;             // container of Track positions
  
 public:

  enum FitStatus : uint8_t { // fit status of crossing hypothesis
    None,                    // no status set (should not be possible!)

    /* Good Conditions */
    Converged, // fit converged
    MaxIter,   // max iterations reached before fit convergence

    /* Error Conditions */
    NoCrossing,      // no reasaonable crossing was found
    RejFiducialVol,  // postion of crossing was not acceptable
    RejTrackZ,       // one candidate track x was below the mimimum required radius
    RejTrackRoughZ,  // rejected by rough cut on tracks Z difference
    RejChi2Max,      // rejected by maximum chi2 cut
    FailProp,        // propagation of at least prong to PCA failed
    FailInvCov,      // inversion of cov.-matrix failed
    FailInvWeight,   // inversion of Ti weight matrix failed
    FailInv2ndDeriv, // inversion of 2nd derivatives failed
    FailCorrTracks,  // correction of tracks to updated x failed
    FailCloserAlt,   // alternative PCA is closer
    //
    NStatusesDefined
  };

  enum BadCovPolicy : uint8_t { // if encountering non-positive defined cov. matrix, the choice is:
    Discard = 0,                // stop evaluation
    Override = 1,               // override correlation coef. to have cov.matrix pos.def and continue
    OverrideAndFlag = 2         // override correlation coef. to have cov.matrix pos.def, set mPropFailed flag of corresponding candidate to true and continue (up to the user to check the flag)
  };

  NA6PDCAFitterN() = default;
  NA6PDCAFitterN(float by, bool useAbsDCA, bool prop2DCA) : mBy(by), mUseAbsDCA(useAbsDCA), mPropagateToPCA(prop2DCA)
  {
    static_assert(N >= NMin && N <= NMax, "N prongs outside of allowed range");
  }
  static constexpr int getNProngs() { return N; }
  template <class... Tr>
  int process(const Tr&... args);
  bool propagateTracksToVertex(int cand = 0) {return false;} // TODO
  bool recalculatePCAWithErrors(int cand = 0) {return false;} // TODO
  void setBadCovPolicy(BadCovPolicy v) { mBadCovPolicy = v; }
  BadCovPolicy getBadCovPolicy() const { return mBadCovPolicy; }

 protected:
  bool calcPCACoefs() { return true;} // TODO;
  void calcPCA() {}; //TODO
  void calcResidDerivatives() {}; //TODO
  void calcTrackResiduals();
  void calcTrackDerivatives();
  void calcChi2Derivatives() {}; //TODO
  double calcChi2() const;
  bool minimizeChi2();
  bool minimizeChi2NoErr() { return true; } // TODO
  bool roughDZCut() const;
  bool closerToAlternative() const {return false;} //TODO
  bool propagateToZ(NA6PTrack& t, float z) { return true;} // TODO

 private:
  VecND mDChi2Dz;                             // 1st derivatives of chi2 over tracks X params
  MatSymND mD2Chi2Dz2;                        // 2nd derivatives of chi2 over tracks Z params (symmetric matrix)
  std::array<const NA6PTrack*, N> mOrigTrPtr;
  std::array<CircleXZ, N> mTrAux;             // Aux track info for each track at each cand. vertex
  std::array<int, MAXHYP> mOrder{0};
  std::array<ArrTrack, MAXHYP> mCandTr;       // tracks at each cond. vertex (Note: Errors are at seed XZ point)
  std::array<ArrTrackCovI, MAXHYP> mTrcEInv;  // errors for each track at each cand. vertex
  std::array<ArrTrDer, MAXHYP> mTrDer;        // Track derivatives
  std::array<ArrTrPos, MAXHYP> mTrPos;        // Track positions
  std::array<ArrTrPos, MAXHYP> mTrRes;        // Track residuals
  std::array<Vec3D, MAXHYP> mPCA;             // PCA for each vertex candidate
  std::array<float, MAXHYP> mChi2 = {0};      // Chi2 at PCA candidate
  std::array<int, MAXHYP> mNIters;            // number of iterations for each seed
  std::array<bool, MAXHYP> mTrPropDone{};     // Flag that the tracks are fully propagated to PCA
  std::array<bool, MAXHYP> mPropFailed{};     // Flag that some propagation failed for this PCA candidate
  std::array<FitStatus, MAXHYP> mFitStatus{}; // fit status of each hypothesis fit
  CrossInfo mCrossings;                       // info on track crossing
  int mCurHyp = 0;
  int mCrossIDCur = 0;
  int mCrossIDAlt = -1;
  BadCovPolicy mBadCovPolicy{BadCovPolicy::Discard}; // what to do in case of non-pos-def. cov. matrix, see BadCovPolicy enum
  bool mAllowAltPreference = true;            // if the fit converges to alternative PCA seed, abandon the current one
  bool mUseAbsDCA = false;                    // use abs. distance minimization rather than chi2
  bool mWeightedFinalPCA = false;             // recalculate PCA as a cov-matrix weighted mean, even if absDCA method was used
  bool mPropagateToPCA = true;                // create tracks version propagated to PCA
  bool mIsCollinear = false;                  // use collinear fits when there 2 crossing points
  int mMaxIter = 20;                          // max number of iterations
  float mBy = 0;                              // by field, to be set by user
  float mMaxVertX = 30.;                      // maximum vertex X coordinate
  float mMinVertZ = -10.;                     // minimum vertex Z coordinate
  float mMaxVertZ = 40.;                      // maximum vertex Z coordinate
  float mMaxDZIni = 4.;                       // reject (if>0) PCA candidate if tracks DZ exceeds threshold
  float mMaxDXZIni = 4.;                      // reject (if>0) PCA candidate if tracks dXY exceeds threshold
  float mMinParamChange = 1e-3;               // stop iterations if largest change of any X is smaller than this
  float mMinRelChi2Change = 0.9;              // stop iterations is chi2/chi2old > this
  float mMaxChi2 = 100;                       // abs cut on chi2 or abs distance
  float mMaxDist2ToMergeSeeds = 1.;           // merge 2 seeds to their average if their distance^2 is below the threshold

  template <class T, class... Tr>
  void assign(int i, const T& t, const Tr&... args)
    {
      mOrigTrPtr[i] = &t;
      assign(i + 1, args...);
    }

  static void setTrackPos(Vec3D& pnt, const NA6PTrack& tr)
  {
    pnt[0] = tr.getX();
    pnt[1] = tr.getY();
    pnt[2] = tr.getZ();
  }

 void clear()
  {
    mCurHyp = 0;
    mAllowAltPreference = true;
    mOrder.fill(0);
    mPropFailed.fill(false);
    mTrPropDone.fill(false);
    mNIters.fill(0);
    mChi2.fill(-1);
    mFitStatus.fill(FitStatus::None);
  }
  
  ClassDefNV(NA6PDCAFitterN, 1);
};

///_________________________________________________________________________
template <int N, typename... Args>
template <class... Tr>
int NA6PDCAFitterN<N, Args...>::process(const Tr&... args)
{
  // This is a main entry point: fit PCA of N tracks
  static_assert(sizeof...(args) == N, "incorrect number of input tracks");
  assign(0, args...);
  clear();
  for (int i = 0; i < N; i++) {
    mTrAux[i].set(*mOrigTrPtr[i], mBy);
  }
  if (!mCrossings.set(mTrAux[0], *mOrigTrPtr[0], mTrAux[1], *mOrigTrPtr[1], mMaxDXZIni, mIsCollinear)) { // even for N>2 it should be enough to test just 1 loop
    mFitStatus[mCurHyp] = FitStatus::NoCrossing;
    return 0;
  }
  if (mCrossings.nDCA == MAXHYP) { // if there are 2 candidates and they are too close, chose their mean as a starting point
    auto dst2 = (mCrossings.xDCA[0] - mCrossings.xDCA[1]) * (mCrossings.xDCA[0] - mCrossings.xDCA[1]) +
                (mCrossings.zDCA[0] - mCrossings.zDCA[1]) * (mCrossings.zDCA[0] - mCrossings.zDCA[1]);
    if (dst2 < mMaxDist2ToMergeSeeds) {
      mCrossings.nDCA = 1;
      mCrossings.xDCA[0] = 0.5 * (mCrossings.xDCA[0] + mCrossings.xDCA[1]);
      mCrossings.zDCA[0] = 0.5 * (mCrossings.zDCA[0] + mCrossings.zDCA[1]);
    }
  }
  // check all crossings
  for (int ic = 0; ic < mCrossings.nDCA; ic++) {
    // check if radius is acceptable
    if (mCrossings.zDCA[ic] > mMaxVertZ ||
        mCrossings.zDCA[ic] < mMinVertZ ||
        std::abs(mCrossings.xDCA[ic]) > mMaxVertX) {
      mFitStatus[mCurHyp] = FitStatus::RejFiducialVol;
      continue;
    }
    mCrossIDCur = ic;
    mCrossIDAlt = (mCrossings.nDCA == 2 && mAllowAltPreference) ? 1 - ic : -1; // works for max 2 crossings
    mPCA[mCurHyp][0] = mCrossings.xDCA[ic];
    mPCA[mCurHyp][2] = mCrossings.zDCA[ic];

    if (mUseAbsDCA ? minimizeChi2NoErr() : minimizeChi2()) {
      mOrder[mCurHyp] = mCurHyp;
      if (mPropagateToPCA && !propagateTracksToVertex(mCurHyp)) {
        continue; // discard candidate if failed to propagate to it
      }
      mCurHyp++;
    }
  }

  for (int i = mCurHyp; i--;) { // order in quality
    for (int j = i; j--;) {
      if (mChi2[mOrder[i]] < mChi2[mOrder[j]]) {
        std::swap(mOrder[i], mOrder[j]);
      }
    }
  }
  if (mUseAbsDCA && mWeightedFinalPCA) {
    for (int i = mCurHyp; i--;) {
      recalculatePCAWithErrors(i);
    }
  }
  return mCurHyp;
}

//___________________________________________________________________
template <int N, typename... Args>
bool NA6PDCAFitterN<N, Args...>::roughDZCut() const
{
  // apply rough cut on DZ between the tracks in the seed point
  bool accept = true;
  for (int i = N; accept && i--;) {
    for (int j = i; j--;) {
      if (std::abs(mCandTr[mCurHyp][i].getZ() - mCandTr[mCurHyp][j].getZ()) > mMaxDZIni) {
        accept = false;
        break;
      }
    }
  }
  return accept;
}

//___________________________________________________________________
template <int N, typename... Args>
void NA6PDCAFitterN<N, Args...>::calcTrackResiduals()
{
  // calculate residuals in the global frame
  for (int i = N; i--;) {
    mTrRes[mCurHyp][i] = mTrPos[mCurHyp][i] - mPCA[mCurHyp];
  }
}
//___________________________________________________________________
template <int N, typename... Args>
void NA6PDCAFitterN<N, Args...>::calcTrackDerivatives()
{
  // calculate track derivatives over Z param
  for (int i = N; i--;) {
    mTrDer[mCurHyp][i].set(mCandTr[mCurHyp][i], mBy);
  }
}
//___________________________________________________________________
template <int N, typename... Args>
double NA6PDCAFitterN<N, Args...>::calcChi2() const
{
  // calculate current chi2
  double chi2 = 0;
  for (int i = N; i--;) {
    const auto& res = mTrRes[mCurHyp][i];
    const auto& covI = mTrcEInv[mCurHyp][i];
    chi2 += res[0] * res[0] * covI.sxx + res[1] * res[1] * covI.syy + res[2] * res[2] * covI.szz + 2. * res[0] * res[1] * covI.sxy;
  }
  return chi2;
}
//___________________________________________________________________
template <int N, typename... Args>
bool NA6PDCAFitterN<N, Args...>::minimizeChi2()
{
  // find best chi2 (weighted DCA) of N tracks in the vicinity of the seed PCA
  for (int i = N; i--;) {
    mCandTr[mCurHyp][i] = *mOrigTrPtr[i];
    auto z = mPCA[mCurHyp][2];
    if (z > mMaxVertZ || z < mMinVertZ) {
      mFitStatus[mCurHyp] = FitStatus::RejTrackZ;
      return false;
    }
    if (!propagateToZ(mCandTr[mCurHyp][i], z)) {
      return false;
    }
    setTrackPos(mTrPos[mCurHyp][i], mCandTr[mCurHyp][i]);             // prepare positions
    if (!mTrcEInv[mCurHyp][i].set(mCandTr[mCurHyp][i], ZerrFactor)) { // prepare inverse cov.matrices at starting point
      mFitStatus[mCurHyp] = FitStatus::FailInvCov;
      if (mBadCovPolicy == Discard) {
        return false;
      } else if (mBadCovPolicy == OverrideAndFlag) {
        mPropFailed[mCurHyp] = true;
      } // otherwise, just use overridden errors w/o flagging
    }
  }

  if (mMaxDZIni > 0 && !roughDZCut()) { // apply rough cut on tracks Z difference
    mFitStatus[mCurHyp] = FitStatus::RejTrackRoughZ;
    return false;
  }

  if (!calcPCACoefs()) { // prepare tracks contribution matrices to the global PCA
    return false;
  }
  calcPCA();            // current PCA
  calcTrackResiduals(); // current track residuals
  float chi2Upd, chi2 = calcChi2();
  do {
    calcTrackDerivatives(); // current track derivatives (1st and 2nd)
    calcResidDerivatives(); // current residals derivatives (1st and 2nd)
    calcChi2Derivatives();  // current chi2 derivatives (1st and 2nd)

    // do Newton-Rapson iteration with corrections = - dchi2/d{x0..xN} * [ d^2chi2/d{x0..xN}^2 ]^-1
    if (!mD2Chi2Dz2.Invert()) {
      mFitStatus[mCurHyp] = FitStatus::FailInv2ndDeriv;
      return false;
    }
    VecND dz = mD2Chi2Dz2 * mDChi2Dz;
    if (!correctTracks(dz)) {
      mFitStatus[mCurHyp] = FitStatus::FailCorrTracks;
      return false;
    }
    calcPCA(); // updated PCA
    if (mCrossIDAlt >= 0 && closerToAlternative()) {
      mFitStatus[mCurHyp] = FitStatus::FailCloserAlt;
      mAllowAltPreference = false;
      return false;
    }
    calcTrackResiduals(); // updated residuals
    chi2Upd = calcChi2(); // updated chi2
    if (getAbsMax(dz) < mMinParamChange || chi2Upd > chi2 * mMinRelChi2Change) {
      chi2 = chi2Upd;
      mFitStatus[mCurHyp] = FitStatus::Converged;
      break; // converged
    }
    chi2 = chi2Upd;
  } while (++mNIters[mCurHyp] < mMaxIter);
  if (mNIters[mCurHyp] == mMaxIter) {
    mFitStatus[mCurHyp] = FitStatus::MaxIter;
  }
  //
  mChi2[mCurHyp] = chi2 * NInv;
  if (mChi2[mCurHyp] >= mMaxChi2) {
    mFitStatus[mCurHyp] = FitStatus::RejChi2Max;
    return false;
  }
  return true;
}

using NA6PDCAFitter2 = NA6PDCAFitterN<2, NA6PTrack>;
using NA6PDCAFitter3 = NA6PDCAFitterN<3, NA6PTrack>;

#endif
