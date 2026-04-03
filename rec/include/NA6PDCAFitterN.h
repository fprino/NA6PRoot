#ifndef _NA6P_DCA_FITTERN_H
#define _NA6P_DCA_FITTERN_H

#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <fairlogger/Logger.h>
#include "NA6PTrack.h"
#include "NA6PHelixHelper.h"

struct TrackCovI {
  float sxx, syy, sxy, szz;

  TrackCovI(const NA6PTrack& trc, float zerrFactor = 1.f) { set(trc, zerrFactor); }
  TrackCovI() = default;
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

struct TrackDeriv {
  float dxdz, dydz, d2xdz2, d2ydz2;
  TrackDeriv() = default;
  TrackDeriv(const NA6PTrack& trc, float by) { set(trc, by); }
  void set(const NA6PTrack& trc, float by)
  {
    auto pxyz = trc.getPXYZ();
    if (std::abs(pxyz[2]) < 1e-9)
      return;
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
  using MatSym3D = ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>>;
  using MatStd3D = ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepStd<double, 3>>;
  using MatSymND = ROOT::Math::SMatrix<double, N, N, ROOT::Math::MatRepSym<double, N>>;
  using ArrTrack = std::array<NA6PTrack, N>;     // container for prongs (tracks) at single vertex cand.
  using ArrTrackCovI = std::array<TrackCovI, N>; // container for inv.cov.matrices at single vertex cand.
  using ArrTrCoef = std::array<MatStd3D, N>;     // container of TrackCoefVtx coefficients at single vertex cand.
  using ArrTrDer = std::array<TrackDeriv, N>;    // container of Track 1st and 2nd derivative over their Z param
  using ArrTrPos = std::array<Vec3D, N>;         // container of Track positions

 public:
  enum BadCovPolicy : uint8_t { // if encountering non-positive defined cov. matrix, the choice is:
    Discard = 0,                // stop evaluation
    Override = 1,               // override correlation coef. to have cov.matrix pos.def and continue
    OverrideAndFlag = 2         // override correlation coef. to have cov.matrix pos.def, set mPropFailed flag of corresponding candidate to true and continue (up to the user to check the flag)
  };

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

  static constexpr int getNProngs() { return N; }

  NA6PDCAFitterN() = default;
  NA6PDCAFitterN(float by, bool useAbsDCA, bool prop2DCA) : mBy(by), mUseAbsDCA(useAbsDCA), mPropagateToPCA(prop2DCA)
  {
    static_assert(N >= NMin && N <= NMax, "N prongs outside of allowed range");
  }

  //=========================================================================
  ///< return PCA candidate, by default best on is provided (no check for the index validity)
  const Vec3D& getPCACandidate(int cand = 0) const { return mPCA[mOrder[cand]]; }
  const auto getPCACandidatePos(int cand = 0) const
  {
    const auto& vd = mPCA[mOrder[cand]];
    return std::array<float, 3>{static_cast<float>(vd[0]), static_cast<float>(vd[1]), static_cast<float>(vd[2])};
  }

  ///< return position of quality-ordered candidate in the internal structures
  int getCandidatePosition(int cand = 0) const { return mOrder[cand]; }

  ///< return Chi2 at PCA candidate (no check for its validity)
  float getChi2AtPCACandidate(int cand = 0) const { return mChi2[mOrder[cand]]; }

  ///< prepare copies of tracks at the V0 candidate (no check for the candidate validity)
  ///  must be called before getTrack(i,cand) query
  bool propagateTracksToVertex(int cand = 0);

  ///< check if propagation of tracks to candidate vertex was done
  bool isPropagateTracksToVertexDone(int cand = 0) const { return mTrPropDone[mOrder[cand]]; }

  ///< check if propagation of tracks to candidate vertex failed
  bool isPropagationFailure(int cand = 0) const { return mPropFailed[mOrder[cand]]; }

  ///< track param propagated to V0 candidate (no check for the candidate validity)
  ///  propagateTracksToVertex must be called in advance
  NA6PTrack& getTrack(int i, int cand = 0)
  {
    if (!mTrPropDone[mOrder[cand]]) {
      LOGP(error, "propagateTracksToVertex was not called yet");
    }
    return mCandTr[mOrder[cand]][i];
  }

  const NA6PTrack& getTrack(int i, int cand = 0) const
  {
    if (!mTrPropDone[mOrder[cand]]) {
      LOGP(error, "propagateTracksToVertex was not called yet");
    }
    return mCandTr[mOrder[cand]][i];
  }

  ///< create parent track param with errors for decay vertex
  NA6PTrack createParentTrackParCov(int cand = 0) const; // TODO: update to use NA6PTRackParCov

  ///< create parent track param w/o errors for decay vertex
  NA6PTrack createParentTrackPar(int cand = 0) const;

  ///< recalculate PCA as a cov-matrix weighted mean, even if absDCA method was used
  bool recalculatePCAWithErrors(int cand = 0);

  MatSym3D calcPCACovMatrix(int cand = 0) const;

  std::array<float, 6> calcPCACovMatrixFlat(int cand = 0) const
  {
    auto m = calcPCACovMatrix(cand);
    return {static_cast<float>(m(0, 0)), static_cast<float>(m(1, 0)), static_cast<float>(m(1, 1)), static_cast<float>(m(2, 0)), static_cast<float>(m(2, 1)), static_cast<float>(m(2, 2))};
  }

  const NA6PTrack* getOrigTrackPtr(int i) const { return mOrigTrPtr[i]; }

  FitStatus getFitStatus(int cand = 0) const noexcept { return mFitStatus[mOrder[cand]]; }

  ///< return number of iterations during minimization (no check for its validity)
  int getNIterations(int cand = 0) const { return mNIters[mOrder[cand]]; }
  void setPropagateToPCA(bool v = true) { mPropagateToPCA = v; }
  void setMaxIter(int n = 20) { mMaxIter = n > 2 ? n : 2; }
  void setMaxDZIni(float d = 4.) { mMaxDZIni = d; }
  void setMaxDXZIni(float d = 4.) { mMaxDXZIni = d > 0 ? d : 1e9; }
  void setMaxChi2(float chi2 = 999.) { mMaxChi2 = chi2; }
  void setBy(float by) { mBy = std::abs(by) > 1.e-12 ? by : 0.f; }
  void setMinParamChange(float x = 1e-3) { mMinParamChange = x > 1e-4 ? x : 1.e-4; }
  void setMinRelChi2Change(float r = 0.9) { mMinRelChi2Change = r > 0.1 ? r : 999.; }
  void setUseAbsDCA(bool v) { mUseAbsDCA = v; }
  void setWeightedFinalPCA(bool v) { mWeightedFinalPCA = v; }
  void setMaxDistance2ToMerge(float v) { mMaxDist2ToMergeSeeds = v; }
  void setUsePropagator(bool v) { mUsePropagator = v; }
  void setRefitWithMatCorr(bool v) { mRefitWithMatCorr = v; }
  void setMaxVertX(float x) { mMaxVertX = x; }
  void setMinVertZ(float z) { mMinVertZ = z; }
  void setMaxVertZ(float z) { mMaxVertZ = z; }
  void setCollinear(bool isCollinear) { mIsCollinear = isCollinear; }

  int getNCandidates() const { return mCurHyp; }
  int getMaxIter() const { return mMaxIter; }
  float getMaxDZIni() const { return mMaxDZIni; }
  float getMaxDXZIni() const { return mMaxDXZIni; }
  float getMaxChi2() const { return mMaxChi2; }
  float getMinParamChange() const { return mMinParamChange; }
  float getBy() const { return mBy; }
  float getMaxDistance2ToMerge() const { return mMaxDist2ToMergeSeeds; }
  bool getUseAbsDCA() const { return mUseAbsDCA; }
  bool getWeightedFinalPCA() const { return mWeightedFinalPCA; }
  bool getPropagateToPCA() const { return mPropagateToPCA; }
  bool getUsePropagator() const { return mUsePropagator; }
  bool getRefitWithMatCorr() const { return mRefitWithMatCorr; }
  float getMaxVertX() const { return mMaxVertX; }
  float getMinVertZ() const { return mMinVertZ; }
  float getMaxVertZ() const { return mMaxVertZ; }

  template <class... Tr>
  int process(const Tr&... args);

 protected:
  bool calcPCACoefs();
  bool calcInverseWeight();
  void calcResidDerivatives();
  void calcResidDerivativesNoErr();
  void calcChi2Derivatives();
  void calcChi2DerivativesNoErr();
  void calcPCA();
  void calcPCANoErr();
  void calcTrackResiduals();
  void calcTrackDerivatives();
  double calcChi2() const;
  double calcChi2NoErr() const;
  bool correctTracks(const VecND& corrZ);
  bool minimizeChi2();
  bool minimizeChi2NoErr();
  bool roughDZCut() const;
  bool closerToAlternative() const;
  bool propagateToZ(NA6PTrack& t, float z) const;
  bool propagateParamToZ(NA6PTrack& t, float z) const;

  static double getAbsMax(const VecND& v);
  ///< track param positions at V0 candidate (no check for the candidate validity)
  const Vec3D& getTrackPos(int i, int cand = 0) const { return mTrPos[mOrder[cand]][i]; }

  void assign(int) {}

  template <class T, class... Tr>
  void assign(int i, const T& t, const Tr&... args)
  {
    mOrigTrPtr[i] = &t;
    assign(i + 1, args...);
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

  static void setTrackPos(Vec3D& pnt, const NA6PTrack& tr)
  {
    pnt[0] = tr.getX();
    pnt[1] = tr.getY();
    pnt[2] = tr.getZ();
  }

  void setBadCovPolicy(BadCovPolicy v) { mBadCovPolicy = v; }
  BadCovPolicy getBadCovPolicy() const { return mBadCovPolicy; }

 private:
  // vectors of 1st derivatives of track local residuals over Z parameters
  std::array<std::array<Vec3D, N>, N> mDResidDz;
  // vectors of 1nd derivatives of track local residuals over Z parameters
  // (cross-derivatives DR/(dx_j*dx_k) = 0 for j!=k, therefore the hessian is diagonal)
  std::array<std::array<Vec3D, N>, N> mD2ResidDz2;
  VecND mDChi2Dz;                             // 1st derivatives of chi2 over tracks X params
  MatSymND mD2Chi2Dz2;                        // 2nd derivatives of chi2 over tracks Z params (symmetric matrix)
  std::array<const NA6PTrack*, N> mOrigTrPtr; // Pointers to original tracks
  std::array<CircleXZ, N> mTrAux;             // Aux track info for each track at each cand. vertex
  MatSym3D mWeightInv;                        // inverse weight of single track, [sum{M^T E M}]^-1 in EQ.T
  std::array<int, MAXHYP> mOrder{0};
  std::array<ArrTrack, MAXHYP> mCandTr;       // tracks at each cond. vertex (Note: Errors are at seed XZ point)
  std::array<ArrTrCoef, MAXHYP> mTrCFVT;      // TrackCoefVtx for each track at each cand. vertex
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
  bool mAllowAltPreference = true;                   // if the fit converges to alternative PCA seed, abandon the current one
  bool mUseAbsDCA = false;                           // use abs. distance minimization rather than chi2
  bool mWeightedFinalPCA = false;                    // recalculate PCA as a cov-matrix weighted mean, even if absDCA method was used
  bool mPropagateToPCA = true;                       // create tracks version propagated to PCA
  bool mUsePropagator = false;                       // use propagator with 3D B-field, set automatically if material correction is requested
  bool mRefitWithMatCorr = false;                    // when doing propagateTracksToVertex, propagate tracks to V0 with material corrections and rerun minimization again
  bool mIsCollinear = false;                         // use collinear fits when there 2 crossing points
  int mMaxIter = 20;                                 // max number of iterations
  float mBy = 0;                                     // by field, to be set by user
  float mMaxVertX = 30.;                             // maximum vertex X coordinate
  float mMinVertZ = -10.;                            // minimum vertex Z coordinate
  float mMaxVertZ = 40.;                             // maximum vertex Z coordinate
  float mMaxDZIni = 4.;                              // reject (if>0) PCA candidate if tracks DZ exceeds threshold
  float mMaxDXZIni = 4.;                             // reject (if>0) PCA candidate if tracks dXY exceeds threshold
  float mMinParamChange = 1e-3;                      // stop iterations if largest change of any X is smaller than this
  float mMinRelChi2Change = 0.9;                     // stop iterations is chi2/chi2old > this
  float mMaxChi2 = 100;                              // abs cut on chi2 or abs distance
  float mMaxDist2ToMergeSeeds = 1.;                  // merge 2 seeds to their average if their distance^2 is below the threshold

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

//__________________________________________________________________________
template <int N, typename... Args>
bool NA6PDCAFitterN<N, Args...>::calcPCACoefs()
{
  //< calculate Ti matrices for global vertex decomposition to V = sum_{0<i<N} Ti pi, see EQ.T in the ref
  if (!calcInverseWeight()) {
    mFitStatus[mCurHyp] = FitStatus::FailInvWeight;
    return false;
  }
  for (int i = N; i--;) { // build Mi*Ei matrix
    const auto& taux = mTrAux[i];
    const auto& tcov = mTrcEInv[mCurHyp][i];
    MatStd3D miei;
    miei[0][0] = tcov.sxx;
    miei[0][1] = tcov.sxy;
    miei[0][2] = 0;
    miei[1][0] = tcov.sxy;
    miei[1][1] = tcov.syy;
    miei[1][2] = 0;
    miei[2][0] = 0;
    miei[2][1] = 0;
    miei[2][2] = tcov.szz;
    mTrCFVT[mCurHyp][i] = mWeightInv * miei;
  }
  return true;
}
//__________________________________________________________________________
template <int N, typename... Args>
bool NA6PDCAFitterN<N, Args...>::calcInverseWeight()
{
  //< calculate [sum_{0<j<N} M_j*E_j*M_j^T]^-1 used for Ti matrices, see EQ.T
  auto* arrmat = mWeightInv.Array();
  memset(arrmat, 0, sizeof(mWeightInv));
  enum { XX,
         XY,
         YY,
         XZ,
         YZ,
         ZZ };
  for (int i = N; i--;) {
    const auto& taux = mTrAux[i];
    const auto& tcov = mTrcEInv[mCurHyp][i];
    arrmat[XX] += tcov.sxx;
    arrmat[XY] += tcov.sxy;
    arrmat[XZ] += 0;
    arrmat[YY] += tcov.syy;
    arrmat[YZ] += 0;
    arrmat[ZZ] += tcov.szz;
  }
  // invert 3x3 symmetrix matrix
  return mWeightInv.Invert();
}
//__________________________________________________________________________
template <int N, typename... Args>
void NA6PDCAFitterN<N, Args...>::calcResidDerivatives()
{
  //< calculate matrix of derivatives for weighted chi2: residual i vs parameter X of track j
  for (int i = N; i--;) {                     // residual being differentiated
    for (int j = N; j--;) {                   // track over which we differentiate
      const auto& matT = mTrCFVT[mCurHyp][j]; // coefficient matrix for track J
      const auto& trDz = mTrDer[mCurHyp][j];  // track point derivs over track Z param
      auto& dr1 = mDResidDz[i][j];
      auto& dr2 = mD2ResidDz2[i][j];

      // calculate DResid_i/Dz_j = (delta_ij - M_i^tr * T_j) * DTrack_k/Dz_k
      dr1[0] = -(matT[0][0] * trDz.dxdz + matT[0][1] * trDz.dydz + matT[0][2]);
      dr1[1] = -(matT[1][0] * trDz.dxdz + matT[1][1] * trDz.dydz + matT[1][2]);
      dr1[2] = -(matT[2][0] * trDz.dxdz + matT[2][1] * trDz.dydz + matT[2][2]);

      // calculate D2Resid_I/(Dz_J Dz_K) = (delta_ijk - M_i^tr * T_j * delta_jk) * D2Track_k/dx_k^2
      dr2[0] = -(matT[0][0] * trDz.d2xdz2 + matT[0][1] * trDz.d2ydz2);
      dr2[1] = -(matT[1][0] * trDz.d2xdz2 + matT[1][1] * trDz.d2ydz2);
      dr2[2] = -(matT[2][0] * trDz.d2xdz2 + matT[2][1] * trDz.d2ydz2);

      if (i == j) {
        dr1[0] += trDz.dxdz;
        dr1[1] += trDz.dydz;
        dr1[2] += 1.;

        dr2[0] += trDz.d2xdz2;
        dr2[1] += trDz.d2ydz2;
        // dr2[2] += 0 (Since d2z/dz2 is identically zero)
      }
    } // track over which we differentiate
  }   // residual being differentiated
}

//__________________________________________________________________________
template <int N, typename... Args>
void NA6PDCAFitterN<N, Args...>::calcResidDerivativesNoErr()
{
  //< calculate matrix of derivatives for absolute distance chi2: residual i vs parameter X of track j
  constexpr double NInv1 = 1. - NInv;       // profit from Rii = I/Ninv
  for (int i = N; i--;) {                   // residual being differentiated
    const auto& trDzi = mTrDer[mCurHyp][i]; // track point derivs over track Z param
    auto& dr1ii = mDResidDz[i][i];
    auto& dr2ii = mD2ResidDz2[i][i];
    dr1ii[0] = NInv1 * trDzi.dxdz;
    dr1ii[1] = NInv1 * trDzi.dydz;
    dr1ii[2] = NInv1;

    dr2ii[0] = NInv1 * trDzi.d2xdz2;
    dr2ii[1] = NInv1 * trDzi.d2ydz2;
    dr2ii[2] = 0;

    for (int j = i; j--;) { // track over which we differentiate
      auto& dr1ij = mDResidDz[i][j];
      auto& dr1ji = mDResidDz[j][i];
      const auto& trDzj = mTrDer[mCurHyp][j]; // track point derivs over track Z param

      // calculate DResid_i/Dz_j = (delta_ij - R_ij) * DTrack_j/Dz_j  for j<i
      dr1ij[0] = -trDzj.dxdz * NInv;
      dr1ij[1] = -trDzj.dydz * NInv;
      dr1ij[2] = -NInv;

      // calculate DResid_j/Dz_i = (delta_ij - R_ji) * DTrack_i/Dz_i  for j<i
      dr1ji[0] = -trDzi.dxdz * NInv;
      dr1ji[1] = -trDzi.dydz * NInv;
      dr1ji[2] = -NInv;

      auto& dr2ij = mD2ResidDz2[i][j];
      auto& dr2ji = mD2ResidDz2[j][i];
      // calculate D2Resid_I/(Dz_J Dz_K) = (delta_ij - Rij) * D2Track_j/dz_j^2 * delta_jk for j<i
      dr2ij[0] = -trDzj.d2xdz2 * NInv;
      dr2ij[1] = -trDzj.d2ydz2 * NInv;
      dr2ij[2] = 0;

      // calculate D2Resid_j/(Dx_i Dx_k) = (delta_ij - Rji) * D2Track_i/dx_i^2 * delta_ik for j<i
      dr2ji[0] = -trDzi.d2xdz2 * NInv;
      dr2ji[1] = -trDzi.d2ydz2 * NInv;
      dr2ji[2] = 0;

    } // track over which we differentiate
  }   // residual being differentiated
}

//__________________________________________________________________________
template <int N, typename... Args>
void NA6PDCAFitterN<N, Args...>::calcChi2Derivatives()
{
  //< calculate 1st and 2nd derivatives of wighted DCA (chi2) over track parameters X, see EQ.Chi2 in the ref
  std::array<std::array<Vec3D, N>, N> covIDrDz; // tempory vectors of covI_j * dres_j/dz_i

  // chi2 1st derivative
  for (int i = N; i--;) {
    auto& dchi1 = mDChi2Dz[i]; // DChi2/Dz_i = sum_j { res_j * covI_j * Dres_j/Dz_i }
    dchi1 = 0;
    for (int j = N; j--;) {
      const auto& res = mTrRes[mCurHyp][j];    // vector of residuals of track j
      const auto& covI = mTrcEInv[mCurHyp][j]; // inverse cov matrix of track j
      const auto& dr1 = mDResidDz[j][i];       // vector of j-th residuals 1st derivative over Z param of track i
      auto& cidr = covIDrDz[i][j];             // vector covI_j * dres_j/dz_i, save for 2nd derivative calculation
      cidr[0] = covI.sxx * dr1[0] + covI.sxy * dr1[1];
      cidr[1] = covI.sxy * dr1[0] + covI.syy * dr1[1];
      cidr[2] = covI.szz * dr1[2];
      // calculate res_i * covI_j * dres_j/dx_i
      dchi1 += ROOT::Math::Dot(res, cidr);
    }
  }
  // chi2 2nd derivative
  for (int i = N; i--;) {
    for (int j = i + 1; j--;) {       // symmetric matrix
      auto& dchi2 = mD2Chi2Dz2[i][j]; // D2Chi2/Dz_i/Dz_j = sum_k { Dres_k/Dz_j * covI_k * Dres_k/Dz_i + res_k * covI_k * D2res_k/Dz_i/Dz_j }
      dchi2 = 0;
      for (int k = N; k--;) {
        const auto& dr1j = mDResidDz[k][j];  // vector of k-th residuals 1st derivative over Z param of track j
        const auto& cidrkj = covIDrDz[i][k]; // vector covI_k * dres_k/dz_i
        dchi2 += ROOT::Math::Dot(dr1j, cidrkj);
        if (k == j) {
          const auto& res = mTrRes[mCurHyp][k];    // vector of residuals of track k
          const auto& covI = mTrcEInv[mCurHyp][k]; // inverse cov matrix of track k
          const auto& dr2ij = mD2ResidDz2[k][j];   // vector of k-th residuals 2nd derivative over Z params of track j
          dchi2 += res[0] * (covI.sxx * dr2ij[0] + covI.sxy * dr2ij[1]) + res[1] * (covI.sxy * dr2ij[0] + covI.syy * dr2ij[1]) + res[2] * covI.szz * dr2ij[2];
        }
      }
    }
  }
}
//__________________________________________________________________________
template <int N, typename... Args>
void NA6PDCAFitterN<N, Args...>::calcChi2DerivativesNoErr()
{
  //< calculate 1st and 2nd derivatives of abs DCA (chi2) over track parameters X, see (6) in the ref
  for (int i = N; i--;) {
    auto& dchi1 = mDChi2Dz[i]; // DChi2/Dz_i = sum_j { res_j * Dres_j/Dz_i }
    dchi1 = 0;                 // chi2 1st derivative
    for (int j = N; j--;) {
      const auto& res = mTrRes[mCurHyp][j]; // vector of residuals of track j
      const auto& dr1 = mDResidDz[j][i];    // vector of j-th residuals 1st derivative over Z param of track i
      dchi1 += ROOT::Math::Dot(res, dr1);
      if (i >= j) { // symmetrix matrix
        // chi2 2nd derivative
        auto& dchi2 = mD2Chi2Dz2[i][j]; // D2Chi2/Dz_i/Dz_j = sum_k { Dres_k/Dz_j * covI_k * Dres_k/Dz_i + res_k * covI_k * D2res_k/Dz_i/Dz_j }
        dchi2 = ROOT::Math::Dot(mTrRes[mCurHyp][i], mD2ResidDz2[i][j]);
        for (int k = N; k--;) {
          dchi2 += ROOT::Math::Dot(mDResidDz[k][i], mDResidDz[k][j]);
        }
      }
    }
  }
}
//___________________________________________________________________
template <int N, typename... Args>
void NA6PDCAFitterN<N, Args...>::calcPCA()
{
  // calculate point of closest approach for N prongs
  mPCA[mCurHyp] = mTrCFVT[mCurHyp][N - 1] * mTrPos[mCurHyp][N - 1];
  for (int i = N - 1; i--;) {
    mPCA[mCurHyp] += mTrCFVT[mCurHyp][i] * mTrPos[mCurHyp][i];
  }
}

//___________________________________________________________________
template <int N, typename... Args>
bool NA6PDCAFitterN<N, Args...>::recalculatePCAWithErrors(int cand)
{
  // recalculate PCA as a cov-matrix weighted mean, even if absDCA method was used
  if (isPropagateTracksToVertexDone(cand) && !propagateTracksToVertex(cand)) {
    return false;
  }
  int saveCurHyp = mCurHyp;
  mCurHyp = mOrder[cand];
  if (mUseAbsDCA) {
    for (int i = N; i--;) {
      if (!mTrcEInv[mCurHyp][i].set(mCandTr[mCurHyp][i], ZerrFactor)) { // prepare inverse cov.matrices at starting point
        mFitStatus[mCurHyp] = FitStatus::FailInvCov;
        if (mBadCovPolicy == Discard) {
          return false;
        } else if (mBadCovPolicy == OverrideAndFlag) {
          mPropFailed[mCurHyp] = true;
        } // otherwise, just use overridden errors w/o flagging
      }
    }
    if (!calcPCACoefs()) {
      mCurHyp = saveCurHyp;
      return false;
    }
  }
  auto oldPCA = mPCA[mOrder[cand]];
  calcPCA();
  mCurHyp = saveCurHyp;
  return true;
}

//___________________________________________________________________
template <int N, typename... Args>
void NA6PDCAFitterN<N, Args...>::calcPCANoErr()
{
  // calculate point of closest approach for N prongs w/o errors
  auto& pca = mPCA[mCurHyp];
  pca[0] = 0.;
  pca[1] = 0.;
  pca[2] = 0.;
  for (int i = 0; i < N; ++i) {
    pca[0] += mTrPos[mCurHyp][i][0];
    pca[1] += mTrPos[mCurHyp][i][1];
    pca[2] += mTrPos[mCurHyp][i][2];
  }
  pca[0] *= NInv;
  pca[1] *= NInv;
  pca[2] *= NInv;
}

//___________________________________________________________________
template <int N, typename... Args>
ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> NA6PDCAFitterN<N, Args...>::calcPCACovMatrix(int cand) const
{
  // calculate covariance matrix for the point of closest approach
  MatSym3D covmSum;
  for (int i = N; i--;) {
    const auto& tcov = mTrcEInv[mCurHyp][i];
    covmSum(0, 0) += tcov.sxx;
    covmSum(1, 0) += tcov.sxy;
    covmSum(1, 1) += tcov.syy;
    covmSum(2, 2) += tcov.szz;
  }
  if (covmSum.Invert()) {
    return covmSum;
  }
  // fall back on identity diagonal matrix with 1 cm error
  covmSum(0, 0) = 1.;
  covmSum(1, 0) = 0.;
  covmSum(1, 1) = 1.;
  covmSum(2, 2) = 1.;
  return covmSum;
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
double NA6PDCAFitterN<N, Args...>::calcChi2NoErr() const
{
  // calculate current chi2 of abs. distance minimization
  double chi2 = 0;
  for (int i = N; i--;) {
    const auto& res = mTrRes[mCurHyp][i];
    chi2 += res[0] * res[0] + res[1] * res[1] + res[2] * res[2];
  }
  return chi2;
}
//___________________________________________________________________
template <int N, typename... Args>
bool NA6PDCAFitterN<N, Args...>::correctTracks(const VecND& corrZ)
{
  // propagate tracks to updated Z
  for (int i = N; i--;) {
    const auto& trDer = mTrDer[mCurHyp][i];
    auto dz2h = 0.5 * corrZ[i] * corrZ[i];
    mTrPos[mCurHyp][i][0] -= trDer.dxdz * corrZ[i] - dz2h * trDer.d2xdz2;
    mTrPos[mCurHyp][i][1] -= trDer.dydz * corrZ[i] - dz2h * trDer.d2ydz2;
    mTrPos[mCurHyp][i][2] -= corrZ[i];
  }
  return true;
}
//___________________________________________________________________
template <int N, typename... Args>
bool NA6PDCAFitterN<N, Args...>::propagateTracksToVertex(int icand)
{
  // propagate tracks to current vertex
  int ord = mOrder[icand];
  if (mTrPropDone[ord]) {
    return true;
  }

  // need to refit taking as a seed already found vertex
  if (mRefitWithMatCorr) {
    int curHypSav = mCurHyp, curCrosIDAlt = mCrossIDAlt; // save
    mCurHyp = ord;
    mCrossIDAlt = -1; // disable alternative check
    auto restore = [this, curHypSav, curCrosIDAlt]() { this->mCurHyp = curHypSav; this->mCrossIDAlt = curCrosIDAlt; };
    if (!(mUseAbsDCA ? minimizeChi2NoErr() : minimizeChi2())) { // do final propagation
      restore();
      return false;
    }
    restore();
  }

  for (int i = N; i--;) {
    if (mUseAbsDCA || mUsePropagator /*|| mMatCorr != o2::base::Propagator::MatCorrType::USEMatCorrNONE*/) { // TODO: update after merging the recNative branch
      mCandTr[ord][i] = *mOrigTrPtr[i];                                                                      // fetch the track again, as mCandTr might have been propagated w/o errors or material corrections might be wrong
    }
    auto z = mPCA[ord][2]; // z of PCA
    if (!propagateToZ(mCandTr[ord][i], z)) {
      return false;
    }
  }

  mTrPropDone[ord] = true;
  return true;
}
//___________________________________________________________________
template <int N, typename... Args>
double NA6PDCAFitterN<N, Args...>::getAbsMax(const VecND& v)
{
  double mx = -1;
  for (int i = N; i--;) {
    auto vai = std::abs(v[i]);
    if (mx < vai) {
      mx = vai;
    }
  }
  return mx;
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
//___________________________________________________________________
template <int N, typename... Args>
bool NA6PDCAFitterN<N, Args...>::minimizeChi2NoErr()
{
  // find best chi2 (absolute DCA) of N tracks in the vicinity of the PCA seed

  for (int i = N; i--;) {
    mCandTr[mCurHyp][i] = *mOrigTrPtr[i];
    auto z = mPCA[mCurHyp][2];
    if (z > mMaxVertZ || z < mMinVertZ) {
      mFitStatus[mCurHyp] = FitStatus::RejTrackZ;
      return false;
    }
    if (!propagateParamToZ(mCandTr[mCurHyp][i], z)) {
      return false;
    }
    setTrackPos(mTrPos[mCurHyp][i], mCandTr[mCurHyp][i]); // prepare positions
  }
  if (mMaxDZIni > 0 && !roughDZCut()) { // apply rough cut on tracks Z difference
    mFitStatus[mCurHyp] = FitStatus::RejTrackRoughZ;
    return false;
  }

  calcPCANoErr();       // current PCA
  calcTrackResiduals(); // current track residuals
  float chi2Upd, chi2 = calcChi2NoErr();
  do {
    calcTrackDerivatives();      // current track derivatives (1st and 2nd)
    calcResidDerivativesNoErr(); // current residals derivatives (1st and 2nd)
    calcChi2DerivativesNoErr();  // current chi2 derivatives (1st and 2nd)

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
    calcPCANoErr(); // updated PCA
    if (mCrossIDAlt >= 0 && closerToAlternative()) {
      mFitStatus[mCurHyp] = FitStatus::FailCloserAlt;
      mAllowAltPreference = false;
      return false;
    }
    calcTrackResiduals();      // updated residuals
    chi2Upd = calcChi2NoErr(); // updated chi2
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
bool NA6PDCAFitterN<N, Args...>::closerToAlternative() const
{
  // check if the point current PCA point is closer to the seeding XY point being tested or to alternative see (if any)
  auto dxCur = mPCA[mCurHyp][0] - mCrossings.xDCA[mCrossIDCur], dzCur = mPCA[mCurHyp][2] - mCrossings.zDCA[mCrossIDCur];
  auto dxAlt = mPCA[mCurHyp][0] - mCrossings.xDCA[mCrossIDAlt], dzAlt = mPCA[mCurHyp][2] - mCrossings.zDCA[mCrossIDAlt];
  return dxCur * dxCur + dzCur * dzCur > dxAlt * dxAlt + dzAlt * dzAlt;
}

//___________________________________________________________________
template <int N, typename... Args>
NA6PTrack NA6PDCAFitterN<N, Args...>::createParentTrackParCov(int cand) const
{
  std::array<float, 21> covV = {0.};
  std::array<float, 3> pvecV = {0.};
  int q = 0;
  for (int it = 0; it < N; it++) {
    const auto& trc = getTrack(it, cand);
    auto pvecT = trc.getPXYZ();
    std::array<float, 21> covT = {0.};
    covT[0] = trc.getSigmaX2(); // TODO: update with new parameterization!
    covT[1] = trc.getSigmaXY();
    covT[2] = trc.getSigmaY2();
    covT[9] = trc.getSigmaPX2();
    covT[14] = trc.getSigmaPY2();
    covT[20] = trc.getSigmaPZ2();
    constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      covV[MomInd[i]] += covT[MomInd[i]];
    }
    for (int i = 0; i < 3; i++) {
      pvecV[i] += pvecT[i];
    }
    q += trc.getCharge();
  }
  auto covVtxV = calcPCACovMatrix(cand);
  covV[0] = covVtxV(0, 0);
  covV[1] = covVtxV(1, 0);
  covV[2] = covVtxV(1, 1);
  covV[3] = covVtxV(2, 0);
  covV[4] = covVtxV(2, 1);
  covV[5] = covVtxV(2, 2);
  NA6PTrack tr; // TODO: update adding covariance matrix
  tr.init(getPCACandidatePos(cand), pvecV, q);
  return tr;
}
//___________________________________________________________________
template <int N, typename... Args>
NA6PTrack NA6PDCAFitterN<N, Args...>::createParentTrackPar(int cand) const
{
  const auto& wvtx = getPCACandidate(cand);
  std::array<float, 3> pvecV = {0.};
  int q = 0;
  for (int it = 0; it < N; it++) {
    const auto& trc = getTrack(it, cand);
    auto pvecT = trc.getPXYZ();
    for (int i = 0; i < 3; i++) {
      pvecV[i] += pvecT[i];
    }
    q += trc.getCharge();
  }
  const std::array<float, 3> vertex = {(float)wvtx[0], (float)wvtx[1], (float)wvtx[2]};
  NA6PTrack tr;
  tr.init(vertex, pvecV, q);
  return tr;
}
//___________________________________________________________________
template <int N, typename... Args>
inline bool NA6PDCAFitterN<N, Args...>::propagateParamToZ(NA6PTrack& t, float z) const
{
  return t.propagateParamToZ(z, mBy);
}
//___________________________________________________________________
template <int N, typename... Args>
inline bool NA6PDCAFitterN<N, Args...>::propagateToZ(NA6PTrack& t, float z) const
{
  return t.propagateToZ(z, mBy);
}

using NA6PDCAFitter2 = NA6PDCAFitterN<2, NA6PTrack>;
using NA6PDCAFitter3 = NA6PDCAFitterN<3, NA6PTrack>;

#endif
