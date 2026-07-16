// NA6PCCopyright
#include <fairlogger/Logger.h>
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "Propagator.h"
#include "NA6PVertexerTracks.h"

ClassImp(NA6PVertexerTracks)

TrackVF::TrackVF(const NA6PTrackParCov& src, int id, float beamX, float beamY)
  : mX(src.getX()), mY(src.getY()), mZ(src.getZ()), trackIndex(id)
{
  // The track is already at its transverse DCA to the beam line. Estimate the
  // uncertainty of that Z PCA only for peak finding. The vertex fit itself
  // uses the two X/Y residuals and the covariance propagated to the seed Z.
  const double tx = src.getTx(), ty = src.getTy();
  const double den = tx * tx + ty * ty;
  const double cosPsi = src.getCosPsi();
  const double cosPsiI = 1. / cosPsi;
  const double dSlopeXdTx = cosPsiI * cosPsiI * cosPsiI;
  const double dSlopeYdTx = ty * tx * dSlopeXdTx;
  const double dSlopeYdTy = cosPsiI;
  mSlopeX = tx * cosPsiI;
  mSlopeY = ty * cosPsiI;

  // Precompute S(dz) = H(dz)*C*H(dz)^T as three quadratic polynomials.
  // Only the [x,y,tx,ty] covariance block contributes in the straight-line model.
  const double cXX = src.getCovMatElem(NA6PTrackPar::kX, NA6PTrackPar::kX);
  const double cXY = src.getCovMatElem(NA6PTrackPar::kX, NA6PTrackPar::kY);
  const double cYY = src.getCovMatElem(NA6PTrackPar::kY, NA6PTrackPar::kY);
  const double cXTx = src.getCovMatElem(NA6PTrackPar::kX, NA6PTrackPar::kTx);
  const double cXTy = src.getCovMatElem(NA6PTrackPar::kX, NA6PTrackPar::kTy);
  const double cYTx = src.getCovMatElem(NA6PTrackPar::kY, NA6PTrackPar::kTx);
  const double cYTy = src.getCovMatElem(NA6PTrackPar::kY, NA6PTrackPar::kTy);
  const double cTxTx = src.getCovMatElem(NA6PTrackPar::kTx, NA6PTrackPar::kTx);
  const double cTxTy = src.getCovMatElem(NA6PTrackPar::kTx, NA6PTrackPar::kTy);
  const double cTyTy = src.getCovMatElem(NA6PTrackPar::kTy, NA6PTrackPar::kTy);

  mSxx[0] = cXX;
  mSxx[1] = 2. * dSlopeXdTx * cXTx;
  mSxx[2] = dSlopeXdTx * dSlopeXdTx * cTxTx;

  mSxy[0] = cXY;
  mSxy[1] = dSlopeYdTx * cXTx + dSlopeYdTy * cXTy + dSlopeXdTx * cYTx;
  mSxy[2] = dSlopeXdTx * (dSlopeYdTx * cTxTx + dSlopeYdTy * cTxTy);

  mSyy[0] = cYY;
  mSyy[1] = 2. * (dSlopeYdTx * cYTx + dSlopeYdTy * cYTy);
  mSyy[2] = dSlopeYdTx * dSlopeYdTx * cTxTx + 2. * dSlopeYdTx * dSlopeYdTy * cTxTy + dSlopeYdTy * dSlopeYdTy * cTyTy;

  const double dx = beamX - src.getX(), dy = beamY - src.getY();
  const double num = tx * dx + ty * dy;
  const double denI = 1. / den, denI2 = denI * denI;

  // Gradient of zPCA = z + cosPsi*(tx*dx + ty*dy)/(tx^2 + ty^2)
  // with respect to the state [x, y, tx, ty, q/pXZ]. NA6PTrackPar::kQ2Pxz is not used in the calculation of zPCA, so its gradient is zero.
  double grad[4];
  grad[NA6PTrackPar::kX] = -cosPsi * tx * denI;
  grad[NA6PTrackPar::kY] = -cosPsi * ty * denI;
  grad[NA6PTrackPar::kTx] = -tx * num * denI / cosPsi + cosPsi * (dx * den - 2. * tx * num) * denI2;
  grad[NA6PTrackPar::kTy] = cosPsi * (dy * den - 2. * ty * num) * denI2;

  double sigmaZ2 = 0.;
  for (int i = 0; i < 4; ++i) {
    sigmaZ2 += grad[i] * grad[i] * src.getCovMatElem(i, i);
    for (int j = 0; j < i; ++j) {
      sigmaZ2 += 2. * grad[i] * grad[j] * src.getCovMatElem(i, j);
    }
  }
  mZSeedWeight = static_cast<float>(1. / sigmaZ2);
}

NA6PVertexerTracks::NA6PVertexerTracks()
{
  configurePeakFinding(mZMin, mZMax, mNBinsForPeakFind);
}

//___________________________________________________________________
void NA6PVertexerTracks::createTracksPool(const std::vector<NA6PTrack>& tracks)
{
  mTracksPool.clear();
  auto ntGlo = tracks.size();
  mTracksPool.reserve(ntGlo);
  auto prop = Propagator::Instance();
  for (uint32_t i = 0; i < ntGlo; i++) {
    NA6PTrackParCov trc = tracks[i];
    if (!prop->propagatePCAToLine(trc, mBeamX, mBeamY, mMaxDCA)) {
      continue;
    }
    auto& tvf = mTracksPool.emplace_back(trc, i, mBeamX, mBeamY);
    if (!tvf.isValid()) {
      mTracksPool.pop_back(); // discard bad track
      continue;
    }
  }
  if (mVerbose)
    LOGP(info, "Number of tracks in pool = {}", mTracksPool.size());
}

void NA6PVertexerTracks::buildAndFillHistoZ()
{
  std::fill(mHistZ.begin(), mHistZ.end(), 0.f);
  mFilledBinsZ.clear();
  for (const auto& tvf : mTracksPool) {
    float z = tvf.mZ;
    if (z >= mZMin && z < mZMax) {
      int bin = int((z - mZMin) / mZBinWidth);
      if (mHistZ[bin] == 0.f)
        mFilledBinsZ.push_back(bin);
      mHistZ[bin] += tvf.mZSeedWeight; // weighted fill
    }
  }
}

//___________________________________________________________________
int NA6PVertexerTracks::findPeakBin()
{
  if (mFilledBinsZ.empty())
    return -1;

  int maxBin = -1, ib = mFilledBinsZ.size(), last = ib;
  float maxv = 0.f;
  while (ib--) {
    auto bin = mFilledBinsZ[ib];
    auto v = mHistZ[bin];
    if (v > maxv) {
      maxv = v;
      maxBin = bin;
    } else if (v <= 0.f) {                     // bin was emptied
      mFilledBinsZ[ib] = mFilledBinsZ[--last]; // move last non-empty bin in place of emptied one
    }
  }
  mFilledBinsZ.resize(last);
  return maxBin;
}

//___________________________________________________________________
int NA6PVertexerTracks::findVertices(std::vector<NA6PVertex>& vertices)
{
  int nfound = 0;
  buildAndFillHistoZ();
  int nTrials = 0;
  while (nfound < mMaxVerticesPerCluster && nTrials < mMaxTrialsPerCluster) {
    int peakBin = findPeakBin();
    if (peakBin < 0) {
      break;
    }
    float zSeed = mZMin + (peakBin + 0.5f) * mZBinWidth;

    NA6PVertex vtx;
    if (fitVertex(zSeed, vtx)) {
      finalizeVertex(vtx, vertices);
      nfound++;
      nTrials = 0;
    } else {
      // suppress failed seeding bin and its proximities
      mHistZ[peakBin] = -1.f;
    }
    nTrials++;
  }
  return nfound;
}

//___________________________________________________________________
bool NA6PVertexerTracks::fitVertex(float zSeed, NA6PVertex& vtx)
{
  VertexSeed vtxSeed;
  vtxSeed.x = mBeamX;
  vtxSeed.y = mBeamY;
  vtxSeed.z = zSeed;
  vtxSeed.setScale(mInitScaleSigma2, mTukey2I);
  vtxSeed.scaleSigma2Prev = mInitScaleSigma2;
  vtx.setChi2(1.e30);

  FitStatus result = FitStatus::IterateFurther;
  while (result == FitStatus::IterateFurther) {
    vtxSeed.resetForNewIteration();
    vtxSeed.nIterations++;
    result = fitIteration(vtxSeed);

    if (result == FitStatus::OK) {
      result = evalIterations(vtxSeed);
    } else if (result == FitStatus::NotEnoughTracks) {
      if (vtxSeed.nIterations <= mMaxIterations && upscaleSigma(vtxSeed)) {
        result = FitStatus::IterateFurther;
        continue; // redo with stronger rescaling
      } else {
        break;
      }
    } else if (result == FitStatus::PoolEmpty || result == FitStatus::Failure) {
      break;
    }
  }

  if (result != FitStatus::OK) {
    vtx.setChi2(vtxSeed.maxScaleSigma2Tested);
    return false;
  }
  vtx.setXYZ(vtxSeed.x, vtxSeed.y, vtxSeed.z);
  vtx.setNContributors(vtxSeed.nContributors);
  vtx.setChi2(vtxSeed.getChi2());
  vtx.setCov(vtxSeed.cxx, vtxSeed.cxy, vtxSeed.cyy, vtxSeed.cxz, vtxSeed.cyz, vtxSeed.czz);
  vtx.setVertexType(NA6PVertex::kTrackPrimaryVertex);

  return true;
}

//___________________________________________________________________
void NA6PVertexerTracks::finalizeVertex(NA6PVertex& vtx, std::vector<NA6PVertex>& vertices)
{
  int lastID = vertices.size();
  for (auto& trc : mTracksPool) {
    if (trc.canAssign()) {
      vtx.addTrackID(trc.trackIndex);
      trc.vtxID = lastID;
    }
  }
  vertices.emplace_back(vtx);
}
//___________________________________________________________________
NA6PVertexerTracks::FitStatus NA6PVertexerTracks::fitIteration(VertexSeed& vtxSeed)
{
  int nTested = 0;
  for (auto& tvf : mTracksPool) {
    if (tvf.canUse()) {
      accountTrack(tvf, vtxSeed);
      nTested++;
    }
  }
  vtxSeed.maxScaleSigma2Tested = vtxSeed.scaleSigma2;
  if (vtxSeed.getNContributors() < mMinTracksPerVtx) {
    return nTested < mMinTracksPerVtx ? FitStatus::PoolEmpty : FitStatus::NotEnoughTracks;
  }
  if (!solveVertex(vtxSeed)) {
    return FitStatus::Failure;
  }
  return FitStatus::OK;
}

//___________________________________________________________________
void NA6PVertexerTracks::accountTrack(TrackVF& trc, VertexSeed& vtxSeed) const
{
  trc.wgh = 0.f;

  const double dz = vtxSeed.z - trc.mZ;
  const double predX = trc.mX + trc.mSlopeX * dz;
  const double predY = trc.mY + trc.mSlopeY * dz;
  const double sxx = std::fma(dz, std::fma(dz, trc.mSxx[2], trc.mSxx[1]), trc.mSxx[0]); // Evaluates Sxx[2]*dz^2 + Sxx[1]*dz + Sxx[0] using fused multiply-add
  const double sxy = std::fma(dz, std::fma(dz, trc.mSxy[2], trc.mSxy[1]), trc.mSxy[0]); // Evaluates Sxy[2]*dz^2 + Sxy[1]*dz + Sxy[0] using fused multiply-add
  const double syy = std::fma(dz, std::fma(dz, trc.mSyy[2], trc.mSyy[1]), trc.mSyy[0]); // Evaluates Syy[2]*dz^2 + Syy[1]*dz + Syy[0] using fused multiply-add
  const double slopeX = trc.mSlopeX;
  const double slopeY = trc.mSlopeY;

  // Residual convention: track prediction minus vertex position.
  double rx = predX - vtxSeed.x;
  double ry = predY - vtxSeed.y;
  double det = sxx * syy - sxy * sxy;
  if (det <= 1e-20f) {
    return;
  }

  double detI = 1. / det;
  double wxx = syy * detI;
  double wyy = sxx * detI;
  double wxy = -sxy * detI;

  constexpr float NDOF2I = 0.5f;
  double chi2T = rx * rx * wxx + 2. * rx * ry * wxy + ry * ry * wyy;
  chi2T *= NDOF2I; // chi2 per residual

  double wghT = 1. - chi2T * vtxSeed.scaleSig2ITuk2I;
  if (wghT < kAlmost0F) {
    return;
  }
  wghT *= wghT;

  double wrx = wxx * rx + wxy * ry;
  double wry = wxy * rx + wyy * ry;

  trc.wgh = wghT;
  vtxSeed.wghSum += wghT;
  vtxSeed.wghChi2 += wghT * chi2T;

  // D = d(track-vertex residual)/d(vx,vy,vz)
  //     = [-1,  0, slopeX]
  //       [ 0, -1, slopeY].
  // Accumulate A = D^T W D and b = -D^T W r, then solve A*delta=b.
  vtxSeed.cxx += wghT * wxx;
  vtxSeed.cyy += wghT * wyy;
  vtxSeed.cxy += wghT * wxy;
  vtxSeed.cxz -= wghT * (wxx * slopeX + wxy * slopeY);
  vtxSeed.cyz -= wghT * (wxy * slopeX + wyy * slopeY);
  vtxSeed.czz += wghT * (slopeX * slopeX * wxx + 2. * slopeX * slopeY * wxy + slopeY * slopeY * wyy);
  vtxSeed.cx0 += wghT * wrx;
  vtxSeed.cy0 += wghT * wry;
  vtxSeed.cz0 -= wghT * (slopeX * wrx + slopeY * wry);

  vtxSeed.nContributors++;
}

//___________________________________________________________________
bool NA6PVertexerTracks::solveVertex(VertexSeed& vtxSeed) const
{
  ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> mat;
  mat(0, 0) = vtxSeed.cxx;
  mat(0, 1) = vtxSeed.cxy;
  mat(0, 2) = vtxSeed.cxz;
  mat(1, 1) = vtxSeed.cyy;
  mat(1, 2) = vtxSeed.cyz;
  mat(2, 2) = vtxSeed.czz;
  if (!mat.InvertFast()) {
    LOGP(error, "Failed to invert matrix");
    return false;
  }
  ROOT::Math::SVector<double, 3> rhs(vtxSeed.cx0, vtxSeed.cy0, vtxSeed.cz0);
  auto sol = mat * rhs;
  // sol gives a correction to the current position, not the absolute position
  vtxSeed.x += sol(0);
  vtxSeed.y += sol(1);
  vtxSeed.z += sol(2);
  vtxSeed.cxx = mat(0, 0);
  vtxSeed.cxy = mat(1, 0);
  vtxSeed.cyy = mat(1, 1);
  vtxSeed.cxz = mat(2, 0);
  vtxSeed.cyz = mat(2, 1);
  vtxSeed.czz = mat(2, 2);

  vtxSeed.setChi2((vtxSeed.getNContributors() - vtxSeed.wghSum) / vtxSeed.scaleSig2ITuk2I); // calculate chi^2
  auto newScale = vtxSeed.wghSum > 0. ? vtxSeed.wghChi2 / vtxSeed.wghSum : mInitScaleSigma2;
  vtxSeed.setScale(newScale < mMinScale2 ? mMinScale2 : newScale, mTukey2I);
  return true;
}

//___________________________________________________________________
NA6PVertexerTracks::FitStatus NA6PVertexerTracks::evalIterations(VertexSeed& vtxSeed) const
{
  // decide if new iteration should be done, prepare next one if needed
  // if scaleSigma2 reached its lower limit stop
  FitStatus result = FitStatus::IterateFurther;
  if (mVerbose)
    LOGP(debug, "evalIteration: after {} iterations, scaleSigma2Prev = {}, scaleSigma2 = {}, nScaleIncrease = {}, nScaleSlowConvergence = {}, chi2 = {}", vtxSeed.nIterations, vtxSeed.scaleSigma2Prev, vtxSeed.scaleSigma2, vtxSeed.nScaleIncrease, vtxSeed.nScaleSlowConvergence, vtxSeed.getChi2() / vtxSeed.getNContributors());

  if (vtxSeed.nIterations > mMaxIterations) {
    result = FitStatus::Failure;
    return result;
  } else if (vtxSeed.scaleSigma2Prev <= mMinScale2 + kAlmost0F) {
    result = FitStatus::OK;
    // } else if (std::abs(vtxSeed.scaleSigma2 - vtxSeed.scaleSigma2Prev) < kScaleStability * vtxSeed.scaleSigma2Prev) {
    //   // scale is oscillating at a fixed point
    //   result = FitStatus::OK;
  }

  if (result == FitStatus::OK) {
    auto chi2Mean = vtxSeed.getChi2() / vtxSeed.getNContributors();
    if (chi2Mean > mMaxChi2Mean) {
      result = FitStatus::Failure;
    } else {
      return result;
    }
  }

  if (vtxSeed.scaleSigma2 > vtxSeed.scaleSigma2Prev) {
    if (++vtxSeed.nScaleIncrease > mMaxNScaleIncreased) {
      result = FitStatus::Failure;
    }
  } else if (vtxSeed.scaleSigma2 > mSlowConvergenceFactor * vtxSeed.scaleSigma2Prev) {
    if (++vtxSeed.nScaleSlowConvergence > mMaxNScaleSlowConvergence) {
      if (vtxSeed.scaleSigma2 < mAcceptableScale2) {
        vtxSeed.setScale(mMinScale2, mTukey2I);
        result = FitStatus::IterateFurther;
      } else {
        result = FitStatus::Failure;
      }
    }
  } else {
    vtxSeed.nScaleSlowConvergence = 0;
  }
  return result;
}

//___________________________________________________________________
bool NA6PVertexerTracks::upscaleSigma(VertexSeed& vtxSeed) const
{
  // scale upward the scaleSigma2 if needed
  if (vtxSeed.scaleSigma2 < mMaxScale2) {
    auto s = vtxSeed.scaleSigma2 * mUpscaleFactor;
    vtxSeed.setScale(s > mMaxScale2 ? mMaxScale2 : s, mTukey2I);
    return true;
  }
  return false;
}
