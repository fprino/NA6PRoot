// NA6PCCopyright

#include <TMath.h>
#include <TRandom.h>
#include <fairlogger/Logger.h>
#include "NA6PGlauber.h"

NA6PGlauber::NA6PGlauber()
{
  mNuclA.setNucleusPb();
  mNuclB.setNucleusPb();
  mSigmamb = 30.0;
}

NA6PGlauber::NA6PGlauber(const NA6PNuclearShape& nucA, const NA6PNuclearShape& nucB, double sigmb) : mNuclA(nucA),
                                                                                                     mNuclB(nucB),
                                                                                                     mSigmamb(sigmb)
{
}

void NA6PGlauber::setIntegrationSteps(double dx, double dy, double dz)
{
  mIntegStep[0] = dx;
  mIntegStep[1] = dy;
  mIntegStep[2] = dz;
}

double NA6PGlauber::getEccentricity() const
{
  if (mNuclA.getA() == 1 || mNuclB.getA() == 1 || mTaTb <= 0.)
    return -999.;
  const double varX = getVarX();
  const double varY = getVarY();
  const double sum = varX + varY;
  return sum > 0. ? (varX - varY) / sum : -999.;
}

void NA6PGlauber::doOpticalIntegrals(double b)
{
  mImpPar = b;
  mPartDens = 0.;
  mPartDensMax = 0.;
  mCollDens = 0.;
  mCollDensMax = 0.;
  mAverNPart = 0.;
  mTAB = 0.;
  mTaTb = 0.;
  mRhoL = 0.;
  mTa = 0.;
  mX = 0.;
  mY = 0.;
  mXX = 0.;
  mYY = 0.;
  mTb = 0.;
  double x, y;
  double dx = mIntegStep[0];
  double dy = mIntegStep[1];
  double dz = mIntegStep[2];
  double dxdy = dx * dy;
  int APro = mNuclA.getA();
  int ATgt = mNuclB.getA();
  if (APro == 1 && ATgt == 1) {
    LOGP(error, "pp case not treated");
    return;
  } else if (APro == 1 || ATgt == 1) {
    doOpticalIntegralsPA(b);
    return;
  }
  double maxrad = mNuclB.getMaxRadius();
  double xmin = -maxrad;
  double ymin = -maxrad;
  int nptx = (2. * maxrad) / dx + 1.00001;
  int npty = (2. * maxrad) / dy + 1.00001;

  for (int ix = 0; ix < nptx; ++ix) {
    x = xmin + ix * dx;
    for (int iy = 0; iy < npty; ++iy) {
      y = ymin + iy * dy;
      double tA = mNuclA.calcSurfNuclDens(x, y - 0.5 * mImpPar, dz);
      double tB = mNuclB.calcSurfNuclDens(x, y + 0.5 * mImpPar, dz);
      mTa += tA * dxdy;
      mTb += tB * dxdy;
      mTAB += tA * tB * dxdy;
      double pard = calcPartDens(tA, tB);
      double colld = calcCollDens(tA, tB);
      if (pard > mPartDensMax)
        mPartDensMax = pard;
      if (colld > mCollDensMax)
        mCollDensMax = colld;
      mPartDens += pard * tA * tB * dxdy;
      mCollDens += colld * tA * tB * dxdy;
      mX += x * tA * tB * dxdy;
      mY += y * tA * tB * dxdy;
      mXX += x * x * tA * tB * dxdy;
      mYY += y * y * tA * tB * dxdy;
      mAverNPart += pard * dxdy;
      mTaTb += tA * tB * dxdy;
      mRhoL += tA * tB * 0.5 * ((APro - 1.) * tA + (ATgt - 1.) * tB) * dxdy;
    }
  }
}

void NA6PGlauber::doOpticalIntegralsPA(double b)
{
  bool protonIsA = (mNuclA.getA() == 1);
  const NA6PNuclearShape& nuc = protonIsA ? mNuclB : mNuclA;
  int An = nuc.getA();
  double sigmaNN = getSigmafm2();
  double tN = nuc.calcSurfNuclDens(b, 0., mIntegStep[2]);
  double pInel = -std::expm1(An * std::log1p(-sigmaNN * tN));
  double nColl = An * sigmaNN * tN;
  double nPart = pInel + nColl;
  double yp = protonIsA ? 0.5 * b : -0.5 * b;
  mTa = 1.;
  mTb = 1.;
  mTAB = tN;
  mTaTb = tN;
  mCollDens = 0.;
  mCollDensMax = 0.;
  mAverNPart = nPart;
  mPartDens = 0.;
  mPartDensMax = 0.;
  mX = 0.;
  mY = yp * tN;
  mXX = 0.;
  mYY = yp * yp * tN;
  mRhoL = 0.5 * (An - 1.) * tN * tN;
}

double NA6PGlauber::calcPartDens(double tA, double tB) const
{
  int APro = mNuclA.getA();
  int ATgt = mNuclB.getA();
  double sigmaA, sigmaB;
  double crossec = getSigmafm2();
  sigmaA = -std::expm1(APro * std::log1p(-crossec * tA));
  sigmaB = -std::expm1(ATgt * std::log1p(-crossec * tB));
  double pard = APro * tA * sigmaB + ATgt * tB * sigmaA;
  return pard;
}

double NA6PGlauber::calcPartDensFromCoord(double x, double y, double b) const
{
  double pard = 0;
  double dz = mIntegStep[2];
  double tA = mNuclA.calcSurfNuclDens(x, y - 0.5 * b, dz);
  double tB = mNuclB.calcSurfNuclDens(x, y + 0.5 * b, dz);
  pard = calcPartDens(tA, tB);
  return pard;
}

double NA6PGlauber::calcCollDens(double tA, double tB) const
{
  int APro = mNuclA.getA();
  int ATgt = mNuclB.getA();
  double crossec = getSigmafm2();
  double colld = APro * ATgt * tA * tB * crossec;
  return colld;
}

double NA6PGlauber::calcCollDensFromCoord(double x, double y, double b) const
{
  double colld = 0;
  double dz = mIntegStep[2];
  double tA = mNuclA.calcSurfNuclDens(x, y - 0.5 * b, dz);
  double tB = mNuclB.calcSurfNuclDens(x, y + 0.5 * b, dz);
  colld = mNuclA.getA() * mNuclB.getA() * tA * tB * getSigmafm2();
  return colld;
}

bool NA6PGlauber::isInside(double x, double y, double b) const
{
  double a1, a2;
  double rPro = mNuclA.getMaxRadius();
  double rTgt = mNuclB.getMaxRadius();
  a1 = rTgt * rTgt - x * x - (y + 0.5 * b) * (y + 0.5 * b);
  a2 = rPro * rPro - x * x - (y - 0.5 * b) * (y - 0.5 * b);
  if (a1 > 0 && a2 > 0)
    return true;
  return false;
}

NA6PGlauber::GlauberEvent NA6PGlauber::generateEvent(double bmin, double bmax, bool storePositions)
{
  double r01 = gRandom->Rndm();
  if (bmin < 0) {
    bmin = 0.;
    bmax = mNuclA.getRadius() + 5. * mNuclA.getSkinDepth() + mNuclB.getRadius() + 5. * mNuclB.getSkinDepth();
  }
  double b = bmin;
  if (bmax > bmin)
    b = std::sqrt(bmin * bmin + r01 * (bmax * bmax - bmin * bmin));
  double phib = gRandom->Rndm() * 2. * M_PI;
  double impparx = b * std::cos(phib);
  double imppary = b * std::sin(phib);
  GlauberEvent event;
  event.impactParameter = b;
  event.phiImpPar = phib;

  const int A = mNuclA.getA();
  const int B = mNuclB.getA();
  mNuclA.setNucleonPositions(0., 0.);
  mNuclB.setNucleonPositions(impparx, imppary);
  std::span<const double> xA = mNuclA.getXpos();
  std::span<const double> yA = mNuclA.getYpos();
  std::span<const double> xB = mNuclB.getXpos();
  std::span<const double> yB = mNuclB.getYpos();
  if (static_cast<int>(xA.size()) != A ||
      static_cast<int>(xB.size()) != B ||
      static_cast<int>(yA.size()) != A ||
      static_cast<int>(yB.size()) != B) {
    LOGP(error, "mismatch in array dimensions");
  }

  double interactionRadiusSquare = getSigmafm2() / M_PI;
  std::array<int, NA6PNuclearShape::kMaxNucleons> isNucleonStruckA;
  std::array<int, NA6PNuclearShape::kMaxNucleons> isNucleonStruckB;
  std::fill(isNucleonStruckA.begin(), isNucleonStruckA.begin() + A, 0);
  std::fill(isNucleonStruckB.begin(), isNucleonStruckB.begin() + B, 0);

  int nColl = 0;
  for (int ia = 0; ia < A; ++ia) {
    for (int ib = 0; ib < B; ++ib) {
      double dx = xA[ia] - xB[ib];
      double dy = yA[ia] - yB[ib];
      double distSquare = dx * dx + dy * dy;
      if (distSquare <= interactionRadiusSquare) {
        ++nColl;
        isNucleonStruckA[ia] += 1;
        isNucleonStruckB[ib] += 1;
        if (storePositions) {
          event.xColl.push_back(0.5 * (xA[ia] + xB[ib]));
          event.yColl.push_back(0.5 * (yA[ia] + yB[ib]));
        }
      }
    }
  }
  int nPartA = 0, nPartB = 0;
  double sumx = 0., sumx2 = 0., sumy = 0., sumy2 = 0., sumxy = 0.;
  for (int ia = 0; ia < A; ++ia) {
    if (isNucleonStruckA[ia] > 0) {
      ++nPartA;
      sumx += xA[ia];
      sumx2 += (xA[ia] * xA[ia]);
      sumy += yA[ia];
      sumy2 += (yA[ia] * yA[ia]);
      sumxy += (xA[ia] * yA[ia]);
      if (storePositions) {
        event.xPart.push_back(xA[ia]);
        event.yPart.push_back(yA[ia]);
      }
    }
  }
  for (int ib = 0; ib < B; ++ib) {
    if (isNucleonStruckB[ib] > 0) {
      ++nPartB;
      sumx += xB[ib];
      sumx2 += (xB[ib] * xB[ib]);
      sumy += yB[ib];
      sumy2 += (yB[ib] * yB[ib]);
      sumxy += (xB[ib] * yB[ib]);
      if (storePositions) {
        event.xPart.push_back(xB[ib]);
        event.yPart.push_back(yB[ib]);
      }
    }
  }
  int nPart = nPartA + nPartB;
  if (nPart > 0) {
    sumx /= static_cast<double>(nPart);
    sumx2 /= static_cast<double>(nPart);
    sumy /= static_cast<double>(nPart);
    sumy2 /= static_cast<double>(nPart);
    sumxy /= static_cast<double>(nPart);
  }
  double sigmax2 = sumx2 - sumx * sumx;
  double sigmay2 = sumy2 - sumy * sumy;
  double sigmaxy = sumxy - sumx * sumy;
  double ecc = -999.;
  if ((sigmax2 + sigmay2) > 0.)
    ecc = std::sqrt((sigmay2 - sigmax2) * (sigmay2 - sigmax2) + 4.0 * sigmaxy * sigmaxy) / (sigmax2 + sigmay2);

  event.nParticipants = nPart;
  event.nSpectatorsA = A - nPartA;
  event.nSpectatorsB = B - nPartB;
  event.nCollisions = nColl;
  event.eccentricity = ecc;
  return event;
}
