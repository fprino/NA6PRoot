// NA6PCCopyright

#include <TRandom.h>
#include <TF1.h>
#include <fairlogger/Logger.h>
#include "NA6PNuclearShape.h"

NA6PNuclearShape::NA6PNuclearShape()
{
  setNucleusPb();
}

NA6PNuclearShape::NA6PNuclearShape(int Z, int A)
{
  setNucleusFromParametrization(Z, A);
}

NA6PNuclearShape::NA6PNuclearShape(int Z, int A, double r0, double wid, double w, int e)
{
  mZ = Z;
  setA(A);
  setShapeParameters(r0, wid, w, e);
}

void NA6PNuclearShape::setNucleusFromParametrization(int Z, int A)
{
  // 2pF with parametrized radius and skin depth
  mZ = Z;
  setA(A);
  double r = 1.12 * std::pow(mA, 1. / 3.) - 0.86 * std::pow(mA, -1. / 3.);
  setShapeParameters(r, 0.54, 0., 1);
}

void NA6PNuclearShape::setA(int A)
{
  if (A < 1 || A > kMaxNucleons) {
    LOGP(error, "A = {} outside [1, {}], clamping", A, kMaxNucleons);
    A = std::clamp(A, 1, kMaxNucleons);
  }
  mA = A;
  normalizeRho0();
  resetNucleonPositions();
}

void NA6PNuclearShape::setShapeParameters(double r0, double wid, double w, int e)
{
  mRadius = r0;
  mSkinDepth = wid;
  mW = w;
  mExpon = e;
  normalizeRho0();
  resetNucleonPositions();
}

void NA6PNuclearShape::setNucleusC()
{
  mZ = 6;
  mA = 12;
  setShapeParameters(2.355, 0.5224, -0.149, 1);
}

void NA6PNuclearShape::setNucleusO()
{
  mZ = 8;
  mA = 16;
  setShapeParameters(2.608, 0.513, -0.051, 1);
}

void NA6PNuclearShape::setNucleusS()
{
  mZ = 16;
  mA = 32;
  setShapeParameters(3.503, 0.633, -0.250, 1);
}

void NA6PNuclearShape::setNucleusAr2pF()
{
  mZ = 18;
  mA = 40;
  setShapeParameters(3.534, 0.542, 0., 1);
}

void NA6PNuclearShape::setNucleusAr3pF()
{
  mZ = 18;
  mA = 40;
  setShapeParameters(3.735, 0.621, -0.194, 1);
}

void NA6PNuclearShape::setNucleusCu2pF()
{
  mZ = 29;
  mA = 63;
  setShapeParameters(4.218, 0.596, 0., 1);
}

void NA6PNuclearShape::setNucleusCu3pG()
{
  mZ = 29;
  mA = 63;
  setShapeParameters(3.707, 2.412, 0.339, 2);
}

void NA6PNuclearShape::setNucleusIn()
{
  mZ = 49;
  mA = 115;
  setShapeParameters(5.357, 0.563, 0., 1);
}

void NA6PNuclearShape::setNucleusAu()
{
  mZ = 79;
  mA = 197;
  setShapeParameters(6.386, 0.535, 0., 1);
}

void NA6PNuclearShape::setNucleusPb()
{
  mZ = 82;
  mA = 208;
  setShapeParameters(6.624, 0.549, 0., 1);
}

void NA6PNuclearShape::setNucleusU()
{
  mZ = 92;
  mA = 238;
  setShapeParameters(6.874, 0.556, 0., 1);
}

double NA6PNuclearShape::integrateRho() const
{
  double integral = 0.;
  for (int jStep = 0; jStep < 1000; jStep++) {
    double r = kRMax / 1000. * jStep;
    integral += radialDensity(r) * kRMax / 1000. * 4.0 * M_PI;
  }
  return integral;
}

double NA6PNuclearShape::calcRho(double r) const
{
  double rho = 0.;
  if (mExpon == 1) {
    rho = mRho0 * (1. + mW * r * r / mRadius / mRadius) / ((exp((r - mRadius) / mSkinDepth) + 1.));
  } else if (mExpon == 2) {
    rho = mRho0 * (1 + mW * r * r / mRadius / mRadius) / ((exp((r * r - mRadius * mRadius) / mSkinDepth / mSkinDepth) + 1.));
  }
  return rho;
}

double NA6PNuclearShape::computeMaxRadius() const
{
  double ref = calcRho(mRadius);
  double r = kRMax;
  while (r > mRadius && std::abs(calcRho(r)) <= kPrecision * ref)
    r -= 0.05;
  return r + 0.05;
}

void NA6PNuclearShape::normalizeRho0()
{
  mRho0 = 1.;
  mRho0 = static_cast<double>(mA) / integrateRho();
  mRadDensMax = estimateRadialMax(kRMin, kRMax);
  mMaxRadius = computeMaxRadius();
}
double NA6PNuclearShape::estimateRadialMax(double rMin, double rMax, int samples)
{
  double maxVal = 0.;
  for (int i = 0; i <= samples; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(samples);
    const double r = rMin + t * (rMax - rMin);
    maxVal = std::max(maxVal, radialDensity(r));
  }
  return maxVal * 1.1;
}

double NA6PNuclearShape::sampleRadius(double rhoMax, double rMin, double rMax) const
{
  while (true) {
    const double r = rMin + static_cast<double>(gRandom->Rndm()) * (rMax - rMin);
    const double u = static_cast<double>(gRandom->Rndm()) * rhoMax;
    if (u <= radialDensity(r)) {
      return r;
    }
  }
}

double NA6PNuclearShape::calcSurfNuclDens(double x, double y, double zstep) const
{
  double zmax = mMaxRadius;
  double tZ = calcSurfNuclDens(x, y, zstep, zmax);
  return tZ;
}

double NA6PNuclearShape::calcSurfNuclDens(double x, double y, double zstep, double zm) const
{
  if (x * x + y * y >= mMaxRadius * mMaxRadius)
    return 0.;
  double tZ = 0.;
  double zmax = zm;
  double zmin = -zmax;
  int nz = static_cast<int>((zmax - zmin) / zstep + 1.0001);
  for (int i = 0; i <= nz; ++i) {
    double z = zmin + i * zstep;
    double r = std::sqrt(x * x + y * y + z * z);
    tZ += calcRho(r);
  }
  tZ *= zstep / static_cast<double>(mA);
  return tZ;
}

TF1 NA6PNuclearShape::getDensity() const
{
  const char* formula = (mExpon == 1) ? "[0]*(1.+[1]*x*x/[2]/[2])/(exp((x-[2])/[3])+1.)" : "[0]*(1+[1]*x*x/[2]/[2])/(exp((x*x-[2]*[2])/[3]/[3])+1.)";
  TF1 rho("rho", formula, kRMin, kRMax);
  rho.SetParameters(mRho0, mW, mRadius, mSkinDepth);
  return rho;
}

void NA6PNuclearShape::setNucleonPositions(double shiftx, double shifty)
{
  if (mA == 1) {
    mXpos[0] = shiftx;
    mYpos[0] = shifty;
    mZpos[0] = 0;
    return;
  }
  double dmin2 = 0.;
  if (mMinNucleonDist > 0)
    dmin2 = mMinNucleonDist * mMinNucleonDist;
  for (int i = 0; i < mA; i++) {
    int nAttempts = 0;
    double x, y, z;
    do {
      double r = sampleRadius(mRadDensMax, kRMin, mMaxRadius);
      double phi = gRandom->Rndm() * 2. * M_PI;
      double costh = gRandom->Rndm() * 2. - 1.;
      double sinth = std::sqrt(std::max(0., 1. - costh * costh));
      x = r * sinth * std::cos(phi);
      y = r * sinth * std::sin(phi);
      z = r * costh;
    } while (dmin2 > 0 && ++nAttempts < kMaxAttempts && !isFarEnough(i, x, y, z, dmin2));
    mXpos[i] = x;
    mYpos[i] = y;
    mZpos[i] = z;
  }
  for (int i = 0; i < mA; ++i) {
    mXpos[i] += shiftx;
    mYpos[i] += shifty;
  }
}

bool NA6PNuclearShape::isFarEnough(int n, double x, double y, double z, double d2) const
{
  for (int j = 0; j < n; ++j) {
    double dx = x - mXpos[j];
    double dy = y - mYpos[j];
    double dz = z - mZpos[j];
    if (dx * dx + dy * dy + dz * dz < d2)
      return false;
  }
  return true;
}

double NA6PNuclearShape::calcSurfNuclDens2(double x, double y, double zstep, double zm) const
{
  double tZ = 0;
  double zmax = zm;
  double zmin = -zmax;
  double r1;
  double r2;
  double z1 = zmin;
  double z2 = zmin + zstep;
  do {
    r1 = sqrt(x * x + y * y + z1 * z1);
    r2 = sqrt(x * x + y * y + z2 * z2);
    tZ += zstep * (calcRho(r1) + calcRho(r2)) / 2. / double(mA);
    z1 += zstep;
    z2 = z1 + zstep;
  } while (z1 < zmax);

  return tZ;
}
