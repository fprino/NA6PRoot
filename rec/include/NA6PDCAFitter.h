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
