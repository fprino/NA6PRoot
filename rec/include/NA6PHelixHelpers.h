struct CircleXZ {

    float xC, zC, rC;   // centre and radius of helix projection in XZ bending plane

    CircleXZ() = default;
    CircleXZ(const NA6PTrack& tr, float by) {
      set(tr, by);
    }

    void set(const NA6PTrack& tr, float by)
    {
      std::array<float, 3> par = tr.getCircleParams(by);
      rC = float(par[0]);
      xC = float(par[1]);
      zC = float(par[2]);
    }
};

struct CrossInfo {
  
  static constexpr float MaxDistXZDef = 10.;

  float xDCA[2] = {};
  float zDCA[2] = {};
  int nDCA = 0;

  int circlesCrossInfo(const CircleXZ& trax0, const CircleXZ& trax1, float maxDistXZ = MaxDistXZDef, bool isCollinear = false)
  {
    const auto& trcA = trax0.rC > trax1.rC ? trax0 : trax1; // designate the largest circle as A
    const auto& trcB = trax0.rC > trax1.rC ? trax1 : trax0;
    nDCA = 0;
    float xDist = trcB.xC - trcA.xC, zDist = trcB.zC - trcA.zC;
    float dist2 = xDist * xDist + zDist * zDist, dist = std::sqrt(dist2), rsum = trcA.rC + trcB.rC;
    if (dist < 1e-12) {
      return nDCA; // circles are concentric?
    }
    if (dist > rsum) { // circles don't touch, chose a point in between
      // the parametric equation of lines connecting the centers is
      // x = x0 + t/dist * (x1-x0), z = z0 + t/dist * (z1-z0)
      if (dist - rsum > maxDistXZ) { // too large distance
        return nDCA;
      }
      notTouchingXZ(dist, xDist, zDist, trcA, trcB.rC, isCollinear);
    } else if (auto dfr = dist + trcB.rC - trcA.rC; dfr < 0.) { // the small circle is nestled into large one w/o touching
      if (dfr > -maxDistXZ) {
        // select the point of closest approach of 2 circles
        notTouchingXZ(dist, xDist, zDist, trcA, -trcB.rC, isCollinear);
      } else {
        return nDCA;
      }
    } else { // 2 intersection points
      if (isCollinear) {
        /// collinear tracks, e.g. electrons from photon conversion
        /// if there are 2 crossings of the circle it is better to take
        /// a weighted average of the crossing points as a radius
        float r2r = trcA.rC + trcB.rC;
        float r1_r = trcA.rC / r2r;
        float r2_r = trcB.rC / r2r;
        xDCA[0] = r2_r * trcA.xC + r1_r * trcB.xC;
        zDCA[0] = r2_r * trcA.zC + r1_r * trcB.zC;
        nDCA = 1;
      } else if (std::abs(xDist) < std::abs(zDist)) {
        // to simplify calculations, we move to new frame x->x+Xc0, z->z+Zc0, so that
        // the 1st one is centered in origin
        float a = (trcA.rC * trcA.rC - trcB.rC * trcB.rC + dist2) / (2. * zDist), b = -xDist / zDist, ab = a * b, bb = b * b;
        float det = ab * ab - (1. + bb) * (a * a - trcA.rC * trcA.rC);
        if (det > 0.) {
          det = std::sqrt(det);
          xDCA[0] = (-ab + det) / (1. + b * b);
          zDCA[0] = a + b * xDCA[0] + trcA.zC;
          xDCA[0] += trcA.xC;
          xDCA[1] = (-ab - det) / (1. + b * b);
          zDCA[1] = a + b * xDCA[1] + trcA.zC;
          xDCA[1] += trcA.xC;
          nDCA = 2;
        } else { // due to the finite precision the det<=0, i.e. the circles are barely touching, fall back to this special case
          notTouchingXZ(dist, xDist, zDist, trcA, trcB.rC);
        }
      } else {
        float a = (trcA.rC * trcA.rC - trcB.rC * trcB.rC + dist2) / (2. * xDist), b = -zDist / xDist, ab = a * b, bb = b * b;
        float det = ab * ab - (1. + bb) * (a * a - trcA.rC * trcA.rC);
        if (det > 0.) {
          det = std::sqrt(det);
          zDCA[0] = (-ab + det) / (1. + bb);
          xDCA[0] = a + b * zDCA[0] + trcA.xC;
          zDCA[0] += trcA.zC;
          zDCA[1] = (-ab - det) / (1. + bb);
          xDCA[1] = a + b * zDCA[1] + trcA.xC;
          zDCA[1] += trcA.zC;
          nDCA = 2;
        } else { // due to the finite precision the det<=0, i.e. the circles are barely touching, fall back to this special case
          notTouchingXZ(dist, xDist, zDist, trcA, trcB.rC);
        }
      }
    }
    return nDCA;
  }

  void notTouchingXZ(float dist, float xDist, float zDist, const CircleXZ& trcA, float rBSign, bool isCollinear = false)
  {
    if (isCollinear) {
      /// for collinear tracks it is better to take
      /// a weighted average of the crossing points as a radius
      float r2r = trcA.rC + rBSign;
      float r1_r = trcA.rC / r2r;
      float r2_r = rBSign / r2r;
      xDCA[0] = r2_r * trcA.xC + r1_r * (xDist + trcA.xC);
      zDCA[0] = r2_r * trcA.zC + r1_r * (zDist + trcA.zC);
    } else {
      // fast method to calculate DCA between 2 circles, assuming that they don't touch each outer:
      // the parametric equation of lines connecting the centers is x = xA + t/dist * xDist, z = zA + t/dist * zDist
      // with xA,zA being the center of the circle A ( = trcA.xC, trcA.zC ), xDist = trcB.xC = trcA.xC ...
      // There are 2 special cases:
      // (a) small circle is inside the large one: provide rBSign as -trcB.rC
      // (b) circle are side by side: provide rBSign as trcB.rC
      auto t2d = (dist + trcA.rC - rBSign) / dist;
      xDCA[0] = trcA.xC + 0.5 * (xDist * t2d);
      zDCA[0] = trcA.zC + 0.5 * (zDist * t2d);
    }
    nDCA = 1;
  }

  int linesCrossInfo(const NA6PTrack& tr0, const NA6PTrack& tr1, float maxDistXZ = MaxDistXZDef)
  {
    /// closest approach of 2 straight lines
    /// exploit NA6PLine
    nDCA = 0;
    NA6PLine line0 = NA6PLine::fromPointAndDirection(tr0.getXYZ().data(), tr0.getPXYZ().data());
    NA6PLine line1 = NA6PLine::fromPointAndDirection(tr1.getXYZ().data(), tr1.getPXYZ().data());
    
    float p1[3], p2[3];
    bool isOk = NA6PLine::getClosestPoints(line0, line1, p1, p2);
    if (!isOk)
      return nDCA;
    
    float dx = p1[0] - p2[0];
    float dz = p1[2] - p2[2];
    if (dx * dx + dz * dz > maxDistXZ * maxDistXZ) {
      return nDCA;   // lines too far apart at closest approach
    }
    xDCA[0] = 0.5f * (p1[0] + p2[0]);
    zDCA[0] = 0.5f * (p1[2] + p2[2]);
    nDCA = 1;

    return nDCA;
  }

  int circleLineCrossInfo(const CircleXZ& circ, const NA6PTrack& tr, float maxDistXZ = MaxDistXZDef)
  {
    /// closest approach of line and circle
    /// find intersection of a straight line with a circle in the XZ plane
    /// line: x = x0 + dx*t,  z = z0 + dz*t  (direction cosines dx, dz)
    /// circle: (x - xC)^2 + (z - zC)^2 = rC^2
    ///  (x0 - xC + dx*t)^2 + (z0 -zC + dz*t)^2 - rC^2 = 0
    ///  define: x0C = x0 - xC;  z0C = z0 - zC
    ///  x0C^2 + (dx*t)^2 + 2*x0C*dx*t + z0C^2 + (dz*t)^2 + 2*z0C*dz*t -rC^2 = 0
    ///  t^2 * (dx ^2 + dz^2) + t * (2*x0C*dx + 2*z0C*dz) + (x0C^2 + z0C^2 - rC^2) = 0
    ///  define: aA = (dx ^2 + dz^2); bB = x0C*dx + z0C*dz; cC = x0C^2 + z0C^2 - rC^2; det = bB^2 - aA*cC
    ///  t_1,2 = -bB / aA +- sqrt(det) / aA
    ///  xCross = x0 + t_1,2 * dx
    ///  zCross = z0 + t_1,2 * dz
    
    nDCA = 0;
    NA6PLine line = NA6PLine::fromPointAndDirection(tr.getXYZ().data(), tr.getPXYZ().data());
    float dx = line.mCosinesDirector[0];
    float dz = line.mCosinesDirector[2];
    float x0C = line.mOriginPoint[0] - circ.xC;
    float z0C = line.mOriginPoint[2] - circ.zC;
    float aA = dx * dx + dz * dz;
    float bB = (x0C * dx + z0C * dz);
    float cC = (x0C * x0C + z0C * z0C - circ.rC * circ.rC);
    float det = bB * bB - aA * cC;
    if (det >= 0.f) {
      float t1 = -bB / aA;
      float t2 = -std::sqrt(det) / aA;
      int nCand = (t2 < 1e-6f * std::abs(t1)) ? 1 : 2;
      float t[2] = {t1 - t2, t1 + t2};
      for (int i = 0; i < nCand; i++) {
        float xCross = line.mOriginPoint[0] + t[i] * dx;
        float zCross = line.mOriginPoint[2] + t[i] * dz;
        xDCA[nDCA] = xCross;
        zDCA[nDCA] = zCross;
        nDCA++;
      }
    } else {
      // there is no crossing, find the point of the closest approach on the line which is closest to the circle center
      float tClosest = -bB / aA;
      float xL = line.mOriginPoint[0] + tClosest * dx;
      float zL = line.mOriginPoint[2] + tClosest * dz;
      float dxc = xL - circ.xC;
      float dzc = zL - circ.zC;
      float dist = std::sqrt(dxc * dxc + dzc * dzc);
      if (dist - circ.rC > maxDistXZ) {
        return nDCA;
      }
      
      float drcf = circ.rC / dist; // radius / distance to circle center
      float xH = circ.xC + dxc * drcf;
      float zH = circ.zC + dzc * drcf;
      xDCA[0] = 0.5f * (xL + xH);
      zDCA[0] = 0.5f * (zL + zH);
      nDCA = 1;
    }

    return nDCA;
  }

  int set(const CircleXZ& trax0, const NA6PTrack& tr0, const CircleXZ& trax1, const NA6PTrack& tr1, float maxDistXZ = MaxDistXZDef, bool isCollinear = false)
  {
    // calculate up to 2 crossings between 2 tracks
    nDCA = 0;
    constexpr float Almost0 = 0x1.0p-126f;   // smallest non-denormal float
    bool curv0 = trax0.rC > Almost0;
    bool curv1 = trax1.rC > Almost0;

    if (curv0 && curv1) { // both are not straight lines
      nDCA = circlesCrossInfo(trax0, trax1, maxDistXZ, isCollinear);
    } else if (!curv0 && !curv1) { // both are straight lines
      nDCA = linesCrossInfo(tr0, tr1, maxDistXZ);
    } else {
      // curved track -> circle
      const auto& circ  = curv0 ? trax0  : trax1;
      // not curved track -> line
      const auto& trLin = curv0 ? tr1 : tr0;
      nDCA = circleLineCrossInfo(circ, trLin, maxDistXZ);
    }
    //
    return nDCA;
  }

  CrossInfo() = default;

  CrossInfo(const CircleXZ& trax0, const NA6PTrack& tr0, const CircleXZ& trax1, const NA6PTrack& tr1, float maxDistXZ = MaxDistXZDef, bool isCollinear = false)
  {
    set(trax0, tr0, trax1, tr1, maxDistXZ, isCollinear);
  }
  ClassDefNV(CrossInfo, 1);
};
