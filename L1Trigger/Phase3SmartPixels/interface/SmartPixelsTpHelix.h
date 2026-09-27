#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsTpHelix_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsTpHelix_h

// Pure helix of a TrackingParticle from its production state in a uniform field along
// z (bz at the origin, as SmartPixelsRecHitProducer's truth angle), and its crossing
// with a module mid-plane.
//
// helixPlaneCrossing is the "hx" truth of the cluster ntuple (hxDist, hxLocalX/Y, ...),
// and hxDist < 1 cm is what makes a TP-linked cluster's truth angle trustworthy for the
// emulator's class T: delta rays that Geant4 did not save inherit their parent's
// trackId, and a looper's second pass keeps the TP of the first, so their helix misses
// the hit by centimetres and its direction there is not the cluster's.

#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace smartpixels {

  // s = TRANSVERSE arc length from the production vertex [cm]
  struct TpHelix {
    bool curved = false;
    double R = 0, sgn = 0, cx = 0, cy = 0, th0 = 0;
    double vx = 0, vy = 0, vz = 0, ux = 0, uy = 0, cotTh = 0;

    TpHelix(const TrackingParticle& tp, double bz) {
      const double pt = std::hypot(tp.px(), tp.py());
      vx = tp.vx();
      vy = tp.vy();
      vz = tp.vz();
      ux = pt > 0 ? tp.px() / pt : 1.;
      uy = pt > 0 ? tp.py() / pt : 0.;
      cotTh = pt > 0 ? tp.pz() / pt : 0.;
      if (tp.charge() != 0 && pt > 1e-6 && std::abs(bz) > 1e-6) {
        curved = true;
        R = pt / (0.29979246 * std::abs(bz)) * 100.0;  // cm
        sgn = (tp.charge() * bz > 0.) ? +1.0 : -1.0;
        cx = vx + sgn * R * uy;
        cy = vy - sgn * R * ux;
        th0 = std::atan2(vy - cy, vx - cx);
      }
    }
    GlobalPoint pos(double s) const {
      if (!curved)
        return GlobalPoint(vx + s * ux, vy + s * uy, vz + s * cotTh);
      const double th = th0 - sgn * s / R;
      return GlobalPoint(cx + R * std::cos(th), cy + R * std::sin(th), vz + s * cotTh);
    }
    GlobalVector dir(double s) const {  // d(pos)/ds; transverse part has unit length
      if (!curved)
        return GlobalVector(ux, uy, cotTh);
      const double th = th0 - sgn * s / R;
      return GlobalVector(sgn * std::sin(th), -sgn * std::cos(th), cotTh);
    }
  };

  struct HelixCrossing {
    int nRoots = 0;
    double s = -1;                                         // < 0: the helix never crosses the plane
    double dist = std::numeric_limits<double>::infinity();  // |crossing - hit| [cm]
  };

  // Crossings of the helix with the module mid-plane over s in [0, min(2 turns, 300 cm)],
  // found by a scan in steps of `step` plus bisection; returns the one nearest `hit`.
  // stopBelow > 0 returns at the first crossing closer than it (enough to decide
  // dist < stopBelow; nRoots and dist are then those of the scan so far).
  inline HelixCrossing helixPlaneCrossing(
      const TpHelix& hx, const Plane& surf, const GlobalPoint& hit, double step, double stopBelow = -1.) {
    const GlobalPoint mc = surf.position();
    const GlobalVector n = surf.normalVector();
    auto f = [&](double s) {
      const GlobalPoint p = hx.pos(s);
      return n.x() * (p.x() - mc.x()) + n.y() * (p.y() - mc.y()) + n.z() * (p.z() - mc.z());
    };
    HelixCrossing out;
    const double sMax = hx.curved ? std::min(4. * M_PI * hx.R, 300.) : 300.;
    double s0 = 0, f0 = f(0);
    while (s0 < sMax) {
      const double s1 = std::min(s0 + step, sMax);
      const double f1 = f(s1);
      if ((f0 <= 0 && f1 > 0) || (f0 >= 0 && f1 < 0)) {
        double a = s0, b = s1, fa = f0;
        for (int it = 0; it < 50; ++it) {
          const double m = 0.5 * (a + b), fm = f(m);
          if ((fa <= 0 && fm > 0) || (fa >= 0 && fm < 0))
            b = m;
          else {
            a = m;
            fa = fm;
          }
        }
        const double sr = 0.5 * (a + b);
        ++out.nRoots;
        const double d = (hx.pos(sr) - hit).mag();
        if (d < out.dist) {
          out.dist = d;
          out.s = sr;
        }
        if (stopBelow > 0 && out.dist < stopBelow)
          return out;
      }
      s0 = s1;
      f0 = f1;
    }
    return out;
  }

}  // namespace smartpixels

#endif
