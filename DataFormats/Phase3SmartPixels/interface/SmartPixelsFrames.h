#ifndef DataFormats_Phase3SmartPixels_SmartPixelsFrames_h
#define DataFormats_Phase3SmartPixels_SmartPixelsFrames_h

// Module-frame <-> global-frame conversion for SmartPixels angles, in ONE place.
//
// WHY A SHARED HEADER. TBPX modules are TILTED: the angle between the module
// normal and the radial direction is on average 2.75/1.41/0.84/0.60 degrees on
// L1-L4 and reaches 16.5 degrees. So the local->global rotation is genuinely
// PER MODULE, not a per-layer constant, and a consumer that assumed otherwise
// would be wrong by up to 16 degrees with nothing to flag it. Any code needing a
// global-frame direction must use the module surface rotation, and there must be
// exactly one implementation of that.
//
// WHY cotAlpha/cotBeta ARE NOT GIVEN A "GLOBAL" VARIANT. They are DEFINED in the
// module frame (PixelAV convention: cotAlpha = p_x_local/p_z_local). A "global
// cotAlpha" is not a rotated version of the same quantity, it is a category
// error. The global frame has its own natural direction variables, and CMSSW
// already names them: phi and theta of the direction vector. So the global
// counterparts are dirPhi and dirCotTheta.

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

#include <cmath>

namespace smartpixels {

  struct GlobalDirection {
    float dirPhi = -999.f;       // azimuth of the direction, global frame [rad]
    float dirCotTheta = -999.f;  // pz/pt of the direction: z = z0 + r*cotTheta
    float sigDirPhi = -999.f;      // propagated from sigAlpha/sigBeta
    float sigDirCotTheta = -999.f;
    bool valid = false;
  };

  // Signed phi difference, wrapped to (-pi, pi]. Needed for the numerical
  // derivatives below: a raw subtraction across the +/-pi seam would give ~2*pi
  // and a wildly wrong uncertainty for every module near that seam.
  inline double dPhiWrap(double a, double b) {
    double d = a - b;
    while (d > M_PI)
      d -= 2 * M_PI;
    while (d <= -M_PI)
      d += 2 * M_PI;
    return d;
  }

  // Module-frame (cotAlpha, cotBeta) -> global direction, using the det's own
  // rotation. `toGlobalFn` must be the surface rotation of the SAME module the
  // angles were measured on; passing another module's is the failure this helper
  // exists to prevent.
  // sigAlpha/sigBeta are propagated to the global variables by a numerical
  // Jacobian, added in quadrature: the alpha and beta estimates come from separate
  // payload corrections and are treated as independent, which is the same
  // assumption the refit makes when it applies them as two independent scalar
  // Kalman updates. Pass negative sigmas to skip the propagation.
  template <typename Det>
  inline GlobalDirection toGlobalDirection(const Det& det, double cotAlpha, double cotBeta,
                                           double sigAlpha = -1., double sigBeta = -1.) {
    GlobalDirection out;
    if (!(std::isfinite(cotAlpha) && std::isfinite(cotBeta)))
      return out;
    // THE MODULE FRAME CARRIES NO DIRECTION SENSE, and this is the subtle part.
    // cotAlpha = lx/lz and cotBeta = ly/lz are BOTH invariant under v -> -v, so
    // (cotAlpha, cotBeta, 1) fixes a LINE, not a direction. Whether local +z (the
    // module NORMAL) points inward or outward is a placement detail that differs
    // module to module, so taking +1 blindly reversed the direction on roughly half
    // of all modules. Because a reversal flips phi by pi and negates cotTheta at the
    // same time, the symptom was: correct MAGNITUDES with scrambled SIGNS. Measured
    // before the fix -- 84.5% of tracks had sign-mixed cotTheta across their own
    // clusters, |cotTheta| spread within a track was only 0.0435 (magnitude fine),
    // and z - r*cotTheta gave 6.24 cm RMS against 3.47 cm for using z alone, i.e.
    // the correction was worse than doing nothing. Forcing the sign consistent took
    // it to 2.62 cm.
    //
    // Supply the sense from PHYSICS: a track crossing TBPX travels AWAY from the
    // beam axis, so the global direction must have a positive radial component at
    // the module. A very low-pT track curling back inward genuinely violates this,
    // and no rule can recover its sense from cotAlpha/cotBeta alone -- that
    // ambiguity is in the measurement, not in this code.
    const auto pos = det.position();
    const auto oriented = [&](double ca, double cb) {
      GlobalVector g = det.toGlobal(LocalVector(ca, cb, 1.0));
      if (g.x() * pos.x() + g.y() * pos.y() < 0.)
        g = GlobalVector(-g.x(), -g.y(), -g.z());
      return g;
    };
    const GlobalVector gv = oriented(cotAlpha, cotBeta);
    const double pt = std::hypot(gv.x(), gv.y());
    if (!(pt > 1e-12))
      return out;
    out.dirPhi = static_cast<float>(std::atan2(gv.y(), gv.x()));
    out.dirCotTheta = static_cast<float>(gv.z() / pt);
    out.valid = true;

    if (sigAlpha > 0. || sigBeta > 0.) {
      double vPhi = 0., vCot = 0.;
      const double base_phi = out.dirPhi, base_cot = out.dirCotTheta;
      for (int k = 0; k < 2; ++k) {
        const double sg = (k == 0) ? sigAlpha : sigBeta;
        if (!(sg > 0.))
          continue;
        // Step by the sigma itself rather than an arbitrary epsilon: the map is
        // smooth here, and this makes the linearization exact at the scale that
        // actually matters instead of at a scale nobody uses.
        // Must go through the SAME orientation rule as the nominal, or a module
        // whose normal points inward would give a perturbed vector pointing the
        // opposite way and a spurious ~pi difference in the numerical derivative.
        const GlobalVector gp = oriented(cotAlpha + (k == 0 ? sg : 0.), cotBeta + (k == 1 ? sg : 0.));
        const double ptp = std::hypot(gp.x(), gp.y());
        if (!(ptp > 1e-12))
          continue;
        const double dphi = dPhiWrap(std::atan2(gp.y(), gp.x()), base_phi);
        const double dcot = gp.z() / ptp - base_cot;
        vPhi += dphi * dphi;
        vCot += dcot * dcot;
      }
      out.sigDirPhi = static_cast<float>(std::sqrt(vPhi));
      out.sigDirCotTheta = static_cast<float>(std::sqrt(vCot));
    }
    return out;
  }

  // Inverse, used ONLY as a closure check: rotate the global direction back and
  // confirm the module-frame angles come out again. This catches a mis-applied
  // rotation, which on a 16-degree-tilted module would otherwise be a large silent
  // error.
  //
  // WHAT A ROUND-TRIP TEST CANNOT CATCH, learned the hard way. The angle test below
  // compares RATIOS -- lx/lz and ly/lz -- and both are invariant under v -> -v. A
  // direction that has been REVERSED therefore closes back perfectly. The earlier
  // version of this helper returned true on every module whose normal points inward
  // while the published phi was off by pi and cotTheta had the wrong sign. So the
  // orientation is now asserted SEPARATELY and FIRST. General lesson: "the inverse
  // transform reproduces the input" is evidence that the transform is invertible,
  // never that the convention is correct.
  template <typename Det>
  inline bool closesBackToModule(const Det& det, const GlobalDirection& g,
                                 double cotAlpha, double cotBeta, double tol = 1e-3) {
    if (!g.valid)
      return true;
    const double st = 1.0 / std::hypot(1.0, static_cast<double>(g.dirCotTheta));
    const GlobalVector gv(std::cos(g.dirPhi) * st,
                          std::sin(g.dirPhi) * st,
                          g.dirCotTheta * st);
    const auto pos = det.position();
    if (gv.x() * pos.x() + gv.y() * pos.y() < 0.)
      return false;  // reversed: invisible to the ratio test below, so checked here
    const LocalVector lv = det.toLocal(gv);
    if (!(std::abs(lv.z()) > 1e-12))
      return true;  // grazing in the module frame: the ratio is ill-conditioned
    return std::abs(lv.x() / lv.z() - cotAlpha) < tol && std::abs(lv.y() / lv.z() - cotBeta) < tol;
  }

}  // namespace smartpixels

#endif
