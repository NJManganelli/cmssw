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
    bool valid = false;
  };

  // Module-frame (cotAlpha, cotBeta) -> global direction, using the det's own
  // rotation. `toGlobalFn` must be the surface rotation of the SAME module the
  // angles were measured on; passing another module's is the failure this helper
  // exists to prevent.
  template <typename Det>
  inline GlobalDirection toGlobalDirection(const Det& det, double cotAlpha, double cotBeta) {
    GlobalDirection out;
    if (!(std::isfinite(cotAlpha) && std::isfinite(cotBeta)))
      return out;
    // Local direction with p_z_local = 1 by construction; magnitude is irrelevant
    // because only the direction is wanted.
    const LocalVector lv(cotAlpha, cotBeta, 1.0);
    const GlobalVector gv = det.toGlobal(lv);
    const double pt = std::hypot(gv.x(), gv.y());
    if (!(pt > 1e-12))
      return out;
    out.dirPhi = static_cast<float>(std::atan2(gv.y(), gv.x()));
    out.dirCotTheta = static_cast<float>(gv.z() / pt);
    out.valid = true;
    return out;
  }

  // Inverse, used ONLY as a closure check: rotate the global direction back and
  // confirm the module-frame angles come out again. This is what catches a
  // mis-applied rotation, which on a 16-degree-tilted module would otherwise be a
  // large silent error.
  template <typename Det>
  inline bool closesBackToModule(const Det& det, const GlobalDirection& g,
                                 double cotAlpha, double cotBeta, double tol = 1e-3) {
    if (!g.valid)
      return true;
    const double st = 1.0 / std::hypot(1.0, static_cast<double>(g.dirCotTheta));
    const GlobalVector gv(std::cos(g.dirPhi) * st,
                          std::sin(g.dirPhi) * st,
                          g.dirCotTheta * st);
    const LocalVector lv = det.toLocal(gv);
    if (!(std::abs(lv.z()) > 1e-12))
      return true;  // grazing in the module frame: the ratio is ill-conditioned
    return std::abs(lv.x() / lv.z() - cotAlpha) < tol && std::abs(lv.y() / lv.z() - cotBeta) < tol;
  }

}  // namespace smartpixels

#endif
