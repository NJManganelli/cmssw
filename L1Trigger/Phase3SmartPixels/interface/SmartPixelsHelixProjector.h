#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsHelixProjector_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsHelixProjector_h

// -*- C++ -*-
//
// Shared helix-propagation + Inner-Tracker BARREL (TBPX) module lookup used by
// BOTH the SmartPixelsPayloadAnalyzer (Phase 1) and the L1SmartPixelsTrackProducer
// digiRefit mode (Phase 2). ONE implementation, no drift.
//
// Carries the two hard-won Phase-1 fixes verbatim:
//   1. law-of-cosines uses |R| (unsigned) — a signed R flips cosArg for
//      negative-charge tracks and sent the crossing to the far side of the
//      circle, silently losing every negative track.
//   2. module containment on staggered/tilted TBPX planks: pick the nearest
//      module center, then require the crossing to be inside the module's
//      TRANSVERSE (local x,y) bounds ONLY (ignore thickness/local-z) — the
//      correct test for a thin sensor.
//
// Propagation is intentionally LIGHT: analytic helix -> mean-radius cylinder
// intersection in the transverse plane, then module lookup + local-frame
// transform. No reco propagators, no material effects (documented, adequate for
// barrel window sizing / interim refit). v1 upgrades (k=2-nearest modules, plane
// re-propagation, endcap TEPX/TFPX) land HERE once, for both consumers.

#include <vector>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

class TrackerGeometry;
class TrackerTopology;
class MagneticField;
class GeomDet;
class PixelGeomDetUnit;

namespace smartpixels {

  // Helix parameters in the CMS convention, all extracted from a TTTrack:
  //   rInv  signed 1/R [cm^-1]
  //   phi0  momentum phi at POCA
  //   x0,y0 POCA transverse position [cm]
  //   z0    longitudinal POCA [cm]
  //   tanL  tan(lambda) = pz/pt
  //   pt    transverse momentum [GeV] (only sets the momentum-vector scale;
  //         angles are scale-invariant)
  struct HelixParams {
    double rInv = 0.;
    double phi0 = 0.;
    double x0 = 0.;
    double y0 = 0.;
    double z0 = 0.;
    double tanL = 0.;
    double pt = 1.;
  };

  // Result of propagating a helix to one TBPX layer and landing it on a module.
  struct Crossing {
    bool valid = false;
    int layer = 0;                       // 1-based TBPX layer
    const PixelGeomDetUnit* det = nullptr;
    unsigned int detId = 0;
    GlobalPoint global;                  // crossing point (global) [cm]
    LocalPoint local;                    // crossing point in module-local frame [cm]
    GlobalVector globalMom;              // track momentum direction at crossing (global)
    LocalVector localMom;                // ... in module-local frame
    double cotAlpha = 0.;                // PixelAV convention: p_x_local / p_z_local
    double cotBeta = 0.;                 //                     p_y_local / p_z_local
    double bLocalX = 0., bLocalY = 0., bLocalZ = 0.;  // local B-field [T]
  };

  class HelixProjector {
  public:
    HelixProjector() = default;

    // Build the per-layer module cache from geometry. Idempotent for a fixed
    // (geometry pointer, nLayers): rebuilds only if those change. Safe to call
    // every event. nLayers clamps to [1,4] (TBPX barrel layers).
    void build(const TrackerGeometry& geom, const TrackerTopology& topo, int nLayers);

    bool ready() const { return built_ && geomTag_ != nullptr; }
    int nLayers() const { return static_cast<int>(layerModules_.size()); }
    double layerRadius(int layer1based) const { return layerRadius_.at(layer1based - 1); }

    // Propagate the helix to TBPX layer `layer1based` (1..nLayers) and find the
    // module it lands on. Returns Crossing{valid=false} if the helix does not
    // reach the layer radius or lands off every module. `field` supplies the
    // local B at the module center.
    Crossing crossLayer(const HelixParams& helix, int layer1based, const MagneticField& field) const;

  private:
    bool built_ = false;
    const TrackerGeometry* geomTag_ = nullptr;  // identity of the cached geometry
    int nLayersReq_ = 0;
    std::vector<std::vector<const GeomDet*>> layerModules_;  // [layer-1] -> units (GeomDetUnit)
    std::vector<double> layerRadius_;                            // [layer-1] mean r [cm]
  };

}  // namespace smartpixels

#endif
