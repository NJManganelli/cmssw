#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsHelixProjector_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsHelixProjector_h

// -*- C++ -*-
//
// Shared helix-propagation + Inner-Tracker module lookup (TBPX layers, TFPX discs) used by
// BOTH the SmartPixelsPayloadAnalyzer (Phase 1) and the L1SmartPixelsTrackProducer
// seedlessRefit mode (Phase 2). ONE implementation, no drift.
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
// Propagation is intentionally LIGHT: analytic helix -> mean-radius cylinder (barrel) or
// mean-|z| plane (disc) intersection, then module lookup (plane re-propagation onto the
// candidate modules) + local-frame transform. No reco propagators, no material effects.

#include <array>
#include <utility>
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

  // Result of propagating a helix to one IT surface and landing it on a module.
  struct Crossing {
    bool valid = false;
    int layer = 0;                       // IT surface code (SmartPixelsSurface.h): 1-4 TBPX, 11-14 ID1-ID4
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

    // Build the per-layer (and per-disc) module cache from geometry. Idempotent for a
    // fixed (geometry pointer, nLayers, nDiscs): rebuilds only if those change. Safe to
    // call every event. nLayers clamps to [1,4] (TBPX barrel layers), nDiscs to [0,4]
    // (TFPX disks 1..nDiscs, both sides = ID1-ID4).
    void build(const TrackerGeometry& geom, const TrackerTopology& topo, int nLayers, int nDiscs = 0);

    bool ready() const { return built_ && geomTag_ != nullptr; }
    int nLayers() const { return static_cast<int>(layerModules_.size()); }
    int nDiscs() const { return static_cast<int>(discZ_.size()); }
    double layerRadius(int layer1based) const { return layerRadius_.at(layer1based - 1); }

    // Transverse arc length from the POCA at which the helix meets the nominal surface
    // (layer cylinder at the mean radius; disc plane at the mean |z| on the side the track
    // heads to), or -1 if it does not within half a turn. crossed is set when that point
    // lies within the surface's module extent (barrel |z|, disc radius) widened by
    // `margin` [cm]: the surfaces the particle traverses, for path ordering and material.
    double pathToSurface(const HelixParams& helix, unsigned surface, double margin, bool& crossed) const;

    // crossLayer for a TBPX surface code (1-4), the disc equivalent for 11-14.
    Crossing crossSurface(const HelixParams& helix, unsigned surface, const MagneticField& field,
                          double predAngleMaxAbs = 0.) const;

    // Propagate the helix to TBPX layer `layer1based` (1..nLayers) and find the
    // module it lands on. Returns Crossing{valid=false} if the helix does not
    // reach the layer radius or lands off every module. `field` supplies the
    // local B at the module center.
    //
    // predAngleMaxAbs (spec §6b, secondary hygiene clamp): if the predicted
    // module-local |cotAlpha| or |cotBeta| exceeds this bound, the crossing is
    // treated as INVALID (valid=false) — a near-grazing predicted crossing whose
    // angle is non-physical. Default <=0 DISABLES the clamp: the analyzer
    // measurement path must stay unclipped so payload fits see the full angle
    // distribution; only the seedlessRefit producer passes a finite bound.
    Crossing crossLayer(const HelixParams& helix, int layer1based, const MagneticField& field,
                        double predAngleMaxAbs = 0.) const;

    // Propagate the helix to disc `disk` (1..nDiscs) on the side it heads to and land it
    // on a module, as crossLayer does for a layer.
    Crossing crossDisc(const HelixParams& helix, int disk, const MagneticField& field,
                       double predAngleMaxAbs = 0.) const;

  private:
    // Newton re-propagation of the helix onto each candidate module plane (starting at
    // transverse arc length sTransverse), containment test, nearest accepted centre.
    Crossing land(const HelixParams& helix,
                  double sTransverse,
                  const std::vector<std::pair<double, const GeomDet*>>& candidates,
                  int surface,
                  const MagneticField& field,
                  double predAngleMaxAbs) const;

    bool built_ = false;
    const TrackerGeometry* geomTag_ = nullptr;  // identity of the cached geometry
    int nLayersReq_ = 0;
    int nDiscsReq_ = 0;
    std::vector<std::vector<const GeomDet*>> layerModules_;  // [layer-1] -> units (GeomDetUnit)
    std::vector<double> layerRadius_;                            // [layer-1] mean r [cm]
    std::vector<double> layerHalfZ_;                             // [layer-1] max |z| of the modules [cm]
    std::vector<std::array<std::vector<const GeomDet*>, 2>> discModules_;  // [disk-1][z < 0, z > 0]
    std::vector<double> discZ_;                                  // [disk-1] mean |z| [cm]
    std::vector<std::pair<double, double>> discR_;               // [disk-1] module radial extent [cm]
    std::vector<double> discMaxHalfDiag_;                        // [disk-1] largest module half-diagonal [cm]
  };

}  // namespace smartpixels

#endif
