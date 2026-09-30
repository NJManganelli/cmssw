#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsHelixProjector.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsSurface.h"

#include <algorithm>
#include <cmath>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "MagneticField/Engine/interface/MagneticField.h"

namespace smartpixels {

  namespace {
    // module corners and edge midpoints (global), for the surface extents
    std::array<GlobalPoint, 8> outline(const GeomDet& det) {
      const double hx = 0.5 * det.surface().bounds().width(), hy = 0.5 * det.surface().bounds().length();
      std::array<GlobalPoint, 8> pts;
      unsigned n = 0;
      for (const double fx : {-1., 0., 1.})
        for (const double fy : {-1., 0., 1.})
          if (fx != 0. || fy != 0.)
            pts[n++] = det.toGlobal(LocalPoint(fx * hx, fy * hy, 0.));
      return pts;
    }

    // Past half a turn the helix heads back inward: not a crossing of an outer surface.
    bool beyondHalfTurn(double rInv, double s) { return std::abs(rInv) * s > M_PI; }

    // candidate modules of a disc whose centre is within the largest half-diagonal plus this
    // slack of the mean-plane crossing: covers the +-0.9 cm z spread of the disc faces
    constexpr double kDiscCandidateSlack = 1.5;  // [cm]
  }  // namespace

  void HelixProjector::build(const TrackerGeometry& geom, const TrackerTopology& topo, int nLayers, int nDiscs) {
    const int nreq = std::max(1, std::min(4, nLayers));
    const int ndreq = std::max(0, std::min(4, nDiscs));
    if (built_ && geomTag_ == &geom && nLayersReq_ == nreq && nDiscsReq_ == ndreq)
      return;  // cache still valid

    built_ = false;
    geomTag_ = &geom;
    nLayersReq_ = nreq;
    nDiscsReq_ = ndreq;
    layerModules_.assign(nreq, {});
    layerRadius_.assign(nreq, 0.);
    layerHalfZ_.assign(nreq, 0.);
    discModules_.assign(ndreq, {});
    discZ_.assign(ndreq, 0.);
    discR_.assign(ndreq, {1e9, 0.});
    discMaxHalfDiag_.assign(ndreq, 0.);
    std::vector<int> nPerLayer(nreq, 0);
    std::vector<int> nPerDisc(ndreq, 0);

    for (const auto* det : geom.detUnits()) {
      const DetId id = det->geographicalId();
      if (id.det() != DetId::Tracker)
        continue;
      if (id.subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap)) {
        const unsigned int disk = topo.pxfDisk(id);
        if (disk < 1 || static_cast<int>(disk) > ndreq)
          continue;
        const int d = disk - 1;
        discModules_[d][det->position().z() > 0. ? 1 : 0].push_back(det);
        discZ_[d] += std::abs(det->position().z());
        nPerDisc[d] += 1;
        for (const auto& p : outline(*det)) {
          discR_[d].first = std::min<double>(discR_[d].first, p.perp());
          discR_[d].second = std::max<double>(discR_[d].second, p.perp());
        }
        discMaxHalfDiag_[d] = std::max<double>(
            discMaxHalfDiag_[d], 0.5 * std::hypot(det->surface().bounds().width(), det->surface().bounds().length()));
        continue;
      }
      if (id.subdetId() != static_cast<int>(PixelSubdetector::PixelBarrel))
        continue;
      const unsigned int lay = topo.pxbLayer(id);
      if (lay < 1 || static_cast<int>(lay) > nreq)
        continue;
      layerModules_[lay - 1].push_back(det);
      layerRadius_[lay - 1] += det->position().perp();
      nPerLayer[lay - 1] += 1;
      for (const auto& p : outline(*det))
        layerHalfZ_[lay - 1] = std::max<double>(layerHalfZ_[lay - 1], std::abs(p.z()));
    }
    for (int l = 0; l < nreq; ++l) {
      if (nPerLayer[l] > 0)
        layerRadius_[l] /= nPerLayer[l];
    }
    for (int d = 0; d < ndreq; ++d) {
      if (nPerDisc[d] > 0)
        discZ_[d] /= nPerDisc[d];
    }
    built_ = true;
  }

  double HelixProjector::pathToSurface(const HelixParams& h, unsigned surface, double margin, bool& crossed) const {
    crossed = false;
    const double cosPhi0 = std::cos(h.phi0);
    const double sinPhi0 = std::sin(h.phi0);
    double s = -1.;
    if (!isDiscSurface(surface)) {
      const int l = static_cast<int>(surface) - 1;
      if (!built_ || l < 0 || l >= static_cast<int>(layerModules_.size()) || layerRadius_[l] <= 0.)
        return -1.;
      const double Rl = layerRadius_[l];
      // the cylinder intersection of crossLayer
      if (std::abs(h.rInv) < 1e-6) {
        const double b = h.x0 * cosPhi0 + h.y0 * sinPhi0;
        const double disc = b * b - (h.x0 * h.x0 + h.y0 * h.y0 - Rl * Rl);
        if (disc < 0)
          return -1.;
        s = -b + std::sqrt(disc);
        if (s < 0)
          return -1.;
      } else {
        const double R = 1.0 / h.rInv;
        const double dcenter = std::hypot(h.x0 + R * sinPhi0, h.y0 - R * cosPhi0);
        const double absR = std::abs(R);
        if (Rl > dcenter + absR || Rl < std::abs(dcenter - absR))
          return -1.;
        const double cosArg = (absR * absR + dcenter * dcenter - Rl * Rl) / (2.0 * absR * dcenter);
        if (cosArg < -1.0 || cosArg > 1.0)
          return -1.;
        s = std::acos(cosArg) * absR;
      }
      crossed = std::abs(h.z0 + h.tanL * s) <= layerHalfZ_[l] + margin;
      return s;
    }
    const int d = static_cast<int>(surface - kDiscSurfaceOffset) - 1;
    if (!built_ || d < 0 || d >= static_cast<int>(discZ_.size()) || std::abs(h.tanL) < 1e-9)
      return -1.;
    s = ((h.tanL > 0. ? 1. : -1.) * discZ_[d] - h.z0) / h.tanL;
    if (!(s > 0.) || beyondHalfTurn(h.rInv, s))
      return -1.;
    double gx, gy;
    if (std::abs(h.rInv) < 1e-6) {
      gx = h.x0 + s * cosPhi0;
      gy = h.y0 + s * sinPhi0;
    } else {
      gx = h.x0 - (std::sin(h.phi0 - h.rInv * s) - sinPhi0) / h.rInv;
      gy = h.y0 + (std::cos(h.phi0 - h.rInv * s) - cosPhi0) / h.rInv;
    }
    const double r = std::hypot(gx, gy);
    crossed = r >= discR_[d].first - margin && r <= discR_[d].second + margin;
    return s;
  }

  Crossing HelixProjector::crossSurface(const HelixParams& h, unsigned surface, const MagneticField& field,
                                        double predAngleMaxAbs) const {
    return isDiscSurface(surface)
               ? crossDisc(h, static_cast<int>(surface - kDiscSurfaceOffset), field, predAngleMaxAbs)
               : crossLayer(h, static_cast<int>(surface), field, predAngleMaxAbs);
  }

  Crossing HelixProjector::crossLayer(const HelixParams& h, int layer1based, const MagneticField& field,
                                      double predAngleMaxAbs) const {
    Crossing out;
    if (!built_ || layer1based < 1 || layer1based > static_cast<int>(layerModules_.size()))
      return out;
    const int l = layer1based - 1;
    const double Rl = layerRadius_[l];
    if (Rl <= 0.)
      return out;

    // TTTrack CONVENTIONS, used as they are (no sign flips anywhere below):
    // rInv > 0 is a positive charge, which turns CLOCKWISE in +Bz, so the
    // direction angle along the transverse arc s is phi(s) = phi0 - rInv * s; the
    // POCA is (d0 sin phi0, -d0 cos phi0); the circle centre sits at
    // (x0 + sin(phi0)/rInv, y0 - cos(phi0)/rInv). An earlier version wrote the
    // helix in the counter-clockwise maths convention and negated rInv on entry
    // ("SIGN FIX", validated on RelVal PU200: seed helix lands on its own OT
    // stubs, per-hit pulls flat L1->L4, matched-TP helix hits the OT stubs); this
    // form is the same helix, bit for bit.
    const double rInv = h.rInv;
    const double phi0 = h.phi0;
    const double x0 = h.x0;
    const double y0 = h.y0;
    const double z0 = h.z0;
    const double tanL = h.tanL;
    const double cosPhi0 = std::cos(phi0);
    const double sinPhi0 = std::sin(phi0);

    // Analytic helix -> cylinder(R=Rl) intersection in the transverse plane.
    double sTransverse;  // transverse path length to the crossing
    if (std::abs(rInv) < 1e-6) {
      // straight line: |(x0,y0) + s*dir| = Rl
      const double b = x0 * cosPhi0 + y0 * sinPhi0;
      const double c = x0 * x0 + y0 * y0 - Rl * Rl;
      const double disc = b * b - c;
      if (disc < 0)
        return out;
      sTransverse = -b + std::sqrt(disc);
      if (sTransverse < 0)
        return out;
    } else {
      // Circle center is offset by 1/rInv perpendicular to the momentum, on the
      // right of the direction for a positive (clockwise) track.
      const double R = 1.0 / rInv;  // signed
      const double cx = x0 + R * sinPhi0;
      const double cy = y0 - R * cosPhi0;
      const double dcenter = std::hypot(cx, cy);
      const double absR = std::abs(R);
      if (Rl > dcenter + absR || Rl < std::abs(dcenter - absR))
        return out;
      // Law of cosines for the turning angle from POCA to the crossing.
      // NB: must use |R| — a signed R flips cosArg for negative charge and
      // sends the crossing to the far side of the circle (psi -> pi - psi).
      const double cosArg = (absR * absR + dcenter * dcenter - Rl * Rl) / (2.0 * absR * dcenter);
      if (cosArg < -1.0 || cosArg > 1.0)
        return out;
      const double delta = std::acos(std::max(-1.0, std::min(1.0, cosArg)));
      sTransverse = delta * absR;  // arc length; the turning SENSE lives in rInv's sign
    }

    // Global crossing point via helix stepping in transverse arc length.
    double gx, gy;
    if (std::abs(rInv) < 1e-6) {
      gx = x0 + sTransverse * cosPhi0;
      gy = y0 + sTransverse * sinPhi0;
    } else {
      const double psi = -rInv * sTransverse;  // turning angle (clockwise for rInv > 0)
      gx = x0 - (std::sin(phi0 + psi) - sinPhi0) / rInv;
      gy = y0 + (std::cos(phi0 + psi) - cosPhi0) / rInv;
    }
    const double gz = z0 + tanL * sTransverse;
    const GlobalPoint gp(gx, gy, gz);

    // ---- v1 module lookup: k-nearest candidates + plane re-propagation ----
    // The cylinder crossing above lives at the MEAN layer radius; staggered
    // TBPX planks sit at r = Rl +- dr, so containment-testing the off-plane
    // point loses modules (v0 symptom: L3/L4 central efficiency 0.15/0.045).
    // Instead: pick the k nearest module centers, RE-PROPAGATE the helix onto
    // each candidate's actual plane (Newton in transverse arc length), and
    // containment-test the on-plane intersection.
    constexpr int kCandidates = 4;
    std::vector<std::pair<double, const GeomDetUnit*>> cand;
    cand.reserve(layerModules_[l].size());
    for (const auto* det : layerModules_[l])
      cand.emplace_back((det->position() - gp).mag2(), det);
    const int nCand = std::min<int>(kCandidates, cand.size());
    std::partial_sort(cand.begin(), cand.begin() + nCand, cand.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
    cand.resize(nCand);
    return land(h, sTransverse, cand, layer1based, field, predAngleMaxAbs);
  }

  Crossing HelixProjector::crossDisc(const HelixParams& h, int disk, const MagneticField& field,
                                     double predAngleMaxAbs) const {
    if (!built_ || disk < 1 || disk > static_cast<int>(discZ_.size()) || std::abs(h.tanL) < 1e-9)
      return Crossing{};
    const int d = disk - 1;
    const int side = h.tanL > 0. ? 1 : 0;
    // the mean disc plane on the side the track heads to; the modules sit on four faces
    // within +-0.9 cm of it, which the plane re-propagation in land() resolves
    const double sTransverse = ((side ? 1. : -1.) * discZ_[d] - h.z0) / h.tanL;
    if (!(sTransverse > 0.) || beyondHalfTurn(h.rInv, sTransverse))
      return Crossing{};
    double gx, gy;
    if (std::abs(h.rInv) < 1e-6) {
      gx = h.x0 + sTransverse * std::cos(h.phi0);
      gy = h.y0 + sTransverse * std::sin(h.phi0);
    } else {
      gx = h.x0 - (std::sin(h.phi0 - h.rInv * sTransverse) - std::sin(h.phi0)) / h.rInv;
      gy = h.y0 + (std::cos(h.phi0 - h.rInv * sTransverse) - std::cos(h.phi0)) / h.rInv;
    }
    const GlobalPoint gp(gx, gy, h.z0 + h.tanL * sTransverse);
    // every module that can contain the crossing (overlapping faces, neighbouring rings),
    // nearest first, instead of the barrel's fixed k nearest
    const double maxDist2 = std::pow(discMaxHalfDiag_[d] + kDiscCandidateSlack, 2);
    std::vector<std::pair<double, const GeomDetUnit*>> cand;
    for (const auto* det : discModules_[d][side]) {
      const double d2 = (det->position() - gp).mag2();
      if (d2 < maxDist2)
        cand.emplace_back(d2, det);
    }
    std::sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    return land(h, sTransverse, cand, static_cast<int>(kDiscSurfaceOffset) + disk, field, predAngleMaxAbs);
  }

  Crossing HelixProjector::land(const HelixParams& h,
                                double sTransverse,
                                const std::vector<std::pair<double, const GeomDet*>>& cand,
                                int surface,
                                const MagneticField& field,
                                double predAngleMaxAbs) const {
    Crossing out;
    const double rInv = h.rInv;
    const double phi0 = h.phi0;
    const double x0 = h.x0;
    const double y0 = h.y0;
    const double z0 = h.z0;
    const double tanL = h.tanL;
    const double cosPhi0 = std::cos(phi0);
    const double sinPhi0 = std::sin(phi0);
    const int nCand = cand.size();

    const auto helixAt = [&](double s, double& px, double& py, double& pz,
                             double& dx, double& dy, double& dz) {
      if (std::abs(rInv) < 1e-6) {
        px = x0 + s * cosPhi0;
        py = y0 + s * sinPhi0;
        dx = cosPhi0;
        dy = sinPhi0;
      } else {
        const double ps = -rInv * s;
        px = x0 - (std::sin(phi0 + ps) - sinPhi0) / rInv;
        py = y0 + (std::cos(phi0 + ps) - cosPhi0) / rInv;
        dx = std::cos(phi0 + ps);
        dy = std::sin(phi0 + ps);
      }
      pz = z0 + tanL * s;
      dz = tanL;
    };

    const GeomDetUnit* best = nullptr;
    double bestDist2 = 1e18;
    double bestS = sTransverse;
    LocalPoint bestLocal;
    GlobalPoint bestGlobal;
    for (int ic = 0; ic < nCand; ++ic) {
      const GeomDetUnit* det = cand[ic].second;
      const GlobalPoint p0 = det->position();
      const auto nvec = det->surface().normalVector();  // global unit normal
      // Newton: solve (helix(s) - p0) . n = 0 starting from the cylinder crossing.
      double s = sTransverse;
      bool converged = false;
      for (int it = 0; it < 5; ++it) {
        double px, py, pz, dx, dy, dz;
        helixAt(s, px, py, pz, dx, dy, dz);
        const double f = (px - p0.x()) * nvec.x() + (py - p0.y()) * nvec.y() + (pz - p0.z()) * nvec.z();
        const double fp = dx * nvec.x() + dy * nvec.y() + dz * nvec.z();
        if (std::abs(f) < 1e-4) {  // on-plane to 1 um
          converged = true;
          break;
        }
        if (std::abs(fp) < 1e-9)
          break;  // grazing incidence: give up on this candidate
        s -= f / fp;
        if (!(s > 0.) || s > 3.0 * sTransverse + 10.0)
          break;  // unphysical step
      }
      if (!converged)
        continue;
      double px, py, pz, dx, dy, dz;
      helixAt(s, px, py, pz, dx, dy, dz);
      const GlobalPoint gpl(px, py, pz);
      const LocalPoint lp = det->toLocal(gpl);
      if (!det->surface().bounds().inside(LocalPoint(lp.x(), lp.y(), 0.)))
        continue;
      const double d2 = (det->position() - gpl).mag2();
      if (d2 < bestDist2) {
        bestDist2 = d2;
        best = det;
        bestS = s;
        bestLocal = lp;
        bestGlobal = gpl;
      }
    }
    if (!best)
      return out;

    // Momentum direction at the ACCEPTED (on-plane) crossing point.
    const double psi = -rInv * bestS;
    const double pt = h.pt;
    const GlobalVector gmom(pt * std::cos(phi0 + psi), pt * std::sin(phi0 + psi), pt * tanL);

    const PixelGeomDetUnit* pixDet = dynamic_cast<const PixelGeomDetUnit*>(best);
    if (!pixDet)
      return out;

    const LocalVector lmom = best->toLocal(gmom);
    const double lpz = (std::abs(lmom.z()) > 1e-9) ? lmom.z() : 1e-9;

    const GlobalVector gB = field.inTesla(best->position());
    const LocalVector lB = best->toLocal(gB);

    out.valid = true;
    out.layer = surface;
    out.det = pixDet;
    out.detId = best->geographicalId().rawId();
    out.global = bestGlobal;
    out.local = bestLocal;
    out.globalMom = gmom;
    out.localMom = lmom;
    out.cotAlpha = lmom.x() / lpz;
    out.cotBeta = lmom.y() / lpz;
    // Predicted-angle grazing clamp (spec §6b, secondary hygiene): a
    // near-grazing predicted crossing yields a non-physical local angle. When a
    // finite bound is supplied (producer only; analyzer passes 0 = disabled),
    // reject the crossing. Investigation: ~18/68k crossings on the PU sample have
    // |cotAlpha| up to 51.8 with windowMult~0; these are not the gate drivers but
    // are physically meaningless. Kept out of the analyzer so payload fits see the
    // full distribution.
    if (predAngleMaxAbs > 0. &&
        (std::abs(out.cotAlpha) > predAngleMaxAbs || std::abs(out.cotBeta) > predAngleMaxAbs))
      return Crossing{};  // invalid
    out.bLocalX = lB.x();
    out.bLocalY = lB.y();
    out.bLocalZ = lB.z();
    return out;
  }

}  // namespace smartpixels
