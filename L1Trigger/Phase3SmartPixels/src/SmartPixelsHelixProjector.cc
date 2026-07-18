#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsHelixProjector.h"

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

  void HelixProjector::build(const TrackerGeometry& geom, const TrackerTopology& topo, int nLayers) {
    const int nreq = std::max(1, std::min(4, nLayers));
    if (built_ && geomTag_ == &geom && nLayersReq_ == nreq)
      return;  // cache still valid

    built_ = false;
    geomTag_ = &geom;
    nLayersReq_ = nreq;
    layerModules_.assign(nreq, {});
    layerRadius_.assign(nreq, 0.);
    std::vector<int> nPerLayer(nreq, 0);

    for (const auto* det : geom.detUnits()) {
      const DetId id = det->geographicalId();
      if (id.det() != DetId::Tracker)
        continue;
      if (id.subdetId() != static_cast<int>(PixelSubdetector::PixelBarrel))
        continue;
      const unsigned int lay = topo.pxbLayer(id);
      if (lay < 1 || static_cast<int>(lay) > nreq)
        continue;
      layerModules_[lay - 1].push_back(det);
      layerRadius_[lay - 1] += det->position().perp();
      nPerLayer[lay - 1] += 1;
    }
    for (int l = 0; l < nreq; ++l) {
      if (nPerLayer[l] > 0)
        layerRadius_[l] /= nPerLayer[l];
    }
    built_ = true;
  }

  Crossing HelixProjector::crossLayer(const HelixParams& h, int layer1based, const MagneticField& field) const {
    Crossing out;
    if (!built_ || layer1based < 1 || layer1based > static_cast<int>(layerModules_.size()))
      return out;
    const int l = layer1based - 1;
    const double Rl = layerRadius_[l];
    if (Rl <= 0.)
      return out;

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
      // Circle center is offset by 1/rInv perpendicular to the momentum.
      const double R = 1.0 / rInv;  // signed
      const double cx = x0 - R * sinPhi0;
      const double cy = y0 + R * cosPhi0;
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
      const double dpsi = (rInv > 0) ? delta : -delta;  // turning direction from charge sign
      sTransverse = std::abs(dpsi) * absR;
    }

    // Global crossing point via helix stepping in transverse arc length.
    double gx, gy;
    if (std::abs(rInv) < 1e-6) {
      gx = x0 + sTransverse * cosPhi0;
      gy = y0 + sTransverse * sinPhi0;
    } else {
      const double psi = rInv * sTransverse;  // turning angle
      gx = x0 + (std::sin(phi0 + psi) - sinPhi0) / rInv;
      gy = y0 - (std::cos(phi0 + psi) - cosPhi0) / rInv;
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

    const auto helixAt = [&](double s, double& px, double& py, double& pz,
                             double& dx, double& dy, double& dz) {
      if (std::abs(rInv) < 1e-6) {
        px = x0 + s * cosPhi0;
        py = y0 + s * sinPhi0;
        dx = cosPhi0;
        dy = sinPhi0;
      } else {
        const double ps = rInv * s;
        px = x0 + (std::sin(phi0 + ps) - sinPhi0) / rInv;
        py = y0 - (std::cos(phi0 + ps) - cosPhi0) / rInv;
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
    const double psi = rInv * bestS;
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
    out.layer = layer1based;
    out.det = pixDet;
    out.detId = best->geographicalId().rawId();
    out.global = bestGlobal;
    out.local = bestLocal;
    out.globalMom = gmom;
    out.localMom = lmom;
    out.cotAlpha = lmom.x() / lpz;
    out.cotBeta = lmom.y() / lpz;
    out.bLocalX = lB.x();
    out.bLocalY = lB.y();
    out.bLocalZ = lB.z();
    return out;
  }

}  // namespace smartpixels
