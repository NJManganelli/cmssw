#ifndef DataFormats_Phase3SmartPixels_SmartPixelsRecHit_h
#define DataFormats_Phase3SmartPixels_SmartPixelsRecHit_h

// A SmartPixels measurement: one per IT pixel CLUSTER, carrying what a
// smart-pixel sensor would emit for that cluster.
//
// NAMING follows the CMSSW convention deliberately. A "Cluster" is the raw
// grouping of adjacent fired pixels with position in pixel units and no error; a
// "RecHit" is a calibrated measurement with a local position, a position error,
// and a Ref back to its cluster (see SiPixelRecHit / SiPixelRecHitCollection).
// This object is a RecHit in exactly that sense, extended with the angle
// estimate a smart pixel adds. Calling it a "hit" without that distinction is
// what previously allowed individual PIXELS to be treated as hits in a refit.
//
// SELF-CONTAINED by design: it carries position and error rather than requiring
// a second dereference of the SiPixelRecHit collection. The ClusterRef is kept
// for provenance, not as the primary access path.
//
// TRUTH-FREE by design. The angle is *derived* from simulation truth today, but
// no truth identifier, no TrackingParticle link and no charge-share information
// lives here. Truth associations are a separate product, per the CMSSW
// convention set by PixelDigiSimLink, ClusterTPAssociation and the
// TTCluster/TTStub/TTTrack association maps. Keeping truth out is what stops it
// leaking into logic that claims to be deployable.

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include <cstdint>

class SmartPixelsRecHit {
public:
  using ClusterRef = edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster>;

  SmartPixelsRecHit() = default;

  SmartPixelsRecHit(const LocalPoint& pos,
                    const LocalError& err,
                    uint32_t detId,
                    uint16_t sizeX,
                    uint16_t sizeY,
                    uint16_t size,
                    float charge,
                    const ClusterRef& cluster)
      : position_(pos),
        error_(err),
        detId_(detId),
        sizeX_(sizeX),
        sizeY_(sizeY),
        size_(size),
        charge_(charge),
        cluster_(cluster) {}

  // --- position, from the pixel CPE (offline Phase-2 reco quality; NOT a model
  //     of what an ASIC computes on-chip) ---
  const LocalPoint& localPosition() const { return position_; }
  const LocalError& localPositionError() const { return error_; }
  uint32_t detId() const { return detId_; }

  // --- cluster shape ---
  uint16_t sizeX() const { return sizeX_; }        // bounding-box extent, local x
  uint16_t sizeY() const { return sizeY_; }        // bounding-box extent, local y
  uint16_t size() const { return size_; }          // FIRED-PIXEL COUNT (not the box)
  float charge() const { return charge_; }         // total cluster charge [ADC]
  ClusterRef cluster() const { return cluster_; }

  // --- the SmartPixels angle estimate -------------------------------------
  // Incidence angles in the module-local frame, PixelAV convention
  //   cotAlpha = p_x_local / p_z_local,  cotBeta = p_y_local / p_z_local
  // evaluated AT THE MODULE. sigAlpha/sigBeta are the estimator's resolution
  // from the PixelAV angle-response payload. has* is false when the sensor would
  // not report that angle at all (payload validity gate, or a grazing incidence
  // beyond the physical bound), which is a real operating mode and NOT a filler
  // value to be silently averaged over.
  float cotAlpha() const { return cotAlpha_; }
  float cotBeta() const { return cotBeta_; }
  float sigAlpha() const { return sigAlpha_; }
  float sigBeta() const { return sigBeta_; }
  bool hasAlpha() const { return hasAlpha_; }
  bool hasBeta() const { return hasBeta_; }

  void setAngles(float cotAlpha, float cotBeta, float sigAlpha, float sigBeta, bool hasA, bool hasB) {
    cotAlpha_ = cotAlpha;
    cotBeta_ = cotBeta;
    sigAlpha_ = sigAlpha;
    sigBeta_ = sigBeta;
    hasAlpha_ = hasA;
    hasBeta_ = hasB;
  }

private:
  LocalPoint position_;
  LocalError error_;
  uint32_t detId_ = 0;
  uint16_t sizeX_ = 0, sizeY_ = 0, size_ = 0;
  float charge_ = 0.f;
  ClusterRef cluster_;

  float cotAlpha_ = -999.f, cotBeta_ = -999.f;
  float sigAlpha_ = -999.f, sigBeta_ = -999.f;
  bool hasAlpha_ = false, hasBeta_ = false;
};

using SmartPixelsRecHitCollection = edmNew::DetSetVector<SmartPixelsRecHit>;

#endif
