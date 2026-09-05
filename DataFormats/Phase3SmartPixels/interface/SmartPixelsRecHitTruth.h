#ifndef DataFormats_Phase3SmartPixels_SmartPixelsRecHitTruth_h
#define DataFormats_Phase3SmartPixels_SmartPixelsRecHitTruth_h

// Truth association for SmartPixelsRecHit, as a SEPARATE product.
//
// CMSSW convention: reco objects carry no truth, and a distinct product maps
// reco -> simulation (PixelDigiSimLink, ClusterTPAssociation,
// TTCluster/TTStub/TTTrackAssociationMap). SmartPixelsRecHit therefore holds no
// truth at all and this collection holds it, so an inference path can consume
// the former and provably not the latter.
//
// ALIGNMENT INVARIANT: produced by the same module as the rechits, with the same
// DetSet ordering and the same order within each DetSet, so entry j of det d here
// describes entry j of det d there. The producer asserts equal total size and
// equal per-DetSet sizes; `idxInDet` is stored so a consumer can verify
// alignment itself rather than trusting it.
//
// The identity is a TrackingParticleRef, not a SimTrack id: a TrackingParticle
// owns several g4Tracks, so comparing SimTrack ids answers the wrong question.

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include <cstdint>

class SmartPixelsRecHitTruth {
public:
  SmartPixelsRecHitTruth() = default;

  // Dominant charge contributor: the TrackingParticle with the largest summed
  // pixel ADC over this cluster. Null ref when no pixel of the cluster carries a
  // simlink (noise-like) or the contributor has no TrackingParticle.
  const TrackingParticleRef& dominantTp() const { return tp_; }
  bool hasTp() const { return tp_.isNonnull(); }

  // Share of the cluster charge held by the dominant contributor. < 1 means the
  // cluster is shared, which biases its position and makes its incidence angle
  // ill-defined; `merged` flags a runner-up above the configured fraction.
  float chargeFrac() const { return chargeFrac_; }
  bool merged() const { return merged_; }

  // TRUE incidence angles at this module, from the dominant TP's helix
  // propagated to the hit -- i.e. the angle the sensor is trying to measure,
  // BEFORE the PixelAV response is applied. Sentinel if unavailable.
  float trueCotAlpha() const { return trueCotAlpha_; }
  float trueCotBeta() const { return trueCotBeta_; }

  uint16_t idxInDet() const { return idxInDet_; }

  void set(const TrackingParticleRef& tp, float frac, bool merged,
           float trueCotAlpha, float trueCotBeta, uint16_t idxInDet) {
    tp_ = tp;
    chargeFrac_ = frac;
    merged_ = merged;
    trueCotAlpha_ = trueCotAlpha;
    trueCotBeta_ = trueCotBeta;
    idxInDet_ = idxInDet;
  }

private:
  TrackingParticleRef tp_;
  float chargeFrac_ = -999.f;
  bool merged_ = false;
  float trueCotAlpha_ = -999.f, trueCotBeta_ = -999.f;
  uint16_t idxInDet_ = 0;
};

using SmartPixelsRecHitTruthCollection = edmNew::DetSetVector<SmartPixelsRecHitTruth>;

#endif
