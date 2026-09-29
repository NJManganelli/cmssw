#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsOTClusterTruth_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsOTClusterTruth_h

#include <map>

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimDataFormats/Associations/interface/TTClusterAssociationMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

namespace smartpixels {

  // TRUTH-ONLY owner of an OT (TTCluster) cluster: its dominant contributor in the
  // TTClusterAssociationMap, i.e. the TrackingParticle with the largest share of the map's own
  // weight (TP pT summed over the cluster's strip x simlink entries, the weight isGenuine() cuts
  // on). A genuine cluster therefore yields exactly findTrackingParticlePtr(); a combinatoric
  // (merged) cluster, for which findTrackingParticlePtr() is null, still yields its dominant TP.
  // Null for an unknown cluster. Ties go to the lowest TP key.
  template <typename T>
  TrackingParticlePtr dominantTrackingParticle(const TTClusterAssociationMap<T>& assoc,
                                               const TTClusterRefT<T>& cluster) {
    std::map<TrackingParticlePtr, double> weight;
    for (const auto& tp : assoc.findTrackingParticlePtrs(cluster))
      if (tp.isNonnull())
        weight[tp] += tp->pt();
    TrackingParticlePtr best;
    double bestWeight = 0.;
    for (const auto& [tp, w] : weight)
      if (w > bestWeight) {
        best = tp;
        bestWeight = w;
      }
    return best;
  }

}  // namespace smartpixels

#endif
