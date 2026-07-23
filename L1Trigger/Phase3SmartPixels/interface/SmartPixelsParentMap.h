#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsParentMap_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsParentMap_h

// Parent-momentum lookup for per-digi angle synthesis, keyed by
// (EncodedEventId.rawId, SimTrack trackId) — the same pair a PixelDigiSimLink
// carries. Built PRIMARILY from TrackingParticles' embedded g4Tracks, which
// cover BOTH the signal event and in-time pileup (the persisted g4SimHits
// SimTrack container is signal-only, so a SimTrack-container-keyed map silently
// loses every PU parent — at PU200 that is most of the in-window digis).
// The signal SimTrack container is kept as a fallback for signal parents whose
// TrackingParticles were pruned. Misses degrade gracefully to position-only
// hits at the call site.

#include <map>
#include <utility>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

namespace smartpixels {

  using ParentMomentumMap = std::map<std::pair<unsigned int, unsigned int>, math::XYZTLorentzVectorD>;

  inline ParentMomentumMap buildParentMomentumMap(const std::vector<TrackingParticle>& tps,
                                                  const edm::SimTrackContainer* signalSimTracks) {
    ParentMomentumMap m;
    for (const auto& tp : tps) {
      const unsigned int evt = tp.eventId().rawId();
      for (const auto& g4 : tp.g4Tracks())
        m.emplace(std::make_pair(evt, g4.trackId()), g4.momentum());
    }
    if (signalSimTracks) {
      const unsigned int sig = EncodedEventId(0, 0).rawId();
      for (const auto& st : *signalSimTracks)
        m.emplace(std::make_pair(sig, st.trackId()), st.momentum());  // no overwrite: TPs win
    }
    return m;
  }

}  // namespace smartpixels

#endif
