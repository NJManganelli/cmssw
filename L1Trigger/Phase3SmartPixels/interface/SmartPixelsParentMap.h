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
  // (EncodedEventId.rawId, SimTrack trackId) -> index into the TrackingParticle
  // collection. Same key as ParentMomentumMap, so a digi simlink resolves to a TP
  // IDENTITY, not just a momentum.
  using ParentTpIndexMap = std::map<std::pair<unsigned int, unsigned int>, int>;

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

  // A TrackingParticle owns SEVERAL g4Tracks, so "does this cluster belong to that
  // track's TP" cannot be answered by comparing SimTrack ids -- every g4Track of the
  // TP must map to the same answer. Mapping every g4Track onto the TP's INDEX gives
  // exactly that, and turns the question into an integer comparison.
  inline ParentTpIndexMap buildParentTpIndexMap(const std::vector<TrackingParticle>& tps) {
    ParentTpIndexMap m;
    for (size_t i = 0; i < tps.size(); ++i) {
      const unsigned int evt = tps[i].eventId().rawId();
      for (const auto& g4 : tps[i].g4Tracks())
        m.emplace(std::make_pair(evt, g4.trackId()), static_cast<int>(i));
    }
    return m;
  }

}  // namespace smartpixels

#endif
