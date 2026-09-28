#ifndef DPGAnalysis_Phase2L1TNanoAOD_L1TPRowLink_h
#define DPGAnalysis_Phase2L1TNanoAOD_L1TPRowLink_h

// tpIdx link of a TrackingParticle into the L1TTP nano table, from the key -> row
// ValueMap that L1TrackingParticleTableProducer puts next to the table (-1 = not selected):
//   >= 0 row in L1TTP; -1 not truth-matched; -2 TrackingParticle product or row map
//   unavailable; <= -3 matched TP not selected into L1TTP, -3 - tpIdx = its index in the
//   unwritten TrackingParticle collection. The L1TTP selection is stated in its table doc.

#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cstdint>

namespace l1tnano {

  inline const char* tpIdxDoc() {
    return "row in L1TTP of the matched TrackingParticle (>= 0); -1 not truth-matched; -2 TrackingParticle "
           "product or L1TTP unavailable; <= -3 matched TP not selected into L1TTP (its selection is stated in "
           "the L1TTP table doc), -3 - tpIdx is its index in the TrackingParticle collection "
           "(mix:MergedTrackTruth), which is NOT written -- unlike L1TTP_genPartIdx, where -3 - x is a row of "
           "the written GenPart table. Two objects share a TP iff their tpIdx are equal";
  }

  inline int32_t tpIdxOfKey(const edm::ValueMap<int>* rowMap, edm::ProductID id, int64_t key) {
    if (key < 0)
      return -1;
    if (rowMap == nullptr)
      return -2;
    const int row = rowMap->get(id, key);
    return row >= 0 ? row : -3 - static_cast<int32_t>(key);
  }

  // tpRef: edm::Ref or edm::Ptr into the TrackingParticle collection the map was built on
  // (a different product throws in the ValueMap lookup).
  template <typename R>
  int32_t tpIdx(const edm::ValueMap<int>* rowMap, const R& tpRef) {
    return tpRef.isNull() ? -1 : tpIdxOfKey(rowMap, tpRef.id(), tpRef.key());
  }

  // For callers holding only a key: the map must cover exactly one TrackingParticle product.
  inline edm::ProductID tpProductOf(const edm::ValueMap<int>& rowMap) {
    if (rowMap.idSize() != 1)
      throw cms::Exception("L1TPRowLink") << "TrackingParticle row map covers " << rowMap.idSize()
                                          << " products, expected 1";
    return rowMap.begin().id();
  }

}  // namespace l1tnano

#endif
