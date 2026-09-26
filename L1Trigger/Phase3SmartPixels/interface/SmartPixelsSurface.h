#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsSurface_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsSurface_h

// IT surface codes, disc scope, and the cluster-table row order, in ONE place.
//
// SURFACE CODE: 1-4 = TBPX layers IL1-IL4; 10 + pxfDisk on the discs, so ID1-ID4
// (TFPX disks 1-4, both sides) = 11-14, TFPX 5-8 = 15-18, TEPX = 19-22. These are
// IT-only codes: the OT stub tables' layer column also uses 10 + disk (11-15 = OT
// disks) and analysis code maps OL1-OL6 to 11-16, so an IT surface must never be
// mixed into an OT or combined IT+OT code space without re-mapping.
//
// ROW ORDER of L1TSmartPixelsCluster: every barrel cluster first (DetSets in
// collection order, layers 1..maxLayer), then every in-scope disc cluster. The
// barrel prefix is therefore the historical barrel-only table, row for row.
// Both the table producer and the refit's selClusterIdx use clusterTableOrder, so
// the two cannot diverge by each re-deriving the order.

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smartpixels {

  constexpr unsigned kDiscSurfaceOffset = 10;

  inline bool isDiscSurface(unsigned surface) { return surface > kDiscSurfaceOffset; }

  // Highest pxfDisk kept: none = no discs, sp = ID1-ID4 (the smart-pixel discs),
  // tfpx = TFPX 1-8, all = TFPX + TEPX.
  inline unsigned discScopeMaxDisk(const std::string& scope) {
    if (scope == "none")
      return 0;
    if (scope == "sp")
      return 4;
    if (scope == "tfpx")
      return 8;
    if (scope == "all")
      return 99;
    throw cms::Exception("Configuration") << "discScope must be none|sp|tfpx|all, got '" << scope << "'";
  }

  // Surface code of a pixel module, or 0 when it is outside the configured scope.
  inline unsigned surfaceCode(DetId did, const TrackerTopology& topo, unsigned maxLayer, unsigned maxDisk) {
    if (did.subdetId() == PixelSubdetector::PixelBarrel) {
      const unsigned lay = topo.pxbLayer(did);
      return (lay >= 1 && lay <= maxLayer) ? lay : 0;
    }
    if (did.subdetId() == PixelSubdetector::PixelEndcap) {
      const unsigned disk = topo.pxfDisk(did);
      return (disk >= 1 && disk <= maxDisk) ? kDiscSurfaceOffset + disk : 0;
    }
    return 0;
  }

  // (DetSet, surface code) in cluster-table row order: barrel pass, then disc pass.
  // Held BY VALUE: an edmNew::DetSetVector iterator yields a lightweight DetSet proxy
  // by value, so its address would dangle.
  template <typename DSV>
  std::vector<std::pair<typename DSV::value_type, unsigned>> clusterTableOrder(
      const DSV& hits, const TrackerTopology& topo, unsigned maxLayer, unsigned maxDisk) {
    std::vector<std::pair<typename DSV::value_type, unsigned>> out;
    for (const bool discPass : {false, true})
      for (const auto& ds : hits) {
        const unsigned s = surfaceCode(DetId(ds.detId()), topo, maxLayer, maxDisk);
        if (s != 0 && isDiscSurface(s) == discPass)
          out.emplace_back(ds, s);
      }
    return out;
  }

  // detId -> row of that module's first cluster in the cluster table. Every hit of
  // an in-scope module is a row, so base + position-within-module is the row.
  template <typename DSV>
  std::unordered_map<uint32_t, unsigned> clusterRowBase(const DSV& hits,
                                                        const TrackerTopology& topo,
                                                        unsigned maxLayer,
                                                        unsigned maxDisk) {
    std::unordered_map<uint32_t, unsigned> base;
    unsigned running = 0;
    for (const auto& [ds, s] : clusterTableOrder(hits, topo, maxLayer, maxDisk)) {
      base[ds.detId()] = running;
      running += ds.size();
    }
    return base;
  }

}  // namespace smartpixels

#endif
