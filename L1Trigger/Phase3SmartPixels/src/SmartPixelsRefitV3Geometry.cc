#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitV3.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

namespace smartpixels::v3 {

  void buildGeometry(Geometry& g, const TrackerGeometry& geom, const TrackerTopology& topo) {
    g = Geometry{};
    for (const auto* det : geom.detUnits()) {
      const DetId id = det->geographicalId();
      if (id.det() != DetId::Tracker)
        continue;
      int S = 0;
      if (id.subdetId() == PixelSubdetector::PixelBarrel)
        S = static_cast<int>(topo.pxbLayer(id));
      else if (id.subdetId() == PixelSubdetector::PixelEndcap)
        S = 10 + static_cast<int>(topo.pxfDisk(id));
      if (!((S >= 1 && S <= 4) || (S >= 11 && S <= 14)))
        continue;
      Module m;
      m.detId = id.rawId();
      m.surface = S;
      const GlobalPoint o = det->toGlobal(LocalPoint(0., 0., 0.));
      const GlobalVector ex = det->toGlobal(LocalVector(1., 0., 0.));
      const GlobalVector ey = det->toGlobal(LocalVector(0., 1., 0.));
      const GlobalVector ez = det->toGlobal(LocalVector(0., 0., 1.));
      m.o = {{o.x(), o.y(), o.z()}};
      m.ex = {{ex.x(), ex.y(), ex.z()}};
      m.ey = {{ey.x(), ey.y(), ey.z()}};
      m.ez = {{ez.x(), ez.y(), ez.z()}};
      m.hw = 0.5 * det->surface().bounds().width();
      m.hl = 0.5 * det->surface().bounds().length();
      g.add(m);
    }
    g.finalize();
  }

}  // namespace smartpixels::v3
