// Writes the inner-tracker pixel module geometry to JSON, once, so offline
// studies can project a track onto ANY module's plane -- including the
// neighbours of the one a helix nominally crosses -- in that module's own local
// frame, with the same transform the CMSSW producers use.
//
// FRAME CONVENTION, stated so no consumer has to guess one: each module stores
// its origin and its three local unit axes AS GLOBAL VECTORS, obtained from
// GeomDet::toGlobal itself, so for a local point (lx, ly, lz)
//     global = origin + lx*ex + ly*ey + lz*ez
// and a global direction d has local components (d.ex, d.ey, d.ez). That is
// exactly GeomDet::toLocal, which is what both the helix projector
// (Crossing::cotAlpha = lmom.x/lmom.z) and SmartPixelsRecHitProducer
// (cotAlpha = plv.x/plv.z) apply -- so angles built from this file are
// directly comparable to the nano's projected and cluster angles, tilted or
// flipped modules included. No rotation-matrix row/column convention is
// involved. ez is the plane normal. halfWidth is along ex, halfLength along ey.
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

class SmartPixelsModuleGeometryDumper : public edm::one::EDAnalyzer<> {
public:
  explicit SmartPixelsModuleGeometryDumper(const edm::ParameterSet& cfg)
      : geomToken_(esConsumes()),
        topoToken_(esConsumes()),
        outputFile_(cfg.getParameter<std::string>("outputFile")) {}

  void analyze(const edm::Event&, const edm::EventSetup& iSetup) override {
    if (done_)
      return;
    done_ = true;
    const TrackerGeometry& geom = iSetup.getData(geomToken_);
    const TrackerTopology& topo = iSetup.getData(topoToken_);

    std::ofstream out(outputFile_);
    if (!out)
      throw cms::Exception("SmartPixelsModuleGeometryDumper") << "cannot open " << outputFile_;
    out << std::setprecision(10);
    out << "{\n \"convention\": \"global = origin + lx*ex + ly*ey + lz*ez; local components of a global "
           "direction d are (d.ex, d.ey, d.ez) == GeomDet::toLocal. cotAlpha = d.ex/d.ez, cotBeta = d.ey/d.ez. "
           "halfWidth along ex, halfLength along ey, units cm.\",\n \"modules\": [\n";

    auto vec = [&out](const char* key, double x, double y, double z) {
      out << "\"" << key << "\": [" << x << ", " << y << ", " << z << "]";
    };
    unsigned int n = 0;
    double worstOrtho = 0.;
    for (const auto* det : geom.detUnits()) {
      const DetId id = det->geographicalId();
      if (id.det() != DetId::Tracker)
        continue;
      const int sub = id.subdetId();
      if (sub != PixelSubdetector::PixelBarrel && sub != PixelSubdetector::PixelEndcap)
        continue;
      const auto* pix = dynamic_cast<const PixelGeomDetUnit*>(det);
      if (!pix)
        continue;
      const GlobalPoint o = det->toGlobal(LocalPoint(0., 0., 0.));
      const GlobalVector ex = det->toGlobal(LocalVector(1., 0., 0.));
      const GlobalVector ey = det->toGlobal(LocalVector(0., 1., 0.));
      const GlobalVector ez = det->toGlobal(LocalVector(0., 0., 1.));
      // a right-handed orthonormal frame is assumed downstream; check, do not trust
      worstOrtho = std::max({worstOrtho,
                             std::abs(static_cast<double>(ex.dot(ey))),
                             std::abs(static_cast<double>(ex.dot(ez))),
                             std::abs(static_cast<double>(ey.dot(ez))),
                             std::abs(static_cast<double>(ex.cross(ey).dot(ez)) - 1.)});
      const Bounds& b = det->surface().bounds();
      const PixelTopology& t = pix->specificTopology();
      const auto pitch = t.pitch();

      out << (n ? ",\n" : "") << "  {\"detId\": " << id.rawId() << ", ";
      if (sub == PixelSubdetector::PixelBarrel) {
        out << "\"subdet\": \"PXB\", \"layer\": " << topo.pxbLayer(id) << ", \"ladder\": " << topo.pxbLadder(id)
            << ", \"module\": " << topo.pxbModule(id) << ", ";
      } else {
        out << "\"subdet\": \"PXF\", \"side\": " << topo.pxfSide(id) << ", \"disk\": " << topo.pxfDisk(id)
            << ", \"blade\": " << topo.pxfBlade(id) << ", \"module\": " << topo.pxfModule(id) << ", ";
      }
      vec("origin", o.x(), o.y(), o.z());
      out << ", ";
      vec("ex", ex.x(), ex.y(), ex.z());
      out << ", ";
      vec("ey", ey.x(), ey.y(), ey.z());
      out << ", ";
      vec("ez", ez.x(), ez.y(), ez.z());
      out << ", \"halfWidth\": " << 0.5 * b.width() << ", \"halfLength\": " << 0.5 * b.length()
          << ", \"thickness\": " << b.thickness() << ", \"pitchX\": " << pitch.first << ", \"pitchY\": " << pitch.second
          << ", \"nrows\": " << t.nrows() << ", \"ncolumns\": " << t.ncolumns() << "}";
      ++n;
    }
    out << "\n ],\n \"n_modules\": " << n << ",\n \"max_frame_orthonormality_violation\": " << worstOrtho << "\n}\n";
    if (worstOrtho > 1e-6)
      throw cms::Exception("SmartPixelsModuleGeometryDumper")
          << "a module frame is not right-handed orthonormal (worst violation " << worstOrtho << ")";
    edm::LogPrint("SmartPixelsModuleGeometryDumper") << "wrote " << n << " pixel modules to " << outputFile_;
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("outputFile", "spix_module_geometry.json");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const std::string outputFile_;
  bool done_ = false;
};

DEFINE_FWK_MODULE(SmartPixelsModuleGeometryDumper);
