// SmartPixelsRecHitProducer: IT pixel rec hits -> SmartPixels rec hits, with the
// sensor's angle estimate attached ONCE, in ONE place.
//
// WHY THIS EXISTS. The angle synthesis used to live inside
// L1SmartPixelsTrackProducer::produce(), four loops deep, on a function-local
// struct that was discarded every iteration. It was therefore impossible to reuse:
// the nano cluster table could only get an angle by copying those lines, which
// would have made a second copy of a contract that can drift silently. The truth
// angle derivation was ALREADY duplicated that way, between the refit producer
// and SmartPixelsPayloadAnalyzer. This producer is the single source; the refit
// and every table consume its output.
//
// THE ANGLE DEFINITION, AND THE BUG IT FIXES.
//
// The old base angle was the dominant contributor's momentum AT ITS PRODUCTION
// VERTEX, rotated into the module frame. That is not the incidence angle: the
// particle bends between the vertex and the module. For pT = 2 GeV in 3.8 T the
// helix radius is ~175 cm, so the direction turns by ~0.017 rad by TBPX L1
// (r = 3 cm) and ~0.091 rad by L4 (r = 16 cm) -- the same size as, or larger
// than, the sensor resolution the PixelAV payload describes (rms ~0.02 in cot).
// So the old "measured" angle was a correct sensor smear applied to a base angle
// taken at the wrong point on the trajectory, and its disagreement with the
// track's projected angle grew with radius as bending does, not as resolution
// does.
//
// Here the base angle is the dominant TP's helix PROPAGATED TO THE HIT. For a
// pure helix in a solenoidal field the transverse direction at any point is
// perpendicular to the vector from the circle centre to that point, and p_z is
// unchanged. So with the TP's charge, production vertex and momentum we get the
// direction at the measured hit position exactly, with no propagator needed.
//
// WHAT THIS STILL NEGLECTS: multiple scattering between the vertex and the
// module. The helix is the no-scattering limit. PSimHit::localDirection() would
// include it, but PSimHits are g4SimHits/"SIM" -- SIGNAL ONLY -- while digis and
// TrackingParticles are post-mixing, so PSimHit angles exist for only ~1.7% of
// PU200 clusters (451 signal clusters/event against 26 479 total). Using PSimHit
// where available and this elsewhere would give signal a correct angle and
// pileup a worse one, manufacturing discrimination out of a truth-access
// asymmetry. So the helix is used for ALL clusters, uniformly, and PSimHit is
// reserved as a validation reference on the signal subset where both exist.
//
// NOISE CLUSTERS (no simlink on any pixel) get an angle from the noiseSet
// inverse CDF, if one is configured. This matters more than it looks: while
// unlinked clusters carried NO angle, "reports an angle" was a PERFECT truth
// proxy -- measured hasAlpha 99.5% for TrackingParticle-linked clusters against
// 0.0% for unlinked ones -- so every angle-weighted selection rule and every
// angle-using MVA was buying "is this cluster real" for free. The chi2 weight
// scan ran to the top of every grid and the per-cluster MVA reached AUC 0.9996
// for that reason alone.
//
// The quantile is NOT drawn from an RNG. It is a deterministic hash of the
// cluster itself (detId, quantized position, charge, plus a per-angle salt), so
// the angle is reproducible, independent of event order, and invariant under job
// splitting -- properties the old CLHEP-engine draw needed a carefully seeded
// per-event engine to achieve. Two salts give alpha and beta independent draws,
// matching the previous behaviour.
//
// ANGLE SOURCE "cnn" (opt-in; the default "payload" is everything above). The
// angle is INFERRED from the cluster's own digi ADC map by a CNN trained on
// TP-linked clusters (ngtagger-train eval_spixel_angles/cluster_cnn), run through
// ONNX. No truth enters: linked, unlinked, barrel and disc clusters take the same
// path, so "has an angle" and the angle distribution carry no link information
// (the payload's truth-keyed draw and the noise CDF are not used in this mode).
// The input window is built exactly as prep.py builds it: 24 x 16 pixels (local
// x rows, local y columns) centred on the bounding box with FLOOR division, value
// = digi ADC + 1 (0 = empty), plus the raw scalars named in cnn*Scalars; the
// ONNX graph applies the training normalisation itself. The sign/magnitude heads
// are decoded as cot = sign(logit) * tan(magnitude) - c, with c the fold offset
// (digitizer tan(theta_L), 0.202 on TBPX L2-L4, 0 on L1 and the discs).
// sigAlpha/sigBeta come from a STATIC table (cnnSigmaSet): the measured
// test-set residual width (q84 - q16)/2 in cot, per layer/disk and |reported cot|
// bin -- a resolution model keyed only on what the sensor reports, like the
// payload's sigma, not a per-cluster network output.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHit.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHitTruth.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsSurface.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include "correction.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
  // sizeY bucketing for the noise-angle CDF. MUST match SIZEY_BINS in
  // ngtagger-train/eval_spixel_angles/derive_noise_angle_payload.py, or the draw
  // is conditioned on a different variable than the one it was binned in.
  constexpr int kSizeYEdges[] = {1, 2, 3, 4, 5, 6, 8, 12};
  inline int sizeYBin(unsigned sy) {
    int b = 0;
    for (int e : kSizeYEdges) {
      if (static_cast<int>(sy) <= e)
        return b;
      ++b;
    }
    return b;
  }

  // splitmix64: deterministic uniform in [0,1) from the cluster's own identity.
  // Chosen over an RNG so the noise angle is reproducible and split-job invariant
  // without needing a seeded per-event engine.
  inline double hashUniform(uint32_t detId, float x, float y, float q, uint64_t salt) {
    uint64_t z = static_cast<uint64_t>(detId) * 0x9E3779B97F4A7C15ull;
    z ^= static_cast<uint64_t>(static_cast<int64_t>(std::llround(x * 1e4))) * 0xBF58476D1CE4E5B9ull;
    z ^= static_cast<uint64_t>(static_cast<int64_t>(std::llround(y * 1e4))) * 0x94D049BB133111EBull;
    z ^= static_cast<uint64_t>(static_cast<int64_t>(std::llround(q))) * 0xD6E8FEB86659FD93ull;
    z += salt;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z = z ^ (z >> 31);
    return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);  // [0,1)
  }

  // CNN input window, MUST match WX, WY in cluster_cnn/prep.py
  constexpr int kWX = 24, kWY = 16;
  inline int floorDiv2(int a) { return a >= 0 ? a / 2 : -((1 - a) / 2); }

  // One CNN scalar input, parsed once from its norm.json column name.
  struct CnnScalar {
    enum Kind { Layer, Disk, SizeX, SizeY, Size, LogQ, AdcSum, TruncX, TruncY, NPxLost } kind;
    unsigned k = 0;
  };
  std::vector<CnnScalar> parseCnnScalars(const std::vector<std::string>& names) {
    static const std::map<std::string, CnnScalar::Kind> fixed = {
        {"sizeX", CnnScalar::SizeX}, {"sizeY", CnnScalar::SizeY}, {"size", CnnScalar::Size},
        {"logq", CnnScalar::LogQ}, {"adcSum", CnnScalar::AdcSum}, {"truncX", CnnScalar::TruncX},
        {"truncY", CnnScalar::TruncY}, {"nPxLost", CnnScalar::NPxLost}};
    std::vector<CnnScalar> out;
    for (const auto& n : names) {
      const auto it = fixed.find(n);
      if (it != fixed.end())
        out.push_back({it->second, 0});
      else if (n.size() > 1 && (n[0] == 'l' || n[0] == 'd') && n.find_first_not_of("0123456789", 1) == std::string::npos)
        out.push_back({n[0] == 'l' ? CnnScalar::Layer : CnnScalar::Disk, static_cast<unsigned>(std::stoul(n.substr(1)))});
      else
        throw cms::Exception("Configuration") << "unknown CNN scalar input '" << n << "'";
    }
    return out;
  }
}  // namespace

// the ONNX sessions, shared by all streams (external linkage: it is a framework member type)
struct SmartPixelsCnnModels {
  std::unique_ptr<cms::Ort::ONNXRuntime> barrel, disc;
};

class SmartPixelsRecHitProducer : public edm::stream::EDProducer<edm::GlobalCache<SmartPixelsCnnModels>> {
public:
  SmartPixelsRecHitProducer(const edm::ParameterSet&, const SmartPixelsCnnModels*);
  static std::unique_ptr<SmartPixelsCnnModels> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const SmartPixelsCnnModels*) {}
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  void setCnnAngles(std::vector<SmartPixelsRecHit>&, const edm::DetSet<PixelDigi>*, unsigned lay, unsigned disk) const;

  const edm::EDGetTokenT<SiPixelRecHitCollection> recHitToken_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> simLinkToken_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;

  const unsigned maxLayer_;
  const unsigned maxDisk_;
  const double clusterMergeFrac_;
  const double measAngleMaxAbs_;
  const bool useCnn_;
  std::string angleSetPath_;

  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> digiToken_;
  std::vector<CnnScalar> cnnBarrelScalars_, cnnDiscScalars_;
  std::vector<float> cnnFoldOffset_;
  std::unique_ptr<correction::CorrectionSet> cnnSigmaSet_;
  correction::Correction::Ref corrCnnAlphaSigma_, corrCnnBetaSigma_;

  std::unique_ptr<correction::CorrectionSet> angleSet_;
  correction::Correction::Ref corrAlphaSigma_, corrBetaSigma_;
  correction::Correction::Ref corrValidProb_, corrValidFlat_;
  correction::CompoundCorrection::Ref corrAlphaShift_, corrBetaShift_;
  std::unique_ptr<correction::CorrectionSet> noiseSet_;
  correction::Correction::Ref corrNoiseCotAlpha_, corrNoiseCotBeta_;
};

std::unique_ptr<SmartPixelsCnnModels> SmartPixelsRecHitProducer::initializeGlobalCache(const edm::ParameterSet& cfg) {
  auto models = std::make_unique<SmartPixelsCnnModels>();
  if (cfg.getParameter<std::string>("angleSource") == "cnn") {
    models->barrel = std::make_unique<cms::Ort::ONNXRuntime>(
        edm::FileInPath(cfg.getParameter<std::string>("cnnBarrelModel")).fullPath());
    if (smartpixels::discScopeMaxDisk(cfg.getParameter<std::string>("discScope")) > 0)
      models->disc = std::make_unique<cms::Ort::ONNXRuntime>(
          edm::FileInPath(cfg.getParameter<std::string>("cnnDiscModel")).fullPath());
  }
  return models;
}

SmartPixelsRecHitProducer::SmartPixelsRecHitProducer(const edm::ParameterSet& cfg, const SmartPixelsCnnModels*)
    : recHitToken_(consumes<SiPixelRecHitCollection>(cfg.getParameter<edm::InputTag>("pixelRecHits"))),
      simLinkToken_(
          consumes<edm::DetSetVector<PixelDigiSimLink>>(cfg.getParameter<edm::InputTag>("pixelDigiSimLink"))),
      tpToken_(consumes<std::vector<TrackingParticle>>(cfg.getParameter<edm::InputTag>("trackingParticles"))),
      geomToken_(esConsumes()),
      topoToken_(esConsumes()),
      fieldToken_(esConsumes()),
      maxLayer_(cfg.getParameter<unsigned>("maxLayer")),
      maxDisk_(smartpixels::discScopeMaxDisk(cfg.getParameter<std::string>("discScope"))),
      clusterMergeFrac_(cfg.getParameter<double>("clusterMergeFrac")),
      measAngleMaxAbs_(cfg.getParameter<double>("measAngleMaxAbs")),
      useCnn_(cfg.getParameter<std::string>("angleSource") == "cnn") {
  const std::string source = cfg.getParameter<std::string>("angleSource");
  if (source != "payload" && source != "cnn")
    throw cms::Exception("Configuration") << "angleSource must be 'payload' or 'cnn', got '" << source << "'";
  const std::string p = cfg.getParameter<std::string>("angleSet");
  const std::string np = cfg.getParameter<std::string>("noiseSet");
  if (useCnn_) {
    // A payload configured alongside the CNN would silently not be used; refuse it.
    if (!p.empty() || !np.empty())
      throw cms::Exception("Configuration")
          << "angleSource 'cnn' takes no angleSet/noiseSet: every cluster gets the CNN angle";
    digiToken_ = consumes<edm::DetSetVector<PixelDigi>>(cfg.getParameter<edm::InputTag>("pixelDigis"));
    cnnBarrelScalars_ = parseCnnScalars(cfg.getParameter<std::vector<std::string>>("cnnBarrelScalars"));
    cnnDiscScalars_ = parseCnnScalars(cfg.getParameter<std::vector<std::string>>("cnnDiscScalars"));
    for (double c : cfg.getParameter<std::vector<double>>("cnnFoldOffset"))
      cnnFoldOffset_.push_back(static_cast<float>(c));
    if (cnnFoldOffset_.size() != maxLayer_ + 1)
      throw cms::Exception("Configuration") << "cnnFoldOffset needs maxLayer + 1 entries (index 0 = discs)";
    cnnSigmaSet_ = correction::CorrectionSet::from_file(
        edm::FileInPath(cfg.getParameter<std::string>("cnnSigmaSet")).fullPath());
    corrCnnAlphaSigma_ = cnnSigmaSet_->at("spix_cnn_alpha_sigma");
    corrCnnBetaSigma_ = cnnSigmaSet_->at("spix_cnn_beta_sigma");
  } else {
    if (p.empty())
      throw cms::Exception("Configuration")
          << "SmartPixelsRecHitProducer requires angleSet (the PixelAV angle-response payload). "
             "Without it there is no sensor response and the product would carry truth angles "
             "labelled as measurements.";
    angleSetPath_ = edm::FileInPath(p).fullPath();
    angleSet_ = correction::CorrectionSet::from_file(angleSetPath_);
    corrAlphaSigma_ = angleSet_->at("spix_angle_alpha_sigma");
    corrBetaSigma_ = angleSet_->at("spix_angle_beta_sigma");
    corrValidProb_ = angleSet_->at("spix_angle_valid_prob");
    corrValidFlat_ = angleSet_->at("spix_angle_valid_flat");
    corrAlphaShift_ = angleSet_->compound().at("spix_angle_alpha_shift");
    corrBetaShift_ = angleSet_->compound().at("spix_angle_beta_shift");
    if (!np.empty()) {
      noiseSet_ = correction::CorrectionSet::from_file(edm::FileInPath(np).fullPath());
      corrNoiseCotAlpha_ = noiseSet_->at("smarthit_noise_cotAlpha");
      corrNoiseCotBeta_ = noiseSet_->at("smarthit_noise_cotBeta");
    }
  }

  produces<SmartPixelsRecHitCollection>();
  produces<SmartPixelsRecHitTruthCollection>();
}

void SmartPixelsRecHitProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& topo = iSetup.getData(topoToken_);
  const auto& geom = iSetup.getData(geomToken_);
  const auto& field = iSetup.getData(fieldToken_);

  const auto& recHits = iEvent.get(recHitToken_);
  const auto& simLinks = iEvent.get(simLinkToken_);
  const edm::DetSetVector<PixelDigi>* digis = useCnn_ ? &iEvent.get(digiToken_) : nullptr;
  edm::Handle<std::vector<TrackingParticle>> tps;
  iEvent.getByToken(tpToken_, tps);

  // (eventId, SimTrackId) -> TP index. A TP owns several g4Tracks; mapping every
  // one of them onto the TP's index makes "same particle?" an integer compare.
  std::map<std::pair<unsigned int, unsigned int>, int> tpIndex;
  for (size_t i = 0; i < tps->size(); ++i) {
    const unsigned int evt = (*tps)[i].eventId().rawId();
    for (const auto& g4 : (*tps)[i].g4Tracks())
      tpIndex.emplace(std::make_pair(evt, g4.trackId()), static_cast<int>(i));
  }

  auto outHits = std::make_unique<SmartPixelsRecHitCollection>();
  auto outTruth = std::make_unique<SmartPixelsRecHitTruthCollection>();

  for (const auto& dsv : recHits) {
    const DetId did(dsv.detId());
    const unsigned surface = smartpixels::surfaceCode(did, topo, maxLayer_, maxDisk_);
    if (surface == 0)
      continue;
    const bool isDisc = smartpixels::isDiscSurface(surface);
    // lay = 0 marks a disc: the PixelAV payload and the noise CDF are TBPX-layer tables.
    const unsigned lay = isDisc ? 0 : surface;
    const unsigned disk = isDisc ? surface - smartpixels::kDiscSurfaceOffset : 0;
    const auto* pdu = dynamic_cast<const PixelGeomDetUnit*>(geom.idToDet(did));
    if (pdu == nullptr)
      continue;

    // Local B field at the MODULE CENTRE -- the same point the payload was
    // indexed on, and independent of both the track and the cluster position, so
    // the smear draw is reproducible outside any fit.
    const double bLocalY = pdu->toLocal(field.inTesla(pdu->position())).y();
    const double bz = field.inTesla(GlobalPoint(0, 0, 0)).z();

    std::map<unsigned int, const PixelDigiSimLink*> linkByChannel;
    const auto dsl = simLinks.find(did);
    if (dsl != simLinks.end())
      for (const auto& lk : *dsl) {
        auto it = linkByChannel.find(lk.channel());
        if (it == linkByChannel.end() || it->second->fraction() < lk.fraction())
          linkByChannel[lk.channel()] = &lk;
      }

    SmartPixelsRecHitCollection::FastFiller hitFill(*outHits, dsv.detId());
    SmartPixelsRecHitTruthCollection::FastFiller truthFill(*outTruth, dsv.detId());
    // held back until the module is complete: the CNN runs once per module
    std::vector<SmartPixelsRecHit> detHits;
    std::vector<SmartPixelsRecHitTruth> detTruth;

    uint16_t idxInDet = 0;
    for (const auto& rh : dsv) {
      const SiPixelCluster* cl = rh.cluster().isNonnull() ? &(*rh.cluster()) : nullptr;
      if (cl == nullptr)
        continue;

      SmartPixelsRecHit hit(rh.localPosition(),
                            rh.localPositionError(),
                            did.rawId(),
                            static_cast<uint16_t>(cl->sizeX()),
                            static_cast<uint16_t>(cl->sizeY()),
                            static_cast<uint16_t>(cl->size()),
                            static_cast<float>(cl->charge()),
                            rh.cluster());
      SmartPixelsRecHitTruth truth;

      // ---- dominant charge contributor, aggregated by TRACKINGPARTICLE ------
      // Keyed on the TP index, not the SimTrack id: a TP whose charge splits
      // across two of its own g4Tracks must not lose to a single-track TP.
      std::map<int, double> qByTp;
      double qTot = 0., qNoTp = 0.;
      for (const auto& px : cl->pixels()) {
        qTot += px.adc;
        const auto lit =
            linkByChannel.find(PixelDigi::pixelToChannel(static_cast<int>(px.x), static_cast<int>(px.y)));
        if (lit == linkByChannel.end()) {
          qNoTp += px.adc;
          continue;
        }
        const auto tit = tpIndex.find({lit->second->eventId().rawId(), lit->second->SimTrackId()});
        if (tit == tpIndex.end()) {
          qNoTp += px.adc;
          continue;
        }
        qByTp[tit->second] += px.adc;  // winner-takes-pixel: measured to matter for <0.7% of clusters
      }
      double qDom = 0., qSecond = 0.;
      int domTp = -1;
      for (const auto& kv : qByTp) {
        if (kv.second > qDom) {
          qSecond = qDom;
          qDom = kv.second;
          domTp = kv.first;
        } else if (kv.second > qSecond) {
          qSecond = kv.second;
        }
      }

      float trueCotA = -999.f, trueCotB = -999.f;
      if (domTp >= 0) {
        const TrackingParticle& tp = (*tps)[domTp];
        // ---- helix propagation to THIS hit -------------------------------
        // Transverse direction at a point on a helix is perpendicular to the
        // vector from the circle centre to that point; p_z is unchanged.
        const double q = tp.charge();
        const double px0 = tp.px(), py0 = tp.py(), pz0 = tp.pz();
        const double ptv = std::hypot(px0, py0);
        const GlobalPoint gh = pdu->toGlobal(rh.localPosition());
        double dirx = px0, diry = py0;
        if (q != 0. && ptv > 1e-6 && std::abs(bz) > 1e-6) {
          // signed radius [cm]; 0.2998 GeV/(T*m) -> /100 for cm
          const double R = ptv / (0.29979246 * std::abs(bz)) * 100.0;
          const double sgn = (q * bz > 0.) ? +1.0 : -1.0;
          // centre is at 90 deg from p, on the side set by the charge/field sign
          const double cx = tp.vx() + sgn * R * (py0 / ptv);
          const double cy = tp.vy() - sgn * R * (px0 / ptv);
          const double rx = gh.x() - cx, ry = gh.y() - cy;
          const double rn = std::hypot(rx, ry);
          if (rn > 1e-6) {
            // rotate the radius vector by -/+90 deg to get the direction of travel
            dirx = -sgn * (-ry) / rn * ptv;
            diry = -sgn * (rx) / rn * ptv;
            // preserve |pT|: components above are already scaled by ptv
          }
        }
        const LocalVector plv = pdu->toLocal(GlobalVector(dirx, diry, pz0));
        const double ppz = (std::abs(plv.z()) > 1e-9) ? plv.z() : 1e-9;
        trueCotA = static_cast<float>(plv.x() / ppz);
        trueCotB = static_cast<float>(plv.y() / ppz);

        // ---- sensor response: PixelAV validity gate + one shift per angle ---
        const std::vector<std::variant<int, double, std::string>> pin = {
            static_cast<int>(lay), static_cast<double>(trueCotA), static_cast<double>(trueCotB), bLocalY};
        if (!useCnn_ && !isDisc && corrValidFlat_->evaluate(pin) < corrValidProb_->evaluate(pin)) {
          const std::vector<std::variant<int, double, std::string>> pinAcc = {
              static_cast<int>(lay), static_cast<double>(trueCotA), static_cast<double>(trueCotB),
              bLocalY, 1.0};
          double cotA = trueCotA + corrAlphaShift_->evaluate(pinAcc);
          double cotB = trueCotB + corrBetaShift_->evaluate(pinAcc);
          // SIGMA IS LOOKED UP AT THE RECONSTRUCTED ANGLES, not the true ones.
          // The validity gate and the shift above must key on truth: they model
          // what a sensor does to a real incident track, and the reco angle does
          // not exist until the shift has been applied. But this sigma is
          // PUBLISHED, and every consumer treats it as a measurement uncertainty
          // -- the refit weights hits by it, and the seeding study derives its z0
          // search window from it as sigma(z0) = r * sigDirCotTheta. Keying it on
          // truth made a published uncertainty depend on information no trigger
          // can ever have, which is exactly the kind of presumption that cannot
          // be put into hardware. A lookup on (layer, cotAlpha, cotBeta, bLocalY)
          // CAN be: the payload tables are 60 values per layer, 240 in total.
          //
          // The noise-cluster branch below already keyed on its drawn angles, so
          // the two paths in this producer previously disagreed about which angle
          // the sigma belonged to.
          //
          // The payload's flow is 'clamp', so a reco angle pushed outside the
          // parametrised range by the shift lands in the edge bin, not a throw.
          const std::vector<std::variant<int, double, std::string>> pinReco = {
              static_cast<int>(lay), cotA, cotB, bLocalY};
          const double sigA = corrAlphaSigma_->evaluate(pinReco);
          const double sigB = corrBetaSigma_->evaluate(pinReco);
          bool hasA = sigA > 0., hasB = sigB > 0.;
          // Grazing clamp lives HERE, not in the fit: a sensor physically cannot
          // report an angle beyond the bound, so it is a property of the hit.
          if (std::abs(cotA) > measAngleMaxAbs_)
            hasA = false;
          if (std::abs(cotB) > measAngleMaxAbs_)
            hasB = false;
          hit.setAngles(static_cast<float>(cotA), static_cast<float>(cotB),
                        static_cast<float>(sigA), static_cast<float>(sigB), hasA, hasB);
        }
      }

      if (domTp < 0 && !useCnn_ && !isDisc && corrNoiseCotAlpha_ && corrNoiseCotBeta_) {
        // No simlink on any pixel: the sensor still sees charge and still emits an
        // angle. Draw from the inclusive per-layer distribution so an unlinked
        // cluster looks like an arbitrary cluster rather than a flagged special
        // case -- the whole point, since "reports an angle" was otherwise a
        // perfect truth proxy.
        const float lx = rh.localPosition().x(), ly = rh.localPosition().y();
        const double qa = hashUniform(did.rawId(), lx, ly, cl->charge(), 0x5CA1AB1Eull);
        const double qb = hashUniform(did.rawId(), lx, ly, cl->charge(), 0xB16B00B5ull);
        // Conditioned on cluster LENGTH, not just layer: a real sensor infers the
        // angle from the charge pattern, so reported angle and shape are physically
        // linked (mean sizeY rises 1.36 -> 6.51 across |cotBeta| bins for real
        // clusters). A layer-only draw breaks that link and is itself a tell.
        const int syb = sizeYBin(rh.cluster()->sizeY());
        double cotA = corrNoiseCotAlpha_->evaluate({static_cast<int>(lay), syb, qa});
        double cotB = corrNoiseCotBeta_->evaluate({static_cast<int>(lay), syb, qb});
        const std::vector<std::variant<int, double, std::string>> pin = {
            static_cast<int>(lay), cotA, cotB, bLocalY};
        const double sigA = corrAlphaSigma_->evaluate(pin);
        const double sigB = corrBetaSigma_->evaluate(pin);
        bool hasA = sigA > 0., hasB = sigB > 0.;
        if (std::abs(cotA) > measAngleMaxAbs_)
          hasA = false;
        if (std::abs(cotB) > measAngleMaxAbs_)
          hasB = false;
        hit.setAngles(static_cast<float>(cotA), static_cast<float>(cotB),
                      static_cast<float>(sigA), static_cast<float>(sigB), hasA, hasB);
      }

      const float frac = (qTot > 0. && qDom > 0.) ? static_cast<float>(qDom / qTot) : -999.f;
      const bool merged = (qTot > 0.) && ((qSecond / qTot) > clusterMergeFrac_);
      truth.set(domTp >= 0 ? TrackingParticleRef(tps, domTp) : TrackingParticleRef(),
                frac, merged, trueCotA, trueCotB, idxInDet);

      detHits.push_back(hit);
      detTruth.push_back(truth);
      ++idxInDet;
    }
    if (useCnn_) {
      const auto ddi = digis->find(did);
      setCnnAngles(detHits, ddi == digis->end() ? nullptr : &(*ddi), lay, disk);
    }
    for (size_t j = 0; j < detHits.size(); ++j) {
      hitFill.push_back(detHits[j]);
      truthFill.push_back(detTruth[j]);
    }
  }

  // The two products must stay index-aligned or every truth join is silently
  // wrong; assert rather than document.
  if (outHits->dataSize() != outTruth->dataSize() || outHits->size() != outTruth->size())
    throw cms::Exception("SmartPixelsRecHitMisaligned")
        << "rechit/truth collections diverged: " << outHits->dataSize() << " vs "
        << outTruth->dataSize() << " entries, " << outHits->size() << " vs " << outTruth->size()
        << " DetSets. The truth join would be meaningless.";

  iEvent.put(std::move(outHits));
  iEvent.put(std::move(outTruth));
}

void SmartPixelsRecHitProducer::setCnnAngles(std::vector<SmartPixelsRecHit>& hits,
                                             const edm::DetSet<PixelDigi>* digis,
                                             unsigned lay,
                                             unsigned disk) const {
  if (hits.empty())
    return;
  const bool isDisc = lay == 0;
  const cms::Ort::ONNXRuntime& model = isDisc ? *globalCache()->disc : *globalCache()->barrel;
  const auto& scal = isDisc ? cnnDiscScalars_ : cnnBarrelScalars_;
  std::unordered_map<unsigned, int> adcByChannel;
  if (digis != nullptr)
    for (const auto& d : *digis)
      adcByChannel[d.channel()] = d.adc();

  const size_t n = hits.size(), ns = scal.size();
  cms::Ort::FloatArrays in(2);
  in[0].assign(n * kWX * kWY, 0.f);
  in[1].assign(n * ns, 0.f);
  for (size_t i = 0; i < n; ++i) {
    const SiPixelCluster& cl = *hits[i].cluster();
    const int offx = floorDiv2(kWX - cl.sizeX()), offy = floorDiv2(kWY - cl.sizeY());
    float adcSum = 0.f;
    int lost = 0;
    for (const auto& px : cl.pixels()) {
      const auto it = adcByChannel.find(PixelDigi::pixelToChannel(px.x, px.y));
      const int a1 = (it == adcByChannel.end() ? -1 : it->second) + 1;  // prep.py: ADC + 1, missing digi -> 0
      adcSum += a1;
      const int r = px.x - cl.minPixelRow() + offx, c = px.y - cl.minPixelCol() + offy;
      if (r >= 0 && r < kWX && c >= 0 && c < kWY)
        in[0][(i * kWX + r) * kWY + c] = static_cast<float>(std::clamp(a1, 1, 255));
      else
        ++lost;
    }
    for (size_t k = 0; k < ns; ++k) {
      float v = 0.f;
      switch (scal[k].kind) {
        case CnnScalar::Layer: v = (lay == scal[k].k); break;
        case CnnScalar::Disk: v = (disk == scal[k].k); break;
        case CnnScalar::SizeX: v = cl.sizeX(); break;
        case CnnScalar::SizeY: v = cl.sizeY(); break;
        case CnnScalar::Size: v = cl.size(); break;
        case CnnScalar::LogQ: v = std::log(std::max(static_cast<float>(cl.charge()), 1.f)); break;
        case CnnScalar::AdcSum: v = adcSum; break;
        case CnnScalar::TruncX: v = cl.sizeX() > kWX; break;
        case CnnScalar::TruncY: v = cl.sizeY() > kWY; break;
        case CnnScalar::NPxLost: v = lost; break;
      }
      in[1][i * ns + k] = v;
    }
  }
  const auto out = model.run({"image", "scalars"},
                             in,
                             {{static_cast<int64_t>(n), kWX, kWY}, {static_cast<int64_t>(n), static_cast<int64_t>(ns)}},
                             {"out"},
                             static_cast<int64_t>(n))[0];
  const int detKey = isDisc ? 100 + static_cast<int>(disk) : static_cast<int>(lay);
  for (size_t i = 0; i < n; ++i) {
    const float* o = &out[i * 6];
    const float cotA = (o[2] > 0.f ? 1.f : -1.f) * std::tan(o[0]) - cnnFoldOffset_[lay];
    const float cotB = (o[3] > 0.f ? 1.f : -1.f) * std::tan(o[1]);
    const double sigA = corrCnnAlphaSigma_->evaluate({detKey, std::abs(static_cast<double>(cotA))});
    const double sigB = corrCnnBetaSigma_->evaluate({detKey, std::abs(static_cast<double>(cotB))});
    hits[i].setAngles(cotA, cotB, static_cast<float>(sigA), static_cast<float>(sigB),
                      sigA > 0. && std::abs(cotA) <= measAngleMaxAbs_, sigB > 0. && std::abs(cotB) <= measAngleMaxAbs_);
  }
}

void SmartPixelsRecHitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pixelRecHits", edm::InputTag("spixPixelRecHits"));
  desc.add<edm::InputTag>("pixelDigiSimLink", edm::InputTag("simSiPixelDigis", "Pixel"));
  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<unsigned>("maxLayer", 4)->setComment("highest TBPX layer processed (SmartPixels scope is 1..4)");
  desc.add<std::string>("discScope", "sp")
      ->setComment("disc clusters emitted (with truth): none | sp = ID1-ID4 (TFPX disks 1-4, the smart-pixel "
                   "discs) | tfpx = TFPX 1-8 | all = TFPX + TEPX. The PixelAV payload and the noise CDF are "
                   "TBPX-only, so disc hits carry no payload angle; the CNN source gives them one");
  desc.add<double>("clusterMergeFrac", 0.1)
      ->setComment("runner-up charge share above which the cluster is flagged merged (TRUTH-ONLY)");
  desc.add<double>("measAngleMaxAbs", 12.0)
      ->setComment("a sensor cannot report |cot| beyond this; clears hasAlpha/hasBeta for that angle "
                   "only. A HIT property, which is why it lives here and not in the fit.");
  desc.add<std::string>("noiseSet", "")
      ->setComment("inverse-CDF payload giving an angle to clusters with NO simlink "
                   "(smarthit_noise_cotAlpha/Beta). Leaving it empty means unlinked clusters "
                   "report no angle at all, which makes 'has an angle' a perfect truth proxy "
                   "(measured 99.5% vs 0.0%) and silently inflates every angle-using study. "
                   "Derive with ngtagger-train/eval_spixel_angles/derive_noise_angle_payload.py.");
  desc.add<std::string>("angleSet", "")
      ->setComment("PixelAV angle-response payload (spix_angle_* corrections); REQUIRED for "
                   "angleSource 'payload', refused for 'cnn'");
  desc.add<std::string>("angleSource", "payload")
      ->setComment("'payload': truth-keyed PixelAV draw (+ noiseSet for unlinked); 'cnn': ONNX CNN on the "
                   "digi ADC window, identical for every cluster, no truth input");
  desc.add<edm::InputTag>("pixelDigis", edm::InputTag("simSiPixelDigis", "Pixel"))
      ->setComment("digis the clusters were built from (cnn: the ADC map)");
  desc.add<std::string>("cnnBarrelModel", "")->setComment("cnn: ONNX model for TBPX clusters (FileInPath)");
  desc.add<std::string>("cnnDiscModel", "")->setComment("cnn: ONNX model for disc clusters (FileInPath)");
  desc.add<std::vector<std::string>>("cnnBarrelScalars", {})
      ->setComment("cnn: scalar inputs of the barrel model, in norm.json 'cols' order");
  desc.add<std::vector<std::string>>("cnnDiscScalars", {})
      ->setComment("cnn: scalar inputs of the disc model, in norm.json 'cols' order");
  desc.add<std::vector<double>>("cnnFoldOffset", {0., 0., 0.202, 0.202, 0.202})
      ->setComment("cnn: tan(theta_L) the alpha head was folded by, index = TBPX layer (0 = discs); "
                   "MUST match models_def.C_LAYER");
  desc.add<std::string>("cnnSigmaSet", "")
      ->setComment("cnn: correctionlib spix_cnn_{alpha,beta}_sigma(detKey = layer | 100 + disk, |cot|), "
                   "the measured residual width used as sigAlpha/sigBeta");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(SmartPixelsRecHitProducer);
