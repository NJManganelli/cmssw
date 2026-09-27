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
// ANGLE SOURCE "emu" (opt-in; the default "payload" is everything above): the
// SmartPixels-OUTPUT EMULATOR, a conservative stand-in for what the SmartPixels chip
// will report (ngtagger-train eval_spixel_angles/cluster_cnn, emulator v4, angle space;
// python reference emu_common.py). One procedure for every cluster, barrel IL1-IL4
// and discs ID1-ID4, linked or not; angles are atan(cot) in DEGREES until the output.
//  1. Pseudo-truth. CLASS T: the dominant TP is charged and its helix crosses the
//     module plane within emuTruthMaxHelixDist of the hit (the ntuple's hxDist; no pT
//     cut); pseudo = its helix-propagated angle (the truth angle above).
//     CLASS C (everything else): pseudo = CNN angle + a residual DRAWN from the A1 model
//     p(rA | z) p(rB | z, rA) (z = the CNN angles and sigmas, 7 cluster observables and
//     the surface; trained on class T of TRAINING events). Both are ONNX graphs per
//     region: the tiny CNN (inputs = the digi ADC window as prep.py builds it + raw
//     scalars, the graph normalises) and A1 (the draw is inside the graph, fed 4 hash
//     uniforms). They run on class-C clusters only. Saturated at +-(90 - eps).
//  2. Response: the emuDumpK nearest rows of the SmartPixels Conv1D_Full-2bit dump
//     (rescaled, CMS frame; emuLibrary) in (alpha, beta) of the pseudo-truth -- on a
//     disc the alpha key is atan(cotAlpha - 0.312), a Lorentz-offset LOOKUP KEY only --
//     give the bias (median residual) and ONE drawn row r, the same for alpha and beta:
//     out = pseudo + bias + emuInflation * (r - bias), saturated (satAlpha/satBeta
//     flag the clamp; never wrapped), reported as cot = tan(out) in the module frame.
//  Every draw is a hash of (run, event, detId, idxInDet) with a salt per draw, so the
//  output is reproducible and split-job invariant. sigAlpha/sigBeta come from a STATIC
//  table (emuSigmaSet) keyed on sensor-reported quantities only: the surface code, the
//  cluster size(s) (as the table declares; the size counts as sensor-reported) and the
//  |reported cot|. EVERY emu cluster reports both angles (hasAlpha = hasBeta = true): a
//  validity flag that depended on the value would tell class C from class T. The class
//  is TRUTH (SmartPixelsRecHitTruth::emuClass), for validation only.
//
// NO |cot| CAP HERE, for either source: hasAlpha/hasBeta only say whether the sensor
// reports the angle at all (payload: the PixelAV validity gate; emu: always). The cap
// is applied by the algorithms that need it (the refit: digiRefitMeasAngleMaxAbs).

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
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsKnn.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsSurface.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsTpHelix.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include "correction.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
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

  // ---- emu: the cluster hash, MUST match emu_common.hash_u -----------------------
  inline uint64_t mix64(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  constexpr uint64_t kSaltDump = 0x5A7E11175A7E1117ull;
  // the A1 draw: (bin, position in bin) for alpha, then for beta
  constexpr uint64_t kSaltA1[4] = {
      0xA2A2A2A2A2A2A201ull, 0xA2A2A2A2A2A2A202ull, 0xA2A2A2A2A2A2A203ull, 0xA2A2A2A2A2A2A204ull};
  // uniform in [0, 1) from the stable cluster key (run, event, detId, idxInDet)
  inline double emuHashU(uint64_t run, uint64_t event, uint32_t detId, uint16_t idxInDet, uint64_t salt) {
    uint64_t h = mix64(((static_cast<uint64_t>(detId) << 16) | idxInDet) + 0x9E3779B97F4A7C15ull);
    h = mix64(h ^ event);
    h = mix64(h ^ run);
    h = mix64(h ^ salt);
    return static_cast<double>(h >> 11) * (1.0 / 9007199254740992.0);
  }
  inline unsigned emuHashIndex(uint64_t run, uint64_t event, uint32_t detId, uint16_t idxInDet, uint64_t salt, unsigned k) {
    return std::min(static_cast<unsigned>(emuHashU(run, event, detId, idxInDet, salt) * k), k - 1);
  }

  constexpr double kDegPerRad = 180.0 / M_PI, kRadPerDeg = M_PI / 180.0;

  // CNN input window, MUST match WX, WY in cluster_cnn/prep.py
  constexpr int kWX = 24, kWY = 16;
  inline int floorDiv2(int a) { return a >= 0 ? a / 2 : -((1 - a) / 2); }

  // One CNN scalar input, parsed once from its norm.json column name.
  struct CnnScalar {
    enum Kind {
      Layer,
      Disk,
      SizeX,
      SizeY,
      Size,
      LogQ,
      AdcSum,
      TruncX,
      TruncY,
      NPxLost,
      Is3D,
      AxXz,
      AxYz,
      AxZz,
      AxXr,
      AxYr,
      AxZr,
      SideZ,
      Pu,
      Pv,
      Pw,
      GlobalR,
      GlobalZ,
      LocalX,
      LocalY
    } kind;
    unsigned k = 0;
  };
  std::vector<CnnScalar> parseCnnScalars(const std::vector<std::string>& names) {
    static const std::map<std::string, CnnScalar::Kind> fixed = {
        {"sizeX", CnnScalar::SizeX},     {"sizeY", CnnScalar::SizeY},     {"size", CnnScalar::Size},
        {"logq", CnnScalar::LogQ},       {"adcSum", CnnScalar::AdcSum},   {"truncX", CnnScalar::TruncX},
        {"truncY", CnnScalar::TruncY},   {"nPxLost", CnnScalar::NPxLost}, {"is3D", CnnScalar::Is3D},
        {"axXz", CnnScalar::AxXz},       {"axYz", CnnScalar::AxYz},       {"axZz", CnnScalar::AxZz},
        {"axXr", CnnScalar::AxXr},       {"axYr", CnnScalar::AxYr},       {"axZr", CnnScalar::AxZr},
        {"sideZ", CnnScalar::SideZ},     {"pu", CnnScalar::Pu},           {"pv", CnnScalar::Pv},
        {"pw", CnnScalar::Pw},           {"globalR", CnnScalar::GlobalR}, {"globalZ", CnnScalar::GlobalZ},
        {"localX", CnnScalar::LocalX},   {"localY", CnnScalar::LocalY}};
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

  // Module orientation / position inputs of the CNN, computed as the training did
  // (ntuple float columns, then v2_common.derived() in numpy float32, no fused
  // multiply-add), so the closure test compares like with like.
  struct ModuleFrame {
    float axXx, axXy, axXz, axYx, axYy, axYz, axZx, axZy, axZz;  // local axes in global coordinates
    float axXr, axYr, axZr;                                     // their module-centre radial projections
  };
  SPX_NO_FP_CONTRACT ModuleFrame moduleFrame(const Plane& surf) {
    SPX_FP_CONTRACT_OFF_BODY
    const GlobalVector ex = surf.toGlobal(LocalVector(1, 0, 0)), ey = surf.toGlobal(LocalVector(0, 1, 0)),
                       ez = surf.toGlobal(LocalVector(0, 0, 1));
    const GlobalPoint mc = surf.position();
    const float modX = mc.x(), modY = mc.y();
    const float r = std::hypot(modX, modY);
    const float rx = modX / r, ry = modY / r;
    ModuleFrame f{ex.x(), ex.y(), ex.z(), ey.x(), ey.y(), ey.z(), ez.x(), ez.y(), ez.z(), 0.f, 0.f, 0.f};
    f.axZr = f.axZx * rx + f.axZy * ry;
    f.axXr = f.axXx * rx + f.axXy * ry;
    f.axYr = f.axYx * rx + f.axYy * ry;
    return f;
  }
  // unit vector origin -> hit in the module frame (pu, pv, pw)
  SPX_NO_FP_CONTRACT std::array<float, 3> promptDir(const ModuleFrame& f, float dx, float dy, float dz) {
    SPX_FP_CONTRACT_OFF_BODY
    const float lx = dx * f.axXx + dy * f.axXy + dz * f.axXz;
    const float ly = dx * f.axYx + dy * f.axYy + dz * f.axYz;
    const float lz = dx * f.axZx + dy * f.axZy + dz * f.axZz;
    const float dn = std::sqrt(dx * dx + dy * dy + dz * dz);
    return {lx / dn, ly / dn, lz / dn};
  }
  // Orientation inputs are rounded to 1e-6 (numpy.round in double, as
  // emu_common.cnn_scalars): the raw values carry geometry rounding residues (~1e-16)
  // that differ across geometry/release/arch and must not reach the network.
  inline float roundOrientation(float v) { return static_cast<float>(std::nearbyint(static_cast<double>(v) * 1e6) / 1e6); }

  // The response in angle space, contraction-free as numpy: the MEDIAN of the dump
  // residuals, then out = pseudo + bias + k (r - bias).
  SPX_NO_FP_CONTRACT double emuRespond(double pseudo, double bias, double r, double k) {
    SPX_FP_CONTRACT_OFF_BODY
    return (pseudo + bias) + k * (r - bias);
  }
  SPX_NO_FP_CONTRACT double emuDiscAlphaKey(double pseudoDeg, double offset) {
    SPX_FP_CONTRACT_OFF_BODY
    return std::atan(std::tan(pseudoDeg * kRadPerDeg) + offset) * kDegPerRad;
  }
}  // namespace

// Everything the emu source loads once and shares read-only between streams, plus
// the job-level counters its summary reports.
struct SmartPixelsEmuCache {
  // A1 inputs: CNN alpha, beta, sigma alpha, sigma beta [deg], sizeX, sizeY, log charge,
  // r, |z|, localX, localY, then the surface one-hot (layer or disk 1-4)
  static constexpr unsigned kNA1 = 15;
  std::unique_ptr<cms::Ort::ONNXRuntime> barrel, disc, a1Barrel, a1Disc;
  std::vector<float> dumpRA, dumpRB;
  smartpixels::KnnTree dumpTree;
  double setupSeconds = 0.;

  // [surface][class C = 0, T = 1]
  mutable std::array<std::array<std::atomic<uint64_t>, 2>, smartpixels::kMaxSurfaceCode + 1> nClusters{}, nSatA{},
      nSatB{}, nPseudoSat{};
  // wall time per step [ns], summed over streams
  enum Step { kHelix, kCnnInput, kCnn, kA1, kDump, kOther, kNSteps };
  mutable std::array<std::atomic<uint64_t>, kNSteps> stepNs{};
};

class SmartPixelsRecHitProducer : public edm::stream::EDProducer<edm::GlobalCache<SmartPixelsEmuCache>> {
public:
  SmartPixelsRecHitProducer(const edm::ParameterSet&, const SmartPixelsEmuCache*);
  static std::unique_ptr<SmartPixelsEmuCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const SmartPixelsEmuCache*);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  // per-cluster emu inputs gathered in the module loop
  struct EmuHit {
    bool classT = false;
    float trueCotA = 0.f, trueCotB = 0.f;
    float hxDist = -1.f;  // -1: not evaluated (no charged dominant TP)
    GlobalPoint gp;
  };
  void setEmuAngles(std::vector<SmartPixelsRecHit>&,
                    std::vector<SmartPixelsRecHitTruth>&,
                    const std::vector<EmuHit>&,
                    const edm::DetSet<PixelDigi>*,
                    const PixelGeomDetUnit&,
                    unsigned surface,
                    unsigned side,
                    const edm::EventID&,
                    std::vector<double>* debug) const;

  const edm::EDGetTokenT<SiPixelRecHitCollection> recHitToken_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> simLinkToken_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;

  const unsigned maxLayer_;
  const unsigned maxDisk_;
  const double clusterMergeFrac_;
  const bool useEmu_;
  std::string angleSetPath_;

  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> digiToken_;
  std::vector<CnnScalar> emuBarrelScalars_, emuDiscScalars_;
  double emuInflation_ = 0., emuDiscAlphaOffset_ = 0., emuSat_ = 0.;
  unsigned emuDumpK_ = 0;
  double emuTruthMaxHelixDist_ = 0., emuHelixStep_ = 0.;
  bool emuDebug_ = false;
  // emuSigmaSet inputs, resolved by name once (the table may key on the cluster sizes)
  enum class SigmaKey { Surface, SizeX, SizeY, AbsCot };
  static std::vector<SigmaKey> sigmaKeys(const correction::Correction&);
  double emuSigma(const correction::Correction::Ref&,
                  const std::vector<SigmaKey>&,
                  unsigned surface,
                  const SmartPixelsRecHit&,
                  float cot) const;
  std::unique_ptr<correction::CorrectionSet> emuSigmaSet_;
  correction::Correction::Ref corrEmuAlphaSigma_, corrEmuBetaSigma_;
  std::vector<SigmaKey> emuAlphaSigmaKeys_, emuBetaSigmaKeys_;

  std::unique_ptr<correction::CorrectionSet> angleSet_;
  correction::Correction::Ref corrAlphaSigma_, corrBetaSigma_;
  correction::Correction::Ref corrValidProb_, corrValidFlat_;
  correction::CompoundCorrection::Ref corrAlphaShift_, corrBetaShift_;
  std::unique_ptr<correction::CorrectionSet> noiseSet_;
  correction::Correction::Ref corrNoiseCotAlpha_, corrNoiseCotBeta_;
};

std::unique_ptr<SmartPixelsEmuCache> SmartPixelsRecHitProducer::initializeGlobalCache(const edm::ParameterSet& cfg) {
  auto cache = std::make_unique<SmartPixelsEmuCache>();
  if (cfg.getParameter<std::string>("angleSource") != "emu")
    return cache;
  const auto t0 = std::chrono::steady_clock::now();
  const auto onnx = [&cfg](const char* name) {
    return std::make_unique<cms::Ort::ONNXRuntime>(edm::FileInPath(cfg.getParameter<std::string>(name)).fullPath());
  };
  cache->barrel = onnx("emuBarrelModel");
  cache->a1Barrel = onnx("emuBarrelA1");
  if (smartpixels::discScopeMaxDisk(cfg.getParameter<std::string>("discScope")) > 0) {
    cache->disc = onnx("emuDiscModel");
    cache->a1Disc = onnx("emuDiscA1");
  }

  const std::string lib = edm::FileInPath(cfg.getParameter<std::string>("emuLibrary")).fullPath();
  std::unique_ptr<TFile> f(TFile::Open(lib.c_str(), "READ"));
  if (!f || f->IsZombie())
    throw cms::Exception("Configuration") << "cannot open emuLibrary " << lib;
  // the SmartPixels dump: lookup key (tA, tB), residuals (rA, rB) [deg]
  auto* t = f->Get<TTree>("dump");
  if (t == nullptr)
    throw cms::Exception("Configuration") << "emuLibrary " << lib << " has no tree 'dump'";
  float tA, tB, rA, rB;
  t->SetBranchAddress("tA", &tA);
  t->SetBranchAddress("tB", &tB);
  t->SetBranchAddress("rA", &rA);
  t->SetBranchAddress("rB", &rB);
  std::vector<double> pts;
  pts.reserve(2 * t->GetEntries());
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    pts.push_back(tA);
    pts.push_back(tB);
    cache->dumpRA.push_back(rA);
    cache->dumpRB.push_back(rB);
  }
  cache->dumpTree = smartpixels::KnnTree(std::move(pts), 2);
  cache->setupSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  edm::LogPrint("SmartPixelsEmu") << "emu setup " << std::fixed << std::setprecision(2) << cache->setupSeconds
                                  << " s: ONNX sessions + " << cache->dumpTree.size() << " dump rows";
  return cache;
}

void SmartPixelsRecHitProducer::globalEndJob(const SmartPixelsEmuCache* cache) {
  uint64_t total = 0;
  for (const auto& s : cache->nClusters)
    total += s[0] + s[1];
  if (total == 0)
    return;
  // the clamp at +-(90 - eps) is a silent failure mode if nobody counts it
  std::ostringstream msg;
  msg << "emu summary (setup " << std::fixed << std::setprecision(2) << cache->setupSeconds << " s)\n"
      << " surface class   clusters  satAlpha  satBeta  pseudoSat\n";
  for (unsigned s = 0; s <= smartpixels::kMaxSurfaceCode; ++s)
    for (unsigned c = 0; c < 2; ++c) {
      const uint64_t n = cache->nClusters[s][c];
      if (n == 0)
        continue;
      msg << std::setw(8) << s << std::setw(6) << (c ? "T" : "C") << std::setw(11) << n << std::setw(10)
          << cache->nSatA[s][c] << std::setw(9) << cache->nSatB[s][c] << std::setw(11) << cache->nPseudoSat[s][c]
          << "\n";
    }
  static const char* steps[] = {
      "helix (class-T test)", "CNN inputs", "CNN (ONNX)", "A1 (inputs + ONNX)", "dump kNN", "response+bookkeeping"};
  msg << " time summed over streams [s]:";
  for (unsigned k = 0; k < SmartPixelsEmuCache::kNSteps; ++k)
    msg << " " << steps[k] << " " << std::setprecision(3) << cache->stepNs[k] * 1e-9 << ";";
  edm::LogPrint("SmartPixelsEmu") << msg.str();
}

SmartPixelsRecHitProducer::SmartPixelsRecHitProducer(const edm::ParameterSet& cfg, const SmartPixelsEmuCache*)
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
      useEmu_(cfg.getParameter<std::string>("angleSource") == "emu") {
  const std::string source = cfg.getParameter<std::string>("angleSource");
  if (source != "payload" && source != "emu")
    throw cms::Exception("Configuration") << "angleSource must be 'payload' or 'emu', got '" << source << "'";
  const std::string p = cfg.getParameter<std::string>("angleSet");
  const std::string np = cfg.getParameter<std::string>("noiseSet");
  if (useEmu_) {
    // A payload configured alongside the emulator would silently not be used; refuse it.
    if (!p.empty() || !np.empty())
      throw cms::Exception("Configuration")
          << "angleSource 'emu' takes no angleSet/noiseSet: every cluster gets the emulated angle";
    // the CNN and A1 models cover IL1-IL4 and ID1-ID4 only
    if (maxLayer_ > 4 || maxDisk_ > 4)
      throw cms::Exception("Configuration") << "angleSource 'emu' supports maxLayer <= 4 and discScope none|sp";
    digiToken_ = consumes<edm::DetSetVector<PixelDigi>>(cfg.getParameter<edm::InputTag>("pixelDigis"));
    emuBarrelScalars_ = parseCnnScalars(cfg.getParameter<std::vector<std::string>>("emuBarrelScalars"));
    emuDiscScalars_ = parseCnnScalars(cfg.getParameter<std::vector<std::string>>("emuDiscScalars"));
    emuInflation_ = cfg.getParameter<double>("emuInflation");
    emuDiscAlphaOffset_ = cfg.getParameter<double>("emuDiscAlphaOffset");
    emuSat_ = 90. - cfg.getParameter<double>("emuSatEpsDeg");
    emuDumpK_ = cfg.getParameter<unsigned>("emuDumpK");
    if (emuDumpK_ == 0 || emuDumpK_ > 64)
      throw cms::Exception("Configuration") << "emuDumpK must be 1..64";
    emuTruthMaxHelixDist_ = cfg.getParameter<double>("emuTruthMaxHelixDist");
    emuHelixStep_ = cfg.getParameter<double>("emuHelixStep");
    emuDebug_ = cfg.getParameter<bool>("emuDebug");
    emuSigmaSet_ =
        correction::CorrectionSet::from_file(edm::FileInPath(cfg.getParameter<std::string>("emuSigmaSet")).fullPath());
    corrEmuAlphaSigma_ = emuSigmaSet_->at("spix_emu_alpha_sigma");
    corrEmuBetaSigma_ = emuSigmaSet_->at("spix_emu_beta_sigma");
    emuAlphaSigmaKeys_ = sigmaKeys(*corrEmuAlphaSigma_);
    emuBetaSigmaKeys_ = sigmaKeys(*corrEmuBetaSigma_);
    if (emuDebug_)
      produces<std::vector<double>>("emuDebug");
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
  const edm::DetSetVector<PixelDigi>* digis = useEmu_ ? &iEvent.get(digiToken_) : nullptr;
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
  auto emuDebug = emuDebug_ ? std::make_unique<std::vector<double>>() : nullptr;
  uint64_t helixNs = 0;

  for (const auto& dsv : recHits) {
    const DetId did(dsv.detId());
    const unsigned surface = smartpixels::surfaceCode(did, topo, maxLayer_, maxDisk_);
    if (surface == 0)
      continue;
    const bool isDisc = smartpixels::isDiscSurface(surface);
    // lay = 0 marks a disc: the PixelAV payload and the noise CDF are TBPX-layer tables.
    const unsigned lay = isDisc ? 0 : surface;
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
    // held back until the module is complete: the emulator's CNN runs once per module
    std::vector<SmartPixelsRecHit> detHits;
    std::vector<SmartPixelsRecHitTruth> detTruth;
    std::vector<EmuHit> emuHits;

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
      EmuHit emu;
      const GlobalPoint gh = pdu->toGlobal(rh.localPosition());

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

        if (useEmu_ && q != 0.) {
          // class T needs the helix to pass the hit: the plane crossing nearest it,
          // exactly the cluster ntuple's hxDist (stopping at the first crossing that
          // decides it, unless the debug product wants the full scan's value)
          const auto t0 = std::chrono::steady_clock::now();
          const smartpixels::TpHelix hx(tp, bz);
          const auto hc = smartpixels::helixPlaneCrossing(
              hx, pdu->surface(), gh, emuHelixStep_, emuDebug_ ? -1. : emuTruthMaxHelixDist_);
          helixNs += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t0).count();
          emu.hxDist = hc.s >= 0 ? static_cast<float>(hc.dist) : -1.f;
          emu.classT = hc.s >= 0 && static_cast<float>(hc.dist) < emuTruthMaxHelixDist_;
        }

        // ---- sensor response: PixelAV validity gate + one shift per angle ---
        const std::vector<std::variant<int, double, std::string>> pin = {
            static_cast<int>(lay), static_cast<double>(trueCotA), static_cast<double>(trueCotB), bLocalY};
        if (!useEmu_ && !isDisc && corrValidFlat_->evaluate(pin) < corrValidProb_->evaluate(pin)) {
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
          hit.setAngles(static_cast<float>(cotA), static_cast<float>(cotB),
                        static_cast<float>(sigA), static_cast<float>(sigB), sigA > 0., sigB > 0.);
        }
      }

      if (domTp < 0 && !useEmu_ && !isDisc && corrNoiseCotAlpha_ && corrNoiseCotBeta_) {
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
        hit.setAngles(static_cast<float>(cotA), static_cast<float>(cotB),
                      static_cast<float>(sigA), static_cast<float>(sigB), sigA > 0., sigB > 0.);
      }

      const float frac = (qTot > 0. && qDom > 0.) ? static_cast<float>(qDom / qTot) : -999.f;
      const bool merged = (qTot > 0.) && ((qSecond / qTot) > clusterMergeFrac_);
      truth.set(domTp >= 0 ? TrackingParticleRef(tps, domTp) : TrackingParticleRef(),
                frac, merged, trueCotA, trueCotB, idxInDet);

      detHits.push_back(hit);
      detTruth.push_back(truth);
      if (useEmu_) {
        emu.trueCotA = trueCotA;
        emu.trueCotB = trueCotB;
        emu.gp = gh;
        emuHits.push_back(emu);
      }
      ++idxInDet;
    }
    if (useEmu_) {
      const auto ddi = digis->find(did);
      setEmuAngles(detHits,
                   detTruth,
                   emuHits,
                   ddi == digis->end() ? nullptr : &(*ddi),
                   *pdu,
                   surface,
                   isDisc ? topo.pxfSide(did) : 0,
                   iEvent.id(),
                   emuDebug.get());
    }
    for (size_t j = 0; j < detHits.size(); ++j) {
      hitFill.push_back(detHits[j]);
      truthFill.push_back(detTruth[j]);
    }
  }
  if (useEmu_)
    globalCache()->stepNs[SmartPixelsEmuCache::kHelix] += helixNs;

  // The two products must stay index-aligned or every truth join is silently
  // wrong; assert rather than document.
  if (outHits->dataSize() != outTruth->dataSize() || outHits->size() != outTruth->size())
    throw cms::Exception("SmartPixelsRecHitMisaligned")
        << "rechit/truth collections diverged: " << outHits->dataSize() << " vs "
        << outTruth->dataSize() << " entries, " << outHits->size() << " vs " << outTruth->size()
        << " DetSets. The truth join would be meaningless.";

  iEvent.put(std::move(outHits));
  iEvent.put(std::move(outTruth));
  if (emuDebug)
    iEvent.put(std::move(emuDebug), "emuDebug");
}

// emuDebug product: kEmuDebugStride doubles per cluster, in rec-hit collection order.
// MUST match kEmuDebugCols in SmartPixelsClusterAdcNtuplizer.
namespace {
  constexpr unsigned kEmuDebugStride = 15;
}

std::vector<SmartPixelsRecHitProducer::SigmaKey> SmartPixelsRecHitProducer::sigmaKeys(const correction::Correction& c) {
  static const std::map<std::string, SigmaKey> known = {{"surface", SigmaKey::Surface},
                                                        {"sizeX", SigmaKey::SizeX},
                                                        {"sizeY", SigmaKey::SizeY},
                                                        {"absCot", SigmaKey::AbsCot}};
  std::vector<SigmaKey> keys;
  for (const auto& v : c.inputs()) {
    const auto it = known.find(v.name());
    if (it == known.end())
      throw cms::Exception("Configuration") << "emuSigmaSet " << c.name() << " has input '" << v.name()
                                            << "'; only surface, sizeX, sizeY, absCot are sensor-reported";
    keys.push_back(it->second);
  }
  return keys;
}

double SmartPixelsRecHitProducer::emuSigma(const correction::Correction::Ref& corr,
                                           const std::vector<SigmaKey>& keys,
                                           unsigned surface,
                                           const SmartPixelsRecHit& hit,
                                           float cot) const {
  std::vector<correction::Variable::Type> in;
  in.reserve(keys.size());
  for (const auto k : keys) {
    switch (k) {
      case SigmaKey::Surface:
        in.emplace_back(static_cast<int>(surface));
        break;
      case SigmaKey::SizeX:
        in.emplace_back(static_cast<int>(hit.sizeX()));
        break;
      case SigmaKey::SizeY:
        in.emplace_back(static_cast<int>(hit.sizeY()));
        break;
      case SigmaKey::AbsCot:
        in.emplace_back(std::abs(static_cast<double>(cot)));
        break;
    }
  }
  return corr->evaluate(in);
}

void SmartPixelsRecHitProducer::setEmuAngles(std::vector<SmartPixelsRecHit>& hits,
                                             std::vector<SmartPixelsRecHitTruth>& truths,
                                             const std::vector<EmuHit>& emu,
                                             const edm::DetSet<PixelDigi>* digis,
                                             const PixelGeomDetUnit& pdu,
                                             unsigned surface,
                                             unsigned side,
                                             const edm::EventID& eid,
                                             std::vector<double>* debug) const {
  if (hits.empty())
    return;
  using Clock = std::chrono::steady_clock;
  const auto ns = [](Clock::time_point a, Clock::time_point b) {
    return static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count());
  };
  const SmartPixelsEmuCache& cache = *globalCache();
  const bool isDisc = smartpixels::isDiscSurface(surface);
  const unsigned lay = isDisc ? 0 : surface;
  const unsigned disk = isDisc ? surface - smartpixels::kDiscSurfaceOffset : 0;
  const unsigned layerOrDisk = isDisc ? disk : lay;
  const cms::Ort::ONNXRuntime& model = isDisc ? *cache.disc : *cache.barrel;
  const cms::Ort::ONNXRuntime& a1 = isDisc ? *cache.a1Disc : *cache.a1Barrel;
  const auto& scal = isDisc ? emuDiscScalars_ : emuBarrelScalars_;
  const ModuleFrame frame = moduleFrame(pdu.surface());
  const bool is3D = !isDisc && lay == 1;
  const float discSideZ = side == 1 ? -1.f : 1.f;  // pxfSide 1 = -z

  const auto t0 = Clock::now();
  // the CNN and A1 are needed only for class C: class T takes the truth angle (the
  // debug product wants them everywhere, for the closure test)
  const size_t n = hits.size(), ns_ = scal.size();
  std::vector<int> cnnRow(n, -1);
  std::vector<size_t> cnnHits;
  for (size_t i = 0; i < n; ++i)
    if (!emu[i].classT || debug != nullptr) {
      cnnRow[i] = cnnHits.size();
      cnnHits.push_back(i);
    }
  std::unordered_map<unsigned, int> adcByChannel;
  if (digis != nullptr && !cnnHits.empty())
    for (const auto& d : *digis)
      adcByChannel[d.channel()] = d.adc();
  const size_t nc = cnnHits.size();
  cms::Ort::FloatArrays in(2);
  in[0].assign(nc * kWX * kWY, 0.f);
  in[1].assign(nc * ns_, 0.f);
  for (size_t ic = 0; ic < nc; ++ic) {
    const size_t i = cnnHits[ic];
    const SiPixelCluster& cl = *hits[i].cluster();
    const int offx = floorDiv2(kWX - cl.sizeX()), offy = floorDiv2(kWY - cl.sizeY());
    float adcSum = 0.f;
    int lost = 0;
    for (const auto& px : cl.pixels()) {
      const auto it = adcByChannel.find(PixelDigi::pixelToChannel(px.x, px.y));
      const int adc1 = (it == adcByChannel.end() ? -1 : it->second) + 1;  // prep.py: ADC + 1, missing digi -> 0
      adcSum += adc1;
      const int r = px.x - cl.minPixelRow() + offx, c = px.y - cl.minPixelCol() + offy;
      if (r >= 0 && r < kWX && c >= 0 && c < kWY)
        in[0][(ic * kWX + r) * kWY + c] = static_cast<float>(std::clamp(adc1, 1, 255));
      else
        ++lost;
    }
    const GlobalPoint& gp = emu[i].gp;
    const auto pd = promptDir(frame, gp.x(), gp.y(), gp.z());
    const float gz = gp.z();
    for (size_t k = 0; k < ns_; ++k) {
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
        case CnnScalar::Is3D: v = is3D; break;
        case CnnScalar::AxXz: v = roundOrientation(frame.axXz); break;
        case CnnScalar::AxYz: v = roundOrientation(frame.axYz); break;
        case CnnScalar::AxZz: v = roundOrientation(frame.axZz); break;
        case CnnScalar::AxXr: v = roundOrientation(frame.axXr); break;
        case CnnScalar::AxYr: v = roundOrientation(frame.axYr); break;
        case CnnScalar::AxZr: v = roundOrientation(frame.axZr); break;
        case CnnScalar::SideZ: v = isDisc ? discSideZ : (gz > 0.f ? 1.f : (gz < 0.f ? -1.f : 0.f)); break;
        case CnnScalar::Pu: v = pd[0]; break;
        case CnnScalar::Pv: v = pd[1]; break;
        case CnnScalar::Pw: v = pd[2]; break;
        case CnnScalar::GlobalR: v = std::hypot(gp.x(), gp.y()); break;
        case CnnScalar::GlobalZ: v = gz; break;
        case CnnScalar::LocalX: v = hits[i].localPosition().x(); break;
        case CnnScalar::LocalY: v = hits[i].localPosition().y(); break;
      }
      in[1][ic * ns_ + k] = v;
    }
  }
  const auto t1 = Clock::now();
  const std::vector<float> out =
      nc == 0 ? std::vector<float>()
              : model.run({"image", "scalars"},
                          in,
                          {{static_cast<int64_t>(nc), kWX, kWY}, {static_cast<int64_t>(nc), static_cast<int64_t>(ns_)}},
                          {"out"},
                          static_cast<int64_t>(nc))[0];
  const auto t2 = Clock::now();

  // CNN decode (atan(cot) [rad], log sigma [rad] -> deg; emu_common.decode), then the
  // A1 residual draw, fed the raw A1 inputs and the 4 hash uniforms
  constexpr unsigned kNA1 = SmartPixelsEmuCache::kNA1;
  std::vector<double> cnnA(nc), cnnB(nc), sigCnnA(nc), sigCnnB(nc);
  cms::Ort::FloatArrays a1in(2);
  a1in[0].assign(nc * kNA1, 0.f);
  a1in[1].assign(nc * 4, 0.f);
  for (size_t ic = 0; ic < nc; ++ic) {
    const size_t i = cnnHits[ic];
    const float* o = &out[ic * 6];
    cnnA[ic] = std::clamp(static_cast<double>(o[0]), -1.5707, 1.5707) * kDegPerRad;
    cnnB[ic] = std::clamp(static_cast<double>(o[1]), -1.5707, 1.5707) * kDegPerRad;
    sigCnnA[ic] = std::exp(std::clamp(static_cast<double>(o[2]), -7., 3.)) * kDegPerRad;
    sigCnnB[ic] = std::exp(std::clamp(static_cast<double>(o[3]), -7., 3.)) * kDegPerRad;
    const SiPixelCluster& cl = *hits[i].cluster();
    const GlobalPoint& gp = emu[i].gp;
    float* z = &a1in[0][ic * kNA1];
    z[0] = cnnA[ic];
    z[1] = cnnB[ic];
    z[2] = sigCnnA[ic];
    z[3] = sigCnnB[ic];
    z[4] = cl.sizeX();
    z[5] = cl.sizeY();
    z[6] = std::log(std::max(static_cast<float>(cl.charge()), 1.f));
    z[7] = std::hypot(gp.x(), gp.y());
    z[8] = std::abs(gp.z());
    z[9] = hits[i].localPosition().x();
    z[10] = hits[i].localPosition().y();
    for (unsigned k = 1; k <= 4; ++k)
      z[10 + k] = (layerOrDisk == k);
    for (unsigned k = 0; k < 4; ++k)
      a1in[1][ic * 4 + k] =
          static_cast<float>(emuHashU(eid.run(), eid.event(), hits[i].detId(), truths[i].idxInDet(), kSaltA1[k]));
  }
  const std::vector<float> draws =
      nc == 0 ? std::vector<float>()
              : a1.run({"z", "u"},
                       a1in,
                       {{static_cast<int64_t>(nc), kNA1}, {static_cast<int64_t>(nc), 4}},
                       {"r"},
                       static_cast<int64_t>(nc))[0];
  const auto t3 = Clock::now();

  const auto sat = [this](double a) { return std::clamp(a, -emuSat_, emuSat_); };
  std::vector<std::pair<double, uint32_t>> nbr;
  nbr.reserve(emuDumpK_);
  std::array<double, 64> rs;
  uint64_t tDump = 0;
  for (size_t i = 0; i < n; ++i) {
    const EmuHit& e = emu[i];
    const int ic = cnnRow[i];
    const uint16_t idx = truths[i].idxInDet();
    const uint32_t detId = hits[i].detId();

    // ---- 1. pseudo-truth [deg]
    double pA, pB;
    bool pseudoSat = false;
    if (e.classT) {
      pA = std::atan(static_cast<double>(e.trueCotA)) * kDegPerRad;
      pB = std::atan(static_cast<double>(e.trueCotB)) * kDegPerRad;
    } else {
      const double rawA = cnnA[ic] + static_cast<double>(draws[ic * 2]);
      const double rawB = cnnB[ic] + static_cast<double>(draws[ic * 2 + 1]);
      pseudoSat = std::abs(rawA) >= emuSat_ || std::abs(rawB) >= emuSat_;
      pA = sat(rawA);
      pB = sat(rawB);
    }

    // ---- 2. SmartPixels response from the dump rows nearest the pseudo-truth
    const auto tb = Clock::now();
    const double key[2] = {isDisc ? emuDiscAlphaKey(pA, emuDiscAlphaOffset_) : pA, pB};
    cache.dumpTree.query(key, emuDumpK_, nbr);
    tDump += ns(tb, Clock::now());
    const uint32_t row = nbr[emuHashIndex(eid.run(), eid.event(), detId, idx, kSaltDump, emuDumpK_)].second;
    const auto median = [&](const std::vector<float>& r) {
      const unsigned m = nbr.size();
      for (unsigned j = 0; j < m; ++j)
        rs[j] = r[nbr[j].second];
      std::sort(rs.begin(), rs.begin() + m);
      return m % 2 ? rs[m / 2] : (rs[m / 2 - 1] + rs[m / 2]) / 2.;
    };
    const double biasA = median(cache.dumpRA), biasB = median(cache.dumpRB);
    const double rawA = emuRespond(pA, biasA, cache.dumpRA[row], emuInflation_);
    const double rawB = emuRespond(pB, biasB, cache.dumpRB[row], emuInflation_);
    const bool satA = std::abs(rawA) >= emuSat_, satB = std::abs(rawB) >= emuSat_;
    const double outA = sat(rawA), outB = sat(rawB);
    const float cotA = static_cast<float>(std::tan(outA * kRadPerDeg));
    const float cotB = static_cast<float>(std::tan(outB * kRadPerDeg));

    const double sigA = emuSigma(corrEmuAlphaSigma_, emuAlphaSigmaKeys_, surface, hits[i], cotA);
    const double sigB = emuSigma(corrEmuBetaSigma_, emuBetaSigmaKeys_, surface, hits[i], cotB);
    if (!(sigA > 0. && sigB > 0.))
      throw cms::Exception("SmartPixelsEmu") << "emuSigmaSet gives a non-positive sigma (" << sigA << ", " << sigB
                                             << ") on surface " << surface << ": every emu angle must be reported";
    hits[i].setAngles(cotA, cotB, static_cast<float>(sigA), static_cast<float>(sigB), true, true);
    hits[i].setSaturation(satA, satB);
    truths[i].setEmuClass(e.classT ? SmartPixelsRecHitTruth::kEmuClassT : SmartPixelsRecHitTruth::kEmuClassC);

    const unsigned c = e.classT ? 1 : 0;
    ++cache.nClusters[surface][c];
    cache.nSatA[surface][c] += satA;
    cache.nSatB[surface][c] += satB;
    cache.nPseudoSat[surface][c] += pseudoSat;
    if (debug != nullptr) {
      const double v[kEmuDebugStride] = {cnnA[ic],
                                         cnnB[ic],
                                         sigCnnA[ic],
                                         sigCnnB[ic],
                                         static_cast<double>(c),
                                         static_cast<double>(e.hxDist),
                                         static_cast<double>(draws[ic * 2]),
                                         static_cast<double>(draws[ic * 2 + 1]),
                                         pA,
                                         pB,
                                         static_cast<double>(row),
                                         biasA,
                                         biasB,
                                         outA,
                                         outB};
      debug->insert(debug->end(), v, v + kEmuDebugStride);
    }
  }
  const auto t4 = Clock::now();
  cache.stepNs[SmartPixelsEmuCache::kCnnInput] += ns(t0, t1);
  cache.stepNs[SmartPixelsEmuCache::kCnn] += ns(t1, t2);
  cache.stepNs[SmartPixelsEmuCache::kA1] += ns(t2, t3);
  cache.stepNs[SmartPixelsEmuCache::kDump] += tDump;
  cache.stepNs[SmartPixelsEmuCache::kOther] += ns(t3, t4) - std::min(ns(t3, t4), tDump);
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
                   "TBPX-only, so disc hits carry no payload angle; the emu source gives them one (none|sp only)");
  desc.add<double>("clusterMergeFrac", 0.1)
      ->setComment("runner-up charge share above which the cluster is flagged merged (TRUTH-ONLY)");
  desc.add<std::string>("noiseSet", "")
      ->setComment("inverse-CDF payload giving an angle to clusters with NO simlink "
                   "(smarthit_noise_cotAlpha/Beta). Leaving it empty means unlinked clusters "
                   "report no angle at all, which makes 'has an angle' a perfect truth proxy "
                   "(measured 99.5% vs 0.0%) and silently inflates every angle-using study. "
                   "Derive with ngtagger-train/eval_spixel_angles/derive_noise_angle_payload.py.");
  desc.add<std::string>("angleSet", "")
      ->setComment("PixelAV angle-response payload (spix_angle_* corrections); REQUIRED for "
                   "angleSource 'payload', refused for 'emu'");
  desc.add<std::string>("angleSource", "payload")
      ->setComment("'payload': truth-keyed PixelAV draw (+ noiseSet for unlinked); 'emu': the SmartPixels-output "
                   "emulator (class T / class C pseudo-truth + dump response), every cluster, discs included");
  desc.add<edm::InputTag>("pixelDigis", edm::InputTag("simSiPixelDigis", "Pixel"))
      ->setComment("digis the clusters were built from (emu: the CNN's ADC map)");
  desc.add<std::string>("emuBarrelModel", "")->setComment("emu: tiny-CNN ONNX model for TBPX clusters (FileInPath)");
  desc.add<std::string>("emuDiscModel", "")->setComment("emu: tiny-CNN ONNX model for disc clusters (FileInPath)");
  desc.add<std::string>("emuBarrelA1", "")
      ->setComment("emu: A1 ONNX model (class-C residual draw: z[15], u[4] -> r[2] deg) for TBPX clusters");
  desc.add<std::string>("emuDiscA1", "")->setComment("emu: A1 ONNX model for disc clusters");
  desc.add<std::vector<std::string>>("emuBarrelScalars", {})
      ->setComment("emu: scalar inputs of the barrel model, in norm.json 'cols' order");
  desc.add<std::vector<std::string>>("emuDiscScalars", {})
      ->setComment("emu: scalar inputs of the disc model, in norm.json 'cols' order");
  desc.add<std::string>("emuLibrary", "")
      ->setComment("emu: ROOT file with the rescaled SmartPixels dump (tree dump); emu_build_library.py");
  desc.add<std::string>("emuSigmaSet", "")
      ->setComment("emu: correctionlib spix_emu_{alpha,beta}_sigma, the emulated residual width used as "
                   "sigAlpha/sigBeta; inputs taken by name from surface, sizeX, sizeY, absCot (|reported cot|)");
  desc.add<double>("emuInflation", 1.5)->setComment("emu: out = pseudo + bias + emuInflation (r - bias)");
  desc.add<double>("emuDiscAlphaOffset", -0.312)
      ->setComment("emu: disc alpha lookup key atan(cotAlpha + offset), the PixelAV Lorentz offset (key only)");
  desc.add<double>("emuSatEpsDeg", 0.01)->setComment("emu: angles saturate at +-(90 - eps) deg, never wrap");
  desc.add<unsigned>("emuDumpK", 32)->setComment("emu: dump rows per lookup (bias = their median residual)");
  desc.add<double>("emuTruthMaxHelixDist", 1.0)
      ->setComment("emu: class T needs the TP helix to cross the module plane within this of the hit [cm]");
  desc.add<double>("emuHelixStep", 0.2)->setComment("emu: helix-plane root scan step [cm], as the cluster ntuple");
  desc.add<bool>("emuDebug", false)
      ->setComment("emu: also put std::vector<double> 'emuDebug' (15 per cluster: CNN angles/sigmas, class, "
                   "hxDist from the full scan, A1 draws, pseudo-truth, dump row, biases, output angles; the CNN "
                   "and A1 then run on every cluster) for closure tests");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(SmartPixelsRecHitProducer);
