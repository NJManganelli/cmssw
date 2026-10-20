#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsRefitV3_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsRefitV3_h

// The seedless refit v3: the port of the Python replay's calibrated arm INN-D with the two-layer hit
// choice (ngtagger-train eval_refitq/replay, the reference; v3cmssw/PLAN.md rev 1). Framework-free so
// that it is unit-testable against the replay's own numbers; the producer builds the geometry and the
// per-event clusters and calls Engine::run per OT track.
//
// Invariants (doc/RefitV3Design.md):
//  - angles enter only locally, as measured minus predicted module-frame cot: no pointing term;
//  - the refit never reads a cluster's class (T or C);
//  - every calibrated number comes from the versioned payload set (data/refit_v3/<set>/).
//
// Units and conventions are the replay's: helix a = (rInv [1/cm], phi0, tanL, z0 [cm], d0 [cm]) in the
// TTTrack convention (phi(s) = phi0 - rInv s, POCA (d0 sin phi0, -d0 cos phi0)), pT = kpt / |rInv|,
// local measurements (x, y, cotAlpha, cotBeta) on the module. Everything is double precision.

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace correction {
  class CorrectionSet;
  class Correction;
}  // namespace correction
class TrackerGeometry;
class TrackerTopology;

namespace smartpixels::v3 {

  using Vec4 = std::array<double, 4>;
  using Vec5 = std::array<double, 5>;
  using Mat4 = std::array<std::array<double, 4>, 4>;
  using Mat5 = std::array<std::array<double, 5>, 5>;
  using Jac = std::array<std::array<double, 5>, 4>;  // d(u, v, cotA, cotB) / d(a)

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
  constexpr double kInf = std::numeric_limits<double>::infinity();
  constexpr int kOTSurface = 5;  // "previously visited surface" before the first visit: the OT

  inline bool isDisc(int surface) { return surface > 10; }

  // ---------------------------------------------------------------- geometry
  struct Module {
    uint32_t detId = 0;
    int surface = 0;  // 1-4 IL1-IL4, 11-14 ID1-ID4
    std::array<double, 3> o{}, ex{}, ey{}, ez{};
    double hw = 0., hl = 0.;
  };

  // The IT modules and the per-surface constants the replay derives (ot_projection.load_geometry):
  // barrel R (mean module-centre radius) and halfZ (largest |z| of the module outlines); disc Z (mean |z|
  // of the centres), radial extent of the outlines and the largest half-diagonal. Modules sorted by detId.
  class Geometry {
  public:
    void add(const Module& m) { modules_.push_back(m); }
    void finalize();
    bool ready() const { return ready_; }

    const std::vector<Module>& modules() const { return modules_; }
    int index(uint32_t detId) const;  // -1 if absent
    bool has(int S) const { return bySurface_.count(S) && !bySurface_.at(S).empty(); }
    const std::vector<int>& bySurface(int S) const { return bySurface_.at(S); }
    const std::vector<int>& bySide(int S, int side) const { return bySide_.at(S * 2 + side); }
    double R(int S) const { return R_.at(S); }  // S = 0 gives 0 (the inner-material bound)
    double halfZ(int S) const { return halfZ_.at(S); }
    double Z(int S) const { return Z_.at(S); }
    std::pair<double, double> rext(int S) const { return rext_.at(S); }
    const std::vector<int>& surfaces() const { return surfaces_; }  // barrel first, then discs

  private:
    std::vector<Module> modules_;
    std::unordered_map<uint32_t, int> byId_;
    std::unordered_map<int, std::vector<int>> bySurface_, bySide_;
    std::unordered_map<int, double> R_, halfZ_, Z_;
    std::unordered_map<int, std::pair<double, double>> rext_;
    std::vector<int> surfaces_;
    bool ready_ = false;
  };

  // The IL1-IL4 and ID1-ID4 modules of the reconstruction geometry, as the geometry dumper writes them (origin and
  // local axes through GeomDet::toGlobal, half-width and half-length of the bounds), converted to double.
  void buildGeometry(Geometry& g, const TrackerGeometry& geom, const TrackerTopology& topo);

  // ---------------------------------------------------------------- helix (ot_projection)
  struct HelixPoint {
    std::array<double, 3> p, d;  // position, direction (cos phi, sin phi, tanL)
  };
  HelixPoint helixAt(const Vec5& a, double s);
  double cylinderS(const Vec5& a, double R);  // NaN if never reached or s <= 0
  double discS(const Vec5& a, double Z);      // the plane on the side the helix heads to
  double surfaceS(const Vec5& a, const Geometry& g, int S);

  struct Projection {
    double u = kNaN, v = kNaN, cotA = kNaN, cotB = kNaN, s = kNaN;
    bool ok = false;
  };
  // project_to_modules (newton "study", refine on discs for the perturbed solves)
  Projection projectToModule(const Vec5& a, double s0, const Module& m, bool refine);
  struct ProjectionCov {
    Projection p;
    Jac J{};    // unguarded (NaN where no side converged)
    Mat4 S{};   // J C J^T (unguarded)
  };
  // project_with_cov, central differences (oneSided = the known-bad: forward differences only)
  ProjectionCov projectWithCov(const Vec5& a, const Mat5& C, double s0, const Module& m, bool oneSided = false);

  // ---------------------------------------------------------------- the visit plan (ORP.surface_path)
  struct VisitPlan {
    std::vector<int> seq;   // visited surfaces in visit order (active, on the path)
    bool discGoing = false; // a disc on the path (before the scenario mask)
  };
  VisitPlan visitPlan(const Vec5& aOT, const Geometry& g, const std::array<bool, 23>& active, bool outsideIn,
                      double discMarginCm = 1.0);

  // ---------------------------------------------------------------- clusters
  struct Cluster {
    double x = 0., y = 0., cotA = 0., cotB = 0.;
    double sigX = 0., sigY = 0., sigA = 0., sigB = 0.;
    bool hasA = false, hasB = false;  // already capped at |cot| <= measAngleMaxAbs
    int sizeX = 0, sizeY = 0;
    int index = -1;                   // the cluster-table index (axis 0 of L1TSmartPixelsCluster)
  };
  // per geometry module index, the clusters in table order
  using EventClusters = std::vector<std::vector<Cluster>>;

  // ---------------------------------------------------------------- payloads
  class SeedModel {
  public:
    void load(const std::string& path);
    Mat5 dC(double absEta, double p) const;  // the correction at |eta| and |p|

  private:
    std::vector<double> edges_;
    std::vector<std::array<double, 36>> theta_;
    std::array<double, 5> scale_{};
    double etaCell_ = 0., etaLast_ = 0.;
    std::vector<double> pLo_, pHi_;
  };

  struct MaterialPlane {
    bool disc = false;
    double r = 0.;              // cylinder radius (barrel planes: <0 means "layer" below)
    int layer = 0;              // barrel plane at a layer's mean radius (1-4), else 0
    double halfZ = kInf;        // cylinder |z| extent (disc-going rule)
    double z = 0., rmin = 0., rmax = 0.;  // disc plane
    double xx0 = 0.;
  };
  struct Material {
    std::vector<MaterialPlane> barrel;      // the calibration planes (barrel-only tracks)
    std::vector<MaterialPlane> discGoing;   // path-ordered, extent-aware planes
    MaterialPlane il1Module, beamPipe;
    double materialScale = 1., discScale = 1., cylScale = 1., otInnerR = 22.;
    double c0 = 0.0136, c1 = 0.038;
    void load(const std::string& path);
  };

  class Payloads {
  public:
    // manifestPath: data/refit_v3/<set>/manifest.json; every file it lists is checked (md5) and loaded
    void load(const std::string& manifestPath);
    bool loaded() const { return loaded_; }
    const std::string& setId() const { return setId_; }

    double rU(int sizeX, double absCotAPred) const;
    double rV(int sizeY, double absCotBPred) const;
    int orientation(uint32_t detId) const;
    // the angle-table entry: mode 0 none, 1 (m - c sigma, w sigma), 2 the angle is unusable
    struct AngleEntry {
      int mode = 0;
      double c = 0., w = 1.;
    };
    AngleEntry barrelAngle(int surface, int sign, double key, bool beta) const;
    AngleEntry discAngle(int surface, double key, bool beta) const;
    double skipPrice(const std::string& dims, int code, bool first) const;

    SeedModel seed;
    Material material;

  private:
    std::shared_ptr<const correction::Correction> rU_, rV_, orient_, angBarrel_, angDisc_, skip_;
    std::vector<std::shared_ptr<const correction::CorrectionSet>> sets_;
    std::string setId_;
    bool loaded_ = false;
  };

  // the incidence keys the angle tables are indexed by (spdefs.key_from_local, spdefs.local_incidence_key)
  inline double barrelKey(int sign, double cotA, double cotB) { return sign * cotB / std::sqrt(1.0 + cotA * cotA); }
  inline double discKey(double cotA, double cotB) {
    return cotB == 0. ? 0. : std::copysign(std::sqrt(cotA * cotA + cotB * cotB), cotB);
  }
  // a table entry applied to a measured angle and its published sigma (mode 2: the angle becomes unusable)
  inline void applyAngleEntry(const Payloads::AngleEntry& e, double& m, double& sig) {
    if (e.mode == 1) {
      m = m - e.c * sig;
      sig = e.w * sig;
    } else if (e.mode == 2) {
      m = kNaN;
      sig = kNaN;
    }
  }

  // ---------------------------------------------------------------- configuration
  struct Config {
    double nsig = 4.0;
    bool tree2 = true;              // the two-layer hit choice; false: K = 1 (best only)
    bool seedSmooth = true;         // + the seed model's dC
    bool highland = true;           // the calibrated material (the only process noise v3 implements)
    bool innerMaterial = true;
    bool rtables = true;            // position errors from the R tables, else the published CPE errors
    bool angleTable = true;         // angle centre / width from the tables, else the published sigma
    bool innovationMetric = true;   // ranking: the innovation chi2 (else the producer's sum of (r / sigma)^2)
    // angle use: [consumer rank|update][angle alpha|beta][region barrel|disc]
    std::array<std::array<std::array<bool, 2>, 2>, 2> use{};
    int angleSteps = 0;             // 0 all, 1 first visit step only, 2 every step but the first
    int maxKFUpdates = 8;
    int minHits = 1;
    int otNPar = 5;
    double chi2UpdateGate = 2.0e6;
    double jacobianMaxAbs = 1.0e4;
    double predAngleMaxAbs = 12.0;
    bool outsideIn = true;
    double discMarginCm = 1.0;
    std::string knownBad;           // test-only planted defects (PLAN §3 item 7); "" in every preset
  };

  // ---------------------------------------------------------------- output
  struct StepRecord {
    int surface = -1;
    bool visited = false, accepted = false;
    uint32_t detId = 0;
    int cluster = -1;
    int ncand = 0;
    int kchoice = -1;               // the committed two-layer choice at its window's first step
    std::array<double, 4> chi2{};   // per-scalar update chi2 (x, y, cotA, cotB)
    std::array<double, 4> pull{{kNaN, kNaN, kNaN, kNaN}};
    std::array<bool, 4> applied{};
    double rankChi2 = kNaN;         // the chosen candidate's rank-consumer innovation chi2
    double skipPrice = kNaN;        // P of the window opened here
    std::array<double, 3> scores{{kNaN, kNaN, kNaN}};  // the committing scores of the window opened here
    Vec4 h0{{kNaN, kNaN, kNaN, kNaN}};                 // the chosen candidate's prediction
    std::array<double, 2> projSig{{kNaN, kNaN}};      // sqrt of the prediction's (u, v) variances
    double selMargin = kNaN;                           // runner-up minus best ranking value (the replay's rankMargin)
    bool fallback = false;                             // the window's nearest-module projection failed: every module
    std::array<bool, 2> has{};                         // the chosen cluster's usable angles
    Vec4 meas{{kNaN, kNaN, kNaN, kNaN}}, sig{{kNaN, kNaN, kNaN, kNaN}};
    // closure (double): the state before Q, the covariance after Q, the state after the update
    Vec5 aPre{};
    Mat5 CPred{};
    Vec5 aPost{};
    Mat5 CPost{};
  };

  struct TrackResult {
    bool covOK = false, refit = false, discGoing = false;
    int nacc = 0;
    Vec5 a{};
    Mat5 C{};
    Mat5 seedC{};
    std::vector<StepRecord> steps;
  };

  class Engine {
  public:
    Engine(const Config& cfg, const Payloads& pay, const Geometry& geo) : cfg_(cfg), pay_(pay), geo_(geo) {}
    // aOT, COT: the OT track's helix and covariance (TTTrack convention); kpt: pT = kpt / |rInv|
    TrackResult run(const Vec5& aOT, const Mat5& COT, double kpt, const std::array<bool, 23>& active,
                    const EventClusters& clusters) const;

    // pieces exposed for the unit tests
    static double highlandTheta2(double xx0, double p, double tanl, double c0, double c1);
    static double innovationChi2(const Vec4& r, const Mat4& S, const Vec4& sig, const std::array<bool, 4>& dims);
    static std::array<double, 3> tree2Scores(const std::array<double, 2>& x1,
                                             const std::array<std::array<double, 2>, 3>* x2, double P);

  private:
    struct Copy;
    struct Cand;
    void qStep(Copy& c, int S, bool discGoing, double kpt) const;
    void addMaterialBarrel(Vec5& a, Mat5& C, int S, int Lv, double kpt, bool inner) const;
    void chargePath(Vec5& a, Mat5& C, double sCur, double sPrev, const std::vector<MaterialPlane>& planes,
                    double kpt) const;
    void innerMaterial(Copy& c, bool discGoing, double kpt) const;
    std::vector<Cand> candidates(const Copy& c, int S, int step, const EventClusters& clusters, bool& fallback) const;
    void measurements(const Cand& k, const Cluster& cl, Vec4& meas, Vec4& sig, std::array<bool, 2>& has) const;
    bool uses(int consumer, int angle, int surface, int step) const;
    double gateStat(const Cand& k, const Cluster& cl, int step, int consumer) const;
    void apply(Copy& c, int step, const Cand& k, const Cluster& cl) const;
    std::string dimsKey(int step) const;

    const Config& cfg_;
    const Payloads& pay_;
    const Geometry& geo_;
  };

}  // namespace smartpixels::v3

#endif
