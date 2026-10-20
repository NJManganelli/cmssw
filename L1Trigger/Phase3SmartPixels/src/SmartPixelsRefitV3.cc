#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitV3.h"

#include "FWCore/Utilities/interface/Digest.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <correction.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>

namespace smartpixels::v3 {

  namespace {
    // ot_projection constants
    constexpr double kBoundSigmaMargin = 1.5;
    constexpr double kBoundSlackCm = 2.0;
    constexpr int kNewtonIters = 8;
    constexpr double kOnPlaneCm = 1e-4;
    constexpr int kRefineSteps = 2;
    constexpr std::array<double, 5> kJacEps{{1e-6, 1e-5, 1e-5, 1e-3, 1e-3}};
    constexpr double kNeutral = 1.0;  // ORP.NEUTRAL: a missing ranking dim's expected chi2
    constexpr double kPredAngleGuard = 12.0;  // ORP.PRED_ANGLE_MAX_ABS (the ellipse update guard)

    double dot3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
      return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    std::string readFile(const std::string& path) {
      std::ifstream in(path, std::ios::binary);
      if (!in)
        throw cms::Exception("SmartPixelsRefitV3") << "cannot open payload " << path;
      std::ostringstream ss;
      ss << in.rdbuf();
      return ss.str();
    }

    // A 4x4 solve with partial pivoting (numpy.linalg.solve's LU); returns NaN on a singular matrix.
    Vec4 solve4(Mat4 M, Vec4 b) {
      for (int c = 0; c < 4; ++c) {
        int piv = c;
        for (int r = c + 1; r < 4; ++r)
          if (std::abs(M[r][c]) > std::abs(M[piv][c]))
            piv = r;
        if (!(std::abs(M[piv][c]) > 0.)) {
          return Vec4{{kNaN, kNaN, kNaN, kNaN}};
        }
        std::swap(M[c], M[piv]);
        std::swap(b[c], b[piv]);
        for (int r = c + 1; r < 4; ++r) {
          const double f = M[r][c] / M[c][c];
          for (int k = c; k < 4; ++k)
            M[r][k] -= f * M[c][k];
          b[r] -= f * b[c];
        }
      }
      Vec4 x{};
      for (int r = 3; r >= 0; --r) {
        double s = b[r];
        for (int k = r + 1; k < 4; ++k)
          s -= M[r][k] * x[k];
        x[r] = s / M[r][r];
      }
      return x;
    }

    Mat4 jcjt(const Jac& J, const Mat5& C) {
      Mat4 S{};
      for (int i = 0; i < 4; ++i)
        for (int l = 0; l < 4; ++l) {
          double s = 0.;
          for (int j = 0; j < 5; ++j)
            for (int k = 0; k < 5; ++k)
              s += J[i][j] * C[j][k] * J[l][k];
          S[i][l] = s;
        }
      return S;
    }

    // engine.guard_zeroed: a column with any non-finite or |J| > jmax entry is zeroed
    Jac guardZeroed(const Jac& J, double jmax) {
      Jac out = J;
      for (int j = 0; j < 5; ++j) {
        bool bad = false;
        for (int r = 0; r < 4; ++r)
          bad |= !std::isfinite(J[r][j]) || std::abs(J[r][j]) > jmax;
        for (int r = 0; r < 4; ++r)
          out[r][j] = bad ? 0. : J[r][j];
      }
      return out;
    }

    double ptOf(const Vec5& a, double kpt) { return kpt / std::max(std::abs(a[0]), 1e-12); }

    void addKink(Mat5& C, double th2, double rs, double tanl) {
      const double sec2 = 1.0 + tanl * tanl, sec = std::sqrt(sec2);
      const Vec5 jT{{0., sec, 0., 0., +rs * sec}};
      const Vec5 jL{{0., 0., sec2, -rs * sec2, 0.}};
      for (int i = 0; i < 5; ++i)
        for (int k = 0; k < 5; ++k)
          C[i][k] += th2 * (jT[i] * jT[k] + jL[i] * jL[k]);
    }
  }  // namespace

  // ---------------------------------------------------------------- geometry
  void Geometry::finalize() {
    std::sort(modules_.begin(), modules_.end(), [](const Module& a, const Module& b) { return a.detId < b.detId; });
    byId_.clear();
    bySurface_.clear();
    bySide_.clear();
    R_.clear();
    halfZ_.clear();
    Z_.clear();
    rext_.clear();
    surfaces_.clear();
    for (int i = 0; i < static_cast<int>(modules_.size()); ++i) {
      byId_[modules_[i].detId] = i;
      bySurface_[modules_[i].surface].push_back(i);
    }
    for (const int S : {1, 2, 3, 4, 11, 12, 13, 14}) {
      if (!bySurface_.count(S))
        continue;
      surfaces_.push_back(S);
      const auto& rows = bySurface_[S];
      double sumR = 0., sumZ = 0., maxZ = 0., rmin = kInf, rmax = 0.;
      for (const int i : rows) {
        const Module& m = modules_[i];
        sumR += std::hypot(m.o[0], m.o[1]);
        sumZ += std::abs(m.o[2]);
        for (const double fx : {-1., 0., 1.})
          for (const double fy : {-1., 0., 1.}) {
            if (fx == 0. && fy == 0.)
              continue;
            std::array<double, 3> p;
            for (int k = 0; k < 3; ++k)
              p[k] = m.o[k] + (fx * m.hw) * m.ex[k] + (fy * m.hl) * m.ey[k];
            maxZ = std::max(maxZ, std::abs(p[2]));
            const double r = std::hypot(p[0], p[1]);
            rmin = std::min(rmin, r);
            rmax = std::max(rmax, r);
          }
        if (isDisc(S))
          bySide_[S * 2 + (m.o[2] > 0 ? 1 : 0)].push_back(i);
      }
      if (isDisc(S)) {
        Z_[S] = sumZ / rows.size();
        rext_[S] = {rmin, rmax};
      } else {
        R_[S] = sumR / rows.size();
        halfZ_[S] = maxZ;
      }
    }
    R_[0] = 0.;
    ready_ = true;
  }

  int Geometry::index(uint32_t detId) const {
    const auto it = byId_.find(detId);
    return it == byId_.end() ? -1 : it->second;
  }

  // ---------------------------------------------------------------- helix
  HelixPoint helixAt(const Vec5& a, double s) {
    const double rinv = a[0], phi0 = a[1], tanl = a[2], z0 = a[3], d0 = a[4];
    const double x0 = d0 * std::sin(phi0), y0 = -d0 * std::cos(phi0);
    HelixPoint h;
    if (std::abs(rinv) < 1e-6) {
      h.p = {{x0 + s * std::cos(phi0), y0 + s * std::sin(phi0), z0 + tanl * s}};
      h.d = {{std::cos(phi0), std::sin(phi0), tanl}};
    } else {
      const double ph = phi0 - rinv * s;
      h.p = {{x0 - (std::sin(ph) - std::sin(phi0)) / rinv, y0 + (std::cos(ph) - std::cos(phi0)) / rinv, z0 + tanl * s}};
      h.d = {{std::cos(ph), std::sin(ph), tanl}};
    }
    return h;
  }

  double cylinderS(const Vec5& a, double R) {
    const double rinv = a[0], phi0 = a[1], d0 = a[4];
    const double x0 = d0 * std::sin(phi0), y0 = -d0 * std::cos(phi0);
    double s;
    if (std::abs(rinv) < 1e-6) {
      const double b = x0 * std::cos(phi0) + y0 * std::sin(phi0);
      const double c = x0 * x0 + y0 * y0 - R * R;
      s = -b + std::sqrt(b * b - c);
    } else {
      const double Rc = 1.0 / rinv;
      const double cx = x0 + Rc * std::sin(phi0), cy = y0 - Rc * std::cos(phi0);
      const double dc = std::hypot(cx, cy), aR = std::abs(Rc);
      const double cosarg = (aR * aR + dc * dc - R * R) / (2.0 * aR * dc);
      const bool ok = (R <= dc + aR) && (R >= std::abs(dc - aR)) && (std::abs(cosarg) <= 1.0);
      s = ok ? std::acos(std::clamp(cosarg, -1.0, 1.0)) * aR : kNaN;
    }
    return s > 0 ? s : kNaN;
  }

  double discS(const Vec5& a, double Z) {
    const double rinv = a[0], tanl = a[2], z0 = a[3];
    if (!(std::abs(tanl) >= 1e-9))
      return kNaN;
    const double s = ((tanl > 0 ? 1.0 : -1.0) * Z - z0) / tanl;
    if (!(s > 0) || std::abs(rinv) * s > M_PI)
      return kNaN;
    return s;
  }

  double surfaceS(const Vec5& a, const Geometry& g, int S) { return isDisc(S) ? discS(a, g.Z(S)) : cylinderS(a, g.R(S)); }

  Projection projectToModule(const Vec5& a, double s0, const Module& m, bool refine) {
    Projection out;
    bool ok = std::isfinite(s0);
    double s = ok ? s0 : 1.0;
    const double sLim = 3.0 * (std::isfinite(s0) ? s0 : 1.0) + 10.0;
    HelixPoint h;
    for (int it = 0; it < kNewtonIters; ++it) {
      h = helixAt(a, s);
      std::array<double, 3> rel{{h.p[0] - m.o[0], h.p[1] - m.o[1], h.p[2] - m.o[2]}};
      const double f = dot3(rel, m.ez), fp = dot3(h.d, m.ez);
      const bool conv = std::abs(f) < kOnPlaneCm;
      const bool stepOk = std::abs(fp) > 1e-9;
      if (!(conv || !stepOk))
        s = s - f / fp;
      ok = ok && (stepOk || conv);
      ok = ok && (s > 0) && (s < sLim);
    }
    h = helixAt(a, s);
    {
      std::array<double, 3> rel{{h.p[0] - m.o[0], h.p[1] - m.o[1], h.p[2] - m.o[2]}};
      ok = ok && std::abs(dot3(rel, m.ez)) < kOnPlaneCm;
    }
    if (refine) {
      for (int it = 0; it < kRefineSteps; ++it) {
        std::array<double, 3> rel{{h.p[0] - m.o[0], h.p[1] - m.o[1], h.p[2] - m.o[2]}};
        const double f = dot3(rel, m.ez), fp = dot3(h.d, m.ez);
        if (ok && std::abs(fp) > 1e-9)
          s = s - f / fp;
        h = helixAt(a, s);
      }
    }
    std::array<double, 3> rel{{h.p[0] - m.o[0], h.p[1] - m.o[1], h.p[2] - m.o[2]}};
    double lz = dot3(h.d, m.ez);
    lz = std::abs(lz) > 1e-9 ? lz : 1e-9;
    out.u = dot3(rel, m.ex);
    out.v = dot3(rel, m.ey);
    out.cotA = dot3(h.d, m.ex) / lz;
    out.cotB = dot3(h.d, m.ey) / lz;
    out.s = s;
    out.ok = ok;
    return out;
  }

  ProjectionCov projectWithCov(const Vec5& a, const Mat5& C, double s0, const Module& m, bool oneSided) {
    ProjectionCov out;
    out.p = projectToModule(a, s0, m, false);
    const bool disc = isDisc(m.surface);
    const std::array<double, 4> nom{{out.p.u, out.p.v, out.p.cotA, out.p.cotB}};
    for (int j = 0; j < 5; ++j) {
      Vec5 ap = a, am = a;
      ap[j] += kJacEps[j];
      am[j] -= kJacEps[j];
      const Projection pp = projectToModule(ap, out.p.s, m, disc);
      const Projection pm = projectToModule(am, out.p.s, m, disc);
      const std::array<double, 4> vp{{pp.u, pp.v, pp.cotA, pp.cotB}}, vm{{pm.u, pm.v, pm.cotA, pm.cotB}};
      const bool both = pp.ok && pm.ok && !oneSided;
      for (int r = 0; r < 4; ++r) {
        const double cen = (vp[r] - vm[r]) / (2 * kJacEps[j]);
        const double fwd = (vp[r] - nom[r]) / kJacEps[j];
        const double bwd = (nom[r] - vm[r]) / kJacEps[j];
        out.J[r][j] = both ? cen : (pp.ok ? fwd : (pm.ok ? bwd : kNaN));
      }
    }
    out.S = jcjt(out.J, C);
    bool fin = true;
    for (const auto& row : out.S)
      for (const double x : row)
        fin = fin && std::isfinite(x);
    out.p.ok = out.p.ok && fin;
    return out;
  }

  // ---------------------------------------------------------------- visit plan
  VisitPlan visitPlan(const Vec5& a, const Geometry& g, const std::array<bool, 23>& active, bool outsideIn,
                      double margin) {
    struct E {
      double key;
      int surface;
      bool inpath;
    };
    std::vector<E> es;
    const double ad0 = std::abs(a[4]);
    VisitPlan vp;
    for (const int S : g.surfaces()) {
      const double s = surfaceS(a, g, S);
      const bool ok = std::isfinite(s);
      const HelixPoint h = helixAt(a, ok ? s : 0.0);
      if (isDisc(S)) {
        const double r = std::hypot(h.p[0], h.p[1]);
        const auto [lo, hi] = g.rext(S);
        const bool inpath = ok && r >= lo - margin && r <= hi + margin;
        es.push_back({inpath ? s : kInf, S, inpath});
        vp.discGoing = vp.discGoing || inpath;
      } else {
        es.push_back({ok ? s : ((g.R(S) < ad0 ? -1e9 : 1e9) + S), S, true});
      }
    }
    std::stable_sort(es.begin(), es.end(), [](const E& x, const E& y) { return x.key < y.key; });
    if (outsideIn)
      std::reverse(es.begin(), es.end());
    for (const E& e : es)
      if (e.inpath && active[e.surface])
        vp.seq.push_back(e.surface);
    return vp;
  }

  // ---------------------------------------------------------------- payloads
  void SeedModel::load(const std::string& path) {
    const auto j = nlohmann::json::parse(readFile(path));
    if (j.at("segments").at("kind").get<std::string>() != "eta")
      throw cms::Exception("SmartPixelsRefitV3") << "seed model: only the |eta| keying is implemented";
    edges_ = j.at("segments").at("edges").get<std::vector<double>>();
    const auto sc = j.at("scale").get<std::vector<double>>();
    if (sc.size() != 5 || j.at("trans").get<std::vector<int>>() != std::vector<int>{0, 1, 4} ||
        j.at("polar").get<std::vector<int>>() != std::vector<int>{2, 3})
      throw cms::Exception("SmartPixelsRefitV3") << "seed model: unexpected parameter layout";
    std::copy(sc.begin(), sc.end(), scale_.begin());
    theta_.clear();
    for (const auto& t : j.at("theta")) {
      const auto v = t.get<std::vector<double>>();
      if (v.size() != 36)
        throw cms::Exception("SmartPixelsRefitV3") << "seed model: 36 parameters per segment expected";
      std::array<double, 36> a;
      std::copy(v.begin(), v.end(), a.begin());
      theta_.push_back(a);
    }
    if (theta_.size() + 1 != edges_.size())
      throw cms::Exception("SmartPixelsRefitV3") << "seed model: segments and parameters disagree";
    const auto& d = j.at("domain");
    etaCell_ = d.at("eta_cell").get<double>();
    etaLast_ = d.at("eta_last").get<double>();
    pLo_ = d.at("p_lo").get<std::vector<double>>();
    pHi_ = d.at("p_hi").get<std::vector<double>>();
  }

  Mat5 SeedModel::dC(double aeta, double p) const {
    // segment index and x (seedmodel.segment_index / segment_x)
    const int nseg = static_cast<int>(edges_.size()) - 1;
    int sidx = static_cast<int>(std::upper_bound(edges_.begin(), edges_.end(), aeta) - edges_.begin()) - 1;
    sidx = std::clamp(sidx, 0, nseg - 1);
    const double x = aeta - 0.5 * (edges_[sidx] + std::min(edges_[sidx + 1], 2.6));
    // p_eff (seedmodel.p_eff)
    int cell = static_cast<int>(std::min(aeta, etaLast_ - 1e-9) / etaCell_);
    cell = std::min(cell, static_cast<int>(pLo_.size()) - 1);
    const double pe = std::clamp(p, pLo_[cell], pHi_[cell]);
    // entries e = theta[:18] + theta[18:] x: A (9) then B (9); each [T: 3 log-sigma, 3 pcorr; L: 2, 1]
    const auto& th = theta_[sidx];
    std::array<double, 18> e;
    for (int i = 0; i < 18; ++i)
      e[i] = th[i] + th[18 + i] * x;
    auto block = [](const double* ent, int nl, int nz, double* out) {
      double s[3];
      for (int i = 0; i < nl; ++i)
        s[i] = std::exp(std::clamp(ent[i], -30.0, 30.0));
      double z[3];
      for (int i = 0; i < nz; ++i)
        z[i] = 0.999 * std::tanh(ent[nl + i]);
      double R[3][3];
      if (nl == 2) {
        R[0][0] = R[1][1] = 1.;
        R[0][1] = R[1][0] = z[0];
      } else {
        const double z12 = z[0], z13 = z[1], z23 = z[2];
        const double L[3][3] = {{1., 0., 0.},
                                {z12, std::sqrt(1 - z12 * z12), 0.},
                                {z13, z23 * std::sqrt(1 - z13 * z13), std::sqrt((1 - z13 * z13) * (1 - z23 * z23))}};
        for (int i = 0; i < 3; ++i)
          for (int k = 0; k < 3; ++k) {
            double v = 0.;
            for (int l = 0; l < 3; ++l)
              v += L[i][l] * L[k][l];
            R[i][k] = v;
          }
      }
      for (int i = 0; i < nl; ++i)
        for (int k = 0; k < nl; ++k)
          out[i * nl + k] = s[i] * R[i][k] * s[k];
    };
    double AT[9], AL[4], BT[9], BL[4];
    block(&e[0], 3, 3, AT);
    block(&e[6], 2, 1, AL);
    block(&e[9], 3, 3, BT);
    block(&e[15], 2, 1, BL);
    const double ip2 = 1.0 / (pe * pe);
    Mat5 M{};
    const int T[3] = {0, 1, 4}, L2[2] = {2, 3};
    for (int i = 0; i < 3; ++i)
      for (int k = 0; k < 3; ++k)
        M[T[i]][T[k]] = AT[i * 3 + k] + BT[i * 3 + k] * ip2;
    for (int i = 0; i < 2; ++i)
      for (int k = 0; k < 2; ++k)
        M[L2[i]][L2[k]] = AL[i * 2 + k] + BL[i * 2 + k] * ip2;
    for (int i = 0; i < 5; ++i)
      for (int k = 0; k < 5; ++k)
        M[i][k] = M[i][k] * scale_[i] * scale_[k];
    return M;
  }

  void Material::load(const std::string& path) {
    const auto j = nlohmann::json::parse(readFile(path));
    if (j.at("schema").get<std::string>() != "spix-refit-v3-material/1")
      throw cms::Exception("SmartPixelsRefitV3") << "material: unknown schema";
    c0 = j.at("highland").at("c0").get<double>();
    c1 = j.at("highland").at("c1").get<double>();
    barrel.clear();
    for (const auto& p : j.at("barrel_planes")) {
      MaterialPlane m;
      m.layer = p.at("layer").get<int>();
      m.r = m.layer ? -1. : p.at("r").get<double>();
      m.xx0 = p.at("xx0").get<double>();
      barrel.push_back(m);
    }
    discGoing.clear();
    for (const auto& p : j.at("disc_going_planes")) {
      MaterialPlane m;
      m.disc = p.at("kind").get<std::string>() == "disc";
      if (m.disc) {
        m.z = p.at("z").get<double>();
        m.rmin = p.at("r_ext").at(0).get<double>();
        m.rmax = p.at("r_ext").at(1).get<double>();
      } else {
        m.r = p.at("r").get<double>();
        m.halfZ = p.at("half_z").get<double>();
      }
      m.xx0 = p.at("xx0").get<double>();
      discGoing.push_back(m);
    }
    il1Module = MaterialPlane{};
    il1Module.layer = j.at("il1_module").at("layer").get<int>();
    il1Module.r = -1.;
    il1Module.xx0 = j.at("il1_module").at("xx0").get<double>();
    beamPipe = MaterialPlane{};
    beamPipe.r = j.at("beam_pipe").at("r").get<double>();
    beamPipe.xx0 = j.at("beam_pipe").at("xx0").get<double>();
    beamPipe.halfZ = kInf;
    materialScale = j.at("material_scale").get<double>();
    cylScale = j.at("cylinder_material_scale").get<double>();
    discScale = j.at("disc_material_scale").get<double>();
    otInnerR = j.at("ot_inner_r_cm").get<double>();
  }

  void Payloads::load(const std::string& manifestPath) {
    const std::string text = readFile(manifestPath);
    const auto man = nlohmann::json::parse(text);
    if (man.at("schema").get<std::string>() != "spix-refit-v3-payload/1")
      throw cms::Exception("SmartPixelsRefitV3") << manifestPath << ": unknown payload manifest schema";
    setId_ = man.at("set_id").get<std::string>();
    const std::filesystem::path dir = std::filesystem::path(manifestPath).parent_path();
    auto file = [&](const std::string& name) {
      const auto& e = man.at("files").at(name);
      const std::string path = (dir / name).string();
      const std::string want = e.at("md5").get<std::string>();
      const std::string got = cms::Digest(readFile(path)).digest().toString();
      if (got != want)
        throw cms::Exception("SmartPixelsRefitV3") << path << ": md5 " << got << " differs from the manifest's " << want;
      return path;
    };
    sets_.clear();
    auto cs = [&](const std::string& name) -> const correction::CorrectionSet& {
      sets_.emplace_back(correction::CorrectionSet::from_file(file(name)));
      return *sets_.back();
    };
    const auto& r = cs("rtables.json");
    rU_ = r.at("R_u");
    rV_ = r.at("R_v");
    const auto& a = cs("angles.json");
    orient_ = a.at("module_orientation");
    angBarrel_ = a.at("angle_barrel");
    angDisc_ = a.at("angle_disc");
    skip_ = cs("skip_price.json").at("skip_price");
    seed.load(file("seed_model.json"));
    material.load(file("material.json"));
    loaded_ = true;
  }

  double Payloads::rU(int sizeX, double c) const { return rU_->evaluate({static_cast<double>(sizeX), c}); }
  double Payloads::rV(int sizeY, double c) const { return rV_->evaluate({static_cast<double>(sizeY), c}); }
  int Payloads::orientation(uint32_t detId) const {
    return static_cast<int>(orient_->evaluate({static_cast<int>(detId)}));
  }

  Payloads::AngleEntry Payloads::barrelAngle(int surface, int sign, double key, bool beta) const {
    AngleEntry e;
    if (sign == 0 || !std::isfinite(key))
      return e;
    const std::string a = beta ? "beta" : "alpha";
    e.mode = static_cast<int>(angBarrel_->evaluate({surface, sign, key, "mode_" + a}));
    e.c = angBarrel_->evaluate({surface, sign, key, "c_" + a});
    e.w = angBarrel_->evaluate({surface, sign, key, "w_" + a});
    return e;
  }

  Payloads::AngleEntry Payloads::discAngle(int surface, double key, bool beta) const {
    AngleEntry e;
    if (!std::isfinite(key))
      return e;
    const std::string a = beta ? "beta" : "alpha";
    e.mode = static_cast<int>(angDisc_->evaluate({surface, key, "mode_" + a}));
    e.c = angDisc_->evaluate({surface, key, "c_" + a});
    e.w = angDisc_->evaluate({surface, key, "w_" + a});
    return e;
  }

  double Payloads::skipPrice(const std::string& dims, int code, bool first) const {
    return skip_->evaluate({dims, code, first ? 1 : 0});
  }

  // ---------------------------------------------------------------- the engine
  struct Engine::Copy {
    Vec5 a{};
    Mat5 C{};
    int nUpd = 0;
    int Lv = kOTSurface;
    int branch = 0;
    bool pend = false;
    double s1 = 0., P1 = 0.;
    int j1 = 0;
    std::vector<StepRecord> rec;
  };

  struct Engine::Cand {
    int module = -1;   // geometry index
    int cl = -1;       // position in the module's cluster list
    Vec4 h0{};
    Jac J{};
    Mat4 Srank{};      // J C J^T (unguarded): ranking
    Mat4 Sz{};         // the guarded one: the consumers' statistics
    double val = kInf;
  };

  double Engine::highlandTheta2(double xx0, double p, double tanl, double c0, double c1) {
    const double X = std::max(xx0 * std::sqrt(1.0 + tanl * tanl), 1e-12);
    const double th = c0 / p * std::sqrt(X) * (1.0 + c1 * std::log(X));
    return th > 0 ? th * th : 0.;
  }

  double Engine::innovationChi2(const Vec4& r, const Mat4& S, const Vec4& sig, const std::array<bool, 4>& dims) {
    Mat4 M = S;
    for (int i = 0; i < 4; ++i)
      M[i][i] += dims[i] ? sig[i] * sig[i] : 1.0;
    for (int i = 0; i < 4; ++i)
      for (int k = 0; k < 4; ++k)
        if (!dims[i] || !dims[k])
          M[i][k] = 0.;
    for (int i = 0; i < 4; ++i)
      if (!dims[i])
        M[i][i] = 1.;
    Vec4 rz;
    for (int i = 0; i < 4; ++i)
      rz[i] = dims[i] ? (std::isfinite(r[i]) ? r[i] : 0.) : 0.;
    const Vec4 x = solve4(M, rz);
    double c = 0.;
    for (int i = 0; i < 4; ++i)
      c += rz[i] * x[i];
    return c;
  }

  std::array<double, 3> Engine::tree2Scores(const std::array<double, 2>& x1,
                                            const std::array<std::array<double, 2>, 3>* x2, double P) {
    std::array<double, 3> s1{{x1[0], x1[1], P}};
    if (!x2)
      return s1;
    std::array<double, 3> out;
    for (int c = 0; c < 3; ++c) {
      const double cand = std::min((*x2)[c][0], (*x2)[c][1]);
      const bool empty = !std::isfinite(cand);
      double L2;
      if (c < 2)
        L2 = empty ? P : std::min(cand, P);
      else
        L2 = empty ? P : cand;
      out[c] = s1[c] + L2;
    }
    return out;
  }

  bool Engine::uses(int consumer, int angle, int surface, int step) const {
    if ((cfg_.angleSteps == 1 && step > 0) || (cfg_.angleSteps == 2 && step == 0))
      return false;
    return cfg_.use[consumer][angle][isDisc(surface) ? 1 : 0];
  }

  std::string Engine::dimsKey(int step) const {
    // arms.Arm.dims_key("rank"): the rank consumer's angle dims on the BARREL, at this step
    std::string k = "xy";
    if (uses(0, 0, 1, step))
      k += "a";
    if (uses(0, 1, 1, step))
      k += "b";
    if (cfg_.knownBad == "skip_price_xyab_always")
      k = "xyab";
    return k;
  }

  void Engine::addMaterialBarrel(Vec5& a, Mat5& C, int S, int Lv, double kpt, bool inner) const {
    const Material& mat = pay_.material;
    const double tanl = a[2];
    const double p = ptOf(a, kpt) * std::sqrt(1.0 + tanl * tanl);
    if (!inner && isDisc(S))
      return;  // a disc visit of a track treated as barrel-only (the no-margin known-bad): no radius rule
    const double rPrev = (Lv >= 1 && Lv <= 4) ? geo_.R(Lv) : mat.otInnerR;
    const double rS = inner ? 0. : geo_.R(S);
    auto charge = [&](const MaterialPlane& pl) {
      const double rs = pl.layer ? geo_.R(pl.layer) : pl.r;
      if (rS < rs && rs <= rPrev + 1e-6)
        addKink(C, highlandTheta2(mat.materialScale * pl.xx0, p, tanl, mat.c0, mat.c1), rs, tanl);
    };
    for (const auto& pl : mat.barrel)
      charge(pl);
    if (inner) {
      charge(mat.il1Module);
      charge(mat.beamPipe);
    }
  }

  void Engine::chargePath(Vec5& a, Mat5& C, double sCur, double sPrev, const std::vector<MaterialPlane>& planes,
                          double kpt) const {
    const Material& mat = pay_.material;
    const double tanl = a[2];
    const double p = ptOf(a, kpt) * std::sqrt(1.0 + tanl * tanl);
    const double discScale = cfg_.knownBad == "disc_scale_1" ? 1.0 : mat.discScale;
    // every crossing on the incoming state first (disc_material.crossings), then the kinks
    struct X {
      double s, r;
      bool on;
    };
    std::vector<X> xs;
    xs.reserve(planes.size());
    for (const auto& pl : planes) {
      double sj, rj;
      bool on;
      if (!pl.disc) {
        sj = cylinderS(a, pl.r);
        const HelixPoint h = helixAt(a, std::isfinite(sj) ? sj : 0.0);
        on = std::isfinite(sj) && std::abs(h.p[2]) <= pl.halfZ;
        rj = pl.r;
      } else {
        sj = discS(a, pl.z);
        const HelixPoint h = helixAt(a, std::isfinite(sj) ? sj : 0.0);
        rj = std::hypot(h.p[0], h.p[1]);
        on = std::isfinite(sj) && rj >= pl.rmin && rj <= pl.rmax;
      }
      xs.push_back({sj, rj, on});
    }
    // (sCur, sPrev]: the previous surface's own plane in, the current surface's out. The replay excludes the
    // current plane by exact equality (its planes and surfaces share one geometry); here the planes come from the
    // payload and the surfaces from the reconstruction geometry, so both ends take the same 1e-6 cm tolerance.
    for (size_t k = 0; k < planes.size(); ++k) {
      const auto& pl = planes[k];
      if (!(xs[k].on && xs[k].s > sCur + 1e-6 && xs[k].s <= sPrev + 1e-6))
        continue;
      const double xx0 = pl.disc ? discScale * pl.xx0 : mat.cylScale * pl.xx0;
      const double path = pl.disc ? 1.0 / std::max(std::abs(tanl), 1e-6) : 1.0;
      addKink(C, highlandTheta2(xx0 * path, p, tanl, mat.c0, mat.c1), xs[k].r, tanl);
    }
  }

  void Engine::qStep(Copy& c, int S, bool discGoing, double kpt) const {
    if (!discGoing) {
      addMaterialBarrel(c.a, c.C, S, c.Lv, kpt, false);
    } else {
      const double sCur = surfaceS(c.a, geo_, S);
      const double sPrev = c.Lv == kOTSurface ? kInf : surfaceS(c.a, geo_, c.Lv);
      chargePath(c.a, c.C, sCur, sPrev, pay_.material.discGoing, kpt);
    }
    c.Lv = S;
  }

  void Engine::innerMaterial(Copy& c, bool discGoing, double kpt) const {
    int Lv = -1;
    for (const auto& r : c.rec)
      if (r.accepted)
        Lv = r.surface;
    if (Lv < 0)
      return;
    const bool byRadius = cfg_.knownBad == "inner_material_by_radius";
    if (!discGoing || byRadius) {
      addMaterialBarrel(c.a, c.C, 0, Lv, kpt, true);
    } else {
      std::vector<MaterialPlane> planes = pay_.material.discGoing;
      planes.push_back(pay_.material.beamPipe);
      chargePath(c.a, c.C, 0., surfaceS(c.a, geo_, Lv), planes, kpt);
    }
  }

  void Engine::measurements(const Cand& k, const Cluster& cl, Vec4& meas, Vec4& sig, std::array<bool, 2>& has) const {
    const Module& m = geo_.modules()[k.module];
    const int S = m.surface;
    meas = {{cl.x, cl.y, cl.cotA, cl.cotB}};
    sig = {{cl.sigX, cl.sigY, cl.sigA, cl.sigB}};
    has = {{cl.hasA, cl.hasB}};
    const double ca = k.h0[2], cb = k.h0[3];
    if (cfg_.rtables) {
      const bool swap = cfg_.knownBad == "r_lookup_swapped";
      sig[0] = pay_.rU(cl.sizeX, std::abs(swap ? cb : ca));
      sig[1] = pay_.rV(cl.sizeY, std::abs(swap ? ca : cb));
    }
    if (!cfg_.angleTable)
      return;
    for (int ia = 0; ia < 2; ++ia) {
      Payloads::AngleEntry e;
      if (!isDisc(S)) {
        const int s = pay_.orientation(m.detId);
        e = pay_.barrelAngle(S, s, barrelKey(s, ca, cb), ia == 1);
        if (cfg_.knownBad == "barrel_centre_global_sign" && ia == 1 && s != 0)
          e.c = s * e.c;  // the centre without the module's orientation sign
      } else {
        double key = discKey(ca, cb);
        if (cfg_.knownBad == "disc_key_facing_mirrored" && m.ez[2] < 0)
          key = -key;
        e = pay_.discAngle(S, key, ia == 1);
      }
      applyAngleEntry(e, meas[2 + ia], sig[2 + ia]);
    }
  }

  double Engine::gateStat(const Cand& k, const Cluster& cl, int step, int consumer) const {
    Vec4 meas, sig;
    std::array<bool, 2> has;
    measurements(k, cl, meas, sig, has);
    const int S = geo_.modules()[k.module].surface;
    const std::array<bool, 4> dims{{true, true, has[0] && uses(consumer, 0, S, step), has[1] && uses(consumer, 1, S, step)}};
    Vec4 r;
    for (int i = 0; i < 4; ++i)
      r[i] = meas[i] - k.h0[i];
    return innovationChi2(r, k.Sz, sig, dims);
  }

  std::vector<Engine::Cand> Engine::candidates(const Copy& c, int S, int step, const EventClusters& clusters,
                                               bool& fallback) const {
    std::vector<Cand> out;
    fallback = false;
    const double nsig = cfg_.knownBad == "window_3sigma" ? 3.0 : cfg_.nsig;
    const bool oneSided = cfg_.knownBad == "one_sided_jacobian";
    const double s0 = surfaceS(c.a, geo_, S);
    if (!std::isfinite(s0))
      return out;
    const HelixPoint pc = helixAt(c.a, s0);
    const auto& rows = geo_.bySurface(S);
    const bool disc = isDisc(S);
    const bool up = c.a[2] > 0;
    std::vector<double> dist(rows.size(), kInf);
    int m0 = -1;
    double d0min = kInf;
    for (size_t i = 0; i < rows.size(); ++i) {
      const Module& m = geo_.modules()[rows[i]];
      if (disc && ((m.o[2] > 0) != up))
        continue;
      dist[i] = std::sqrt((pc.p[0] - m.o[0]) * (pc.p[0] - m.o[0]) + (pc.p[1] - m.o[1]) * (pc.p[1] - m.o[1]) +
                          (pc.p[2] - m.o[2]) * (pc.p[2] - m.o[2]));
      if (dist[i] < d0min) {
        d0min = dist[i];
        m0 = rows[i];
      }
    }
    if (m0 < 0)
      return out;
    const ProjectionCov p0 = projectWithCov(c.a, c.C, s0, geo_.modules()[m0], oneSided);
    double rEll = kInf;  // a failed nearest-module projection: every module of the surface (side)
    fallback = !p0.p.ok;
    if (p0.p.ok) {
      const double tr = (std::isfinite(p0.S[0][0]) ? p0.S[0][0] : 0.) + (std::isfinite(p0.S[1][1]) ? p0.S[1][1] : 0.);
      rEll = kBoundSigmaMargin * nsig * std::sqrt(std::max(tr, 0.));
    }
    const double tl2 = c.a[2] * c.a[2];
    const double slack = kBoundSlackCm * std::sqrt(1.0 + (disc ? 1.0 / std::max(tl2, 1e-12) : tl2));
    for (size_t i = 0; i < rows.size(); ++i) {
      const Module& m = geo_.modules()[rows[i]];
      if (!(dist[i] <= rEll + slack + std::hypot(m.hw, m.hl)))
        continue;
      const ProjectionCov pj = projectWithCov(c.a, c.C, s0, m, oneSided);
      if (!pj.p.ok)
        continue;
      const double su = std::sqrt(std::max(pj.S[0][0], 0.)), sv = std::sqrt(std::max(pj.S[1][1], 0.));
      if (!(std::abs(pj.p.u) <= m.hw + nsig * su && std::abs(pj.p.v) <= m.hl + nsig * sv))
        continue;
      const auto& cls = clusters[rows[i]];
      const double det = pj.S[0][0] * pj.S[1][1] - pj.S[0][1] * pj.S[0][1];
      const Jac Jz = guardZeroed(pj.J, cfg_.jacobianMaxAbs);
      const Mat4 Sz = jcjt(Jz, c.C);
      for (size_t q = 0; q < cls.size(); ++q) {
        const double du = cls[q].x - pj.p.u, dv = cls[q].y - pj.p.v;
        const double m2 =
            det > 0 ? (pj.S[1][1] * du * du - 2 * pj.S[0][1] * du * dv + pj.S[0][0] * dv * dv) / det : kInf;
        if (!(m2 < nsig * nsig))
          continue;
        Cand k;
        k.module = rows[i];
        k.cl = static_cast<int>(q);
        k.h0 = {{pj.p.u, pj.p.v, pj.p.cotA, pj.p.cotB}};
        k.J = pj.J;
        k.Srank = pj.S;
        k.Sz = Sz;
        // the ranking value (arms.rank_values)
        Vec4 meas, sig;
        std::array<bool, 2> has;
        measurements(k, cls[q], meas, sig, has);
        Vec4 r;
        for (int d = 0; d < 4; ++d)
          r[d] = meas[d] - k.h0[d];
        if (!has[0])
          r[2] = kNaN;
        if (!has[1])
          r[3] = kNaN;
        const std::array<bool, 4> rdims{{true, true, has[0] && uses(0, 0, S, step), has[1] && uses(0, 1, S, step)}};
        if (cfg_.innovationMetric) {
          std::array<bool, 4> have;
          int nMiss = 0;
          for (int d = 0; d < 4; ++d) {
            have[d] = std::isfinite(r[d]) && sig[d] > 0 && rdims[d];
            nMiss += !have[d];
          }
          k.val = innovationChi2(r, k.Srank, sig, have) + kNeutral * nMiss;
        } else {
          double v = 0.;
          for (int d = 0; d < 4; ++d)
            if (rdims[d]) {
              const double z = (std::isfinite(r[d]) ? r[d] : 0.) / sig[d];
              v += z * z;
            }
          k.val = v;
        }
        if (!std::isfinite(k.val))
          k.val = kInf;
        out.push_back(k);
      }
    }
    return out;
  }

  void Engine::apply(Copy& c, int step, const Cand& k, const Cluster& cl) const {
    StepRecord& rec = c.rec[step];
    if (!(std::abs(k.h0[2]) <= kPredAngleGuard && std::abs(k.h0[3]) <= kPredAngleGuard))
      return;
    Vec4 meas, sig;
    std::array<bool, 2> has;
    measurements(k, cl, meas, sig, has);
    const int S = geo_.modules()[k.module].surface;
    const std::array<bool, 4> ud{{true, true, has[0] && uses(1, 0, S, step), has[1] && uses(1, 1, S, step)}};
    // engine.kf_update: sequential scalars with relinearisation; guarded Jacobian columns zeroed
    const Jac J = guardZeroed(k.J, cfg_.jacobianMaxAbs);
    const Vec5 aLin = c.a;
    rec.aPre = c.a;
    rec.CPred = c.C;
    for (int s = 0; s < 4; ++s) {
      const auto& H = J[s];
      Vec5 v{};
      for (int i = 0; i < 5; ++i) {
        double x = 0.;
        for (int j = 0; j < 5; ++j)
          x += c.C[i][j] * H[j];
        v[i] = x;
      }
      double Sv = 0.;
      for (int i = 0; i < 5; ++i)
        Sv += H[i] * v[i];
      Sv += sig[s] * sig[s];
      double lin = 0.;
      for (int i = 0; i < 5; ++i)
        lin += H[i] * (c.a[i] - aLin[i]);
      const double r = meas[s] - k.h0[s] - lin;
      bool ok = ud[s] && Sv > 0 && std::isfinite(r);
      const double c2 = ok ? r * r / Sv : 0.;
      ok = ok && !(c2 > cfg_.chi2UpdateGate);
      if (!ok)
        continue;
      rec.pull[s] = r / std::sqrt(Sv);
      rec.chi2[s] = c2;
      rec.applied[s] = true;
      Vec5 K;
      for (int i = 0; i < 5; ++i)
        K[i] = v[i] / Sv;
      for (int i = 0; i < 5; ++i)
        c.a[i] += K[i] * r;
      for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
          c.C[i][j] -= K[i] * v[j];
    }
    rec.aPost = c.a;
    rec.CPost = c.C;
    rec.accepted = true;
    rec.detId = geo_.modules()[k.module].detId;
    rec.cluster = cl.index;
    rec.h0 = k.h0;
    rec.projSig = {{std::sqrt(std::max(k.Srank[0][0], 0.)), std::sqrt(std::max(k.Srank[1][1], 0.))}};
    rec.has = has;
    rec.meas = meas;
    rec.sig = sig;
    ++c.nUpd;
  }

  TrackResult Engine::run(const Vec5& aOT, const Mat5& COT, double kpt, const std::array<bool, 23>& active,
                          const EventClusters& clusters) const {
    TrackResult res;
    res.a = aOT;
    res.C = COT;
    double tr = 0.;
    for (int i = 0; i < 5; ++i)
      tr += COT[i][i];
    res.covOK = tr > 0;
    if (!res.covOK)
      return res;
    const VisitPlan vp = visitPlan(aOT, geo_, active, cfg_.outsideIn, cfg_.discMarginCm);
    res.discGoing = vp.discGoing;
    const bool discGoing = cfg_.knownBad == "disc_flag_no_margin" ? visitPlan(aOT, geo_, active, cfg_.outsideIn, 0.).discGoing
                                                                  : vp.discGoing;
    const int ns = static_cast<int>(vp.seq.size());

    // seed (arms.seed): the OT helix; the d0 floor; + the seed model
    Copy c0;
    c0.a = aOT;
    c0.C = COT;
    if (cfg_.otNPar == 4 || !(c0.C[4][4] > 0))
      c0.C[4][4] = std::max(c0.C[4][4], 0.25);
    if (cfg_.seedSmooth && cfg_.knownBad != "seed_model_off") {
      const double tanl = aOT[2];
      const double p = ptOf(aOT, kpt) * std::sqrt(1.0 + tanl * tanl);
      const Mat5 dC = pay_.seed.dC(std::abs(std::asinh(tanl)), p);
      for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
          c0.C[i][j] += dC[i][j];
    }
    res.seedC = c0.C;
    c0.rec.assign(ns, StepRecord{});
    std::vector<Copy> copies{c0};
    const bool tree2 = cfg_.tree2 && cfg_.knownBad != "k1_not_two_layer";

    // the commit of the pending windows (arms._tree2_commit); x2 per copy, or none (a one-layer window)
    auto commit = [&](const std::vector<int>& idx, const std::vector<std::array<double, 2>>* x2) {
      std::array<double, 2> X1{{kInf, kInf}};
      std::array<std::array<double, 2>, 3> X2;
      for (auto& r : X2)
        r = {{kInf, kInf}};
      double P = 0.;
      for (size_t q = 0; q < idx.size(); ++q) {
        const Copy& c = copies[idx[q]];
        if (c.branch < 2)
          X1[c.branch] = c.s1;
        P = c.P1;
        if (x2)
          X2[c.branch] = (*x2)[q];
      }
      const auto sc = tree2Scores(X1, x2 ? &X2 : nullptr, P);
      const int win = static_cast<int>(std::min_element(sc.begin(), sc.end()) - sc.begin());
      std::vector<bool> keep(idx.size());
      for (size_t q = 0; q < idx.size(); ++q) {
        Copy& c = copies[idx[q]];
        keep[q] = c.branch == win;
        if (keep[q]) {
          c.rec[c.j1].kchoice = win;
          c.rec[c.j1].scores = sc;
        }
        c.pend = false;
      }
      return keep;
    };

    for (int j = 0; j < ns; ++j) {
      const int S = vp.seq[j];
      std::vector<int> t;
      for (int q = 0; q < static_cast<int>(copies.size()); ++q)
        if (copies[q].nUpd < cfg_.maxKFUpdates)
          t.push_back(q);
      if (t.empty())
        continue;
      for (const int q : t) {
        qStep(copies[q], S, discGoing, kpt);
        copies[q].rec[j].visited = true;
        copies[q].rec[j].surface = S;
      }
      // candidates, best and runner-up per visiting copy
      std::vector<std::vector<Cand>> cands(t.size());
      std::vector<int> kb(t.size(), -1), kr(t.size(), -1);
      std::vector<double> xb(t.size(), kInf), xr(t.size(), kInf);
      for (size_t q = 0; q < t.size(); ++q) {
        cands[q] = candidates(copies[t[q]], S, j, clusters, copies[t[q]].rec[j].fallback);
        copies[t[q]].rec[j].ncand = static_cast<int>(cands[q].size());
        std::vector<int> o(cands[q].size());
        std::iota(o.begin(), o.end(), 0);
        std::stable_sort(o.begin(), o.end(), [&](int x, int y) { return cands[q][x].val < cands[q][y].val; });
        if (!o.empty())
          kb[q] = o[0];
        if (o.size() > 1) {
          kr[q] = o[1];
          copies[t[q]].rec[j].selMargin = cands[q][o[1]].val - cands[q][o[0]].val;
        }
        for (auto [kk, xx] : {std::pair<int, double*>{kb[q], &xb[q]}, std::pair<int, double*>{kr[q], &xr[q]}})
          if (kk >= 0) {
            const Cand& k = cands[q][kk];
            *xx = gateStat(k, clusters[k.module][k.cl], j, 0);
          }
      }
      if (!tree2) {
        for (size_t q = 0; q < t.size(); ++q)
          if (kb[q] >= 0) {
            const Cand& k = cands[q][kb[q]];
            copies[t[q]].rec[j].rankChi2 = xb[q];
            apply(copies[t[q]], j, k, clusters[k.module][k.cl]);
          }
        continue;
      }
      // commit the windows opened at the previous visit
      std::vector<bool> dropped(copies.size(), false);
      bool anyPend = false;
      for (const int q : t)
        anyPend = anyPend || copies[q].pend;
      if (anyPend) {
        std::vector<int> allp;
        std::vector<std::array<double, 2>> x2;
        for (int q = 0; q < static_cast<int>(copies.size()); ++q)
          if (copies[q].pend) {
            allp.push_back(q);
            const auto it = std::find(t.begin(), t.end(), q);
            if (it == t.end())
              x2.push_back({{kInf, kInf}});  // a copy not at this visit: an empty layer
            else {
              const size_t p = it - t.begin();
              x2.push_back({{xb[p], xr[p]}});
            }
          }
        const auto keep = commit(allp, &x2);
        for (size_t q = 0; q < allp.size(); ++q)
          dropped[allp[q]] = !keep[q];
      }
      // open a window on each surviving visiting copy with a candidate
      std::vector<Copy> added;
      for (size_t q = 0; q < t.size(); ++q) {
        if (dropped[t[q]] || kb[q] < 0)
          continue;
        const bool first = j == 0;
        const int code = (discGoing && S <= 4) ? S + 100 : S;
        const double P = cfg_.knownBad == "skip_price_star_only" ? pay_.skipPrice("xyab", 0, false)
                                                                 : pay_.skipPrice(dimsKey(j), code, first);
        Copy& c = copies[t[q]];
        c.rec[j].skipPrice = P;
        Copy ru = c, sk = c;
        c.branch = 0;
        c.pend = true;
        c.s1 = xb[q];
        c.P1 = P;
        c.j1 = j;
        sk.branch = 2;
        sk.pend = true;
        sk.s1 = 0.;
        sk.P1 = P;
        sk.j1 = j;
        if (kr[q] >= 0) {
          ru.branch = 1;
          ru.pend = true;
          ru.s1 = xr[q];
          ru.P1 = P;
          ru.j1 = j;
          const Cand& k = cands[q][kr[q]];
          ru.rec[j].rankChi2 = xr[q];
          apply(ru, j, k, clusters[k.module][k.cl]);
          added.push_back(std::move(ru));
        }
        added.push_back(std::move(sk));
        const Cand& k = cands[q][kb[q]];
        c.rec[j].rankChi2 = xb[q];
        apply(c, j, k, clusters[k.module][k.cl]);
      }
      std::vector<Copy> next;
      for (int q = 0; q < static_cast<int>(copies.size()); ++q)
        if (!dropped[q])
          next.push_back(std::move(copies[q]));
      for (auto& c : added)
        next.push_back(std::move(c));
      copies = std::move(next);
    }
    // windows still open: one layer (arms.finish)
    {
      std::vector<int> p;
      for (int q = 0; q < static_cast<int>(copies.size()); ++q)
        if (copies[q].pend)
          p.push_back(q);
      if (!p.empty()) {
        const auto keep = commit(p, nullptr);
        std::vector<Copy> next;
        for (int q = 0; q < static_cast<int>(copies.size()); ++q) {
          const auto it = std::find(p.begin(), p.end(), q);
          if (it == p.end() || keep[it - p.begin()])
            next.push_back(std::move(copies[q]));
        }
        copies = std::move(next);
      }
    }
    if (copies.size() != 1)
      throw cms::Exception("SmartPixelsRefitV3") << "internal: " << copies.size() << " copies survive a track";
    Copy& c = copies[0];
    if (cfg_.innerMaterial)
      innerMaterial(c, discGoing, kpt);
    res.steps = c.rec;
    res.nacc = 0;
    for (const auto& r : c.rec)
      res.nacc += r.accepted;
    res.refit = res.nacc >= cfg_.minHits;
    if (res.refit) {
      res.a = c.a;
      res.C = c.C;
    }
    return res;
  }

}  // namespace smartpixels::v3
