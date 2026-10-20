// The refit v3 engine against the replay's own numbers (test/refit_v3_unit_reference.json, written by
// ngtagger-train eval_refitq/v3cmssw/export_cmssw.py from the replay's code at the reference commit):
// every table at every index along each axis, the seed model on a grid through its edges, the two-layer
// scores; v3cmssw/PLAN.md rev 1 §3 item 1.
#include "catch2/catch_all.hpp"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitV3.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace smartpixels::v3;

namespace {
  std::string base() {
    const char* b = std::getenv("CMSSW_BASE");
    REQUIRE(b != nullptr);
    return std::string(b) + "/src/L1Trigger/Phase3SmartPixels/";
  }

  const nlohmann::json& reference() {
    static const nlohmann::json j = [] {
      std::ifstream in(base() + "test/refit_v3_unit_reference.json");
      REQUIRE(in.good());
      return nlohmann::json::parse(in);
    }();
    return j;
  }

  const Payloads& payloads() {
    static const Payloads p = [] {
      Payloads q;
      q.load(base() + "data/refit_v3/" + reference().at("set_id").get<std::string>() + "/manifest.json");
      return q;
    }();
    return p;
  }

  // the reference doubles went through two JSON writers and correctionlib's reader: 1 ulp
  bool same(double a, double b, double rel = 4e-16, double abs = 0.) {
    return a == b || std::abs(a - b) <= std::max(rel * std::max(std::abs(a), std::abs(b)), abs);
  }
}  // namespace

TEST_CASE("refit v3 R tables", "[refitV3]") {
  for (const auto& e : reference().at("rtables")) {
    const double size = e.at("size").get<double>(), cot = std::abs(e.at("cotPred").get<double>());
    const bool u = e.at("table").get<std::string>() == "u";
    // the producer passes integer cluster sizes; the reference also probes half-integers and 0
    const double got = u ? payloads().rU(static_cast<int>(size), cot) : payloads().rV(static_cast<int>(size), cot);
    if (size == std::floor(size))
      REQUIRE(same(got, e.at("width").get<double>()));
  }
}

TEST_CASE("refit v3 module orientation", "[refitV3]") {
  for (const auto& e : reference().at("orientation"))
    REQUIRE(payloads().orientation(e.at("detId").get<uint32_t>()) == e.at("sign").get<int>());
}

TEST_CASE("refit v3 angle tables", "[refitV3]") {
  for (const auto& e : reference().at("angles")) {
    const int S = e.at("surface").get<int>();
    const double ca = e.at("cotA").get<double>(), cb = e.at("cotB").get<double>();
    const bool barrel = e.at("region").get<std::string>() == "barrel";
    int sign = 0;
    if (barrel) {
      sign = payloads().orientation(e.at("detId").get<uint32_t>());
      REQUIRE(sign == e.at("sign").get<int>());
    }
    for (int ia = 0; ia < 2; ++ia) {
      const auto& want = e.at(ia ? "beta" : "alpha");
      double m = e.at("m").get<double>(), sig = e.at("sigma").get<double>();
      const auto ent = barrel ? payloads().barrelAngle(S, sign, barrelKey(sign, ca, cb), ia == 1)
                              : payloads().discAngle(S, discKey(ca, cb), ia == 1);
      applyAngleEntry(ent, m, sig);
      if (!want.at("usable").get<bool>()) {
        REQUIRE(std::isnan(m));
        continue;
      }
      INFO("surface " << S << " bin " << e.at("bin") << " angle " << ia);
      REQUIRE(same(m, want.at("m").get<double>(), 1e-15));
      REQUIRE(same(sig, want.at("sigma").get<double>(), 1e-15));
    }
  }
}

TEST_CASE("refit v3 skip price", "[refitV3]") {
  for (const auto& e : reference().at("skip_price")) {
    const auto dims = e.at("dims").get<std::string>();
    const int code = e.at("code").get<int>();
    const bool first = e.at("first").get<int>() == 1;
    REQUIRE(same(payloads().skipPrice(dims, code, first), e.at("value").get<double>()));
    REQUIRE(payloads().skipPriceLevel(dims, code, first) == e.at("level").get<int>());
  }
}

TEST_CASE("refit v3 seed model", "[refitV3]") {
  for (const auto& e : reference().at("seed")) {
    const Mat5 got = payloads().seed.dC(e.at("absEta").get<double>(), e.at("p").get<double>());
    const auto& want = e.at("dC");
    double mx = 0.;
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j)
        mx = std::max(mx, std::abs(want.at(i).at(j).get<double>()));
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j) {
        INFO("|eta| " << e.at("absEta") << " p " << e.at("p") << " (" << i << "," << j << ")");
        REQUIRE(same(got[i][j], want.at(i).at(j).get<double>(), 1e-12, 1e-12 * mx));
      }
  }
}

TEST_CASE("refit v3 two-layer scores", "[refitV3]") {
  auto val = [](const nlohmann::json& x) { return x.get<double>() >= 1e29 ? kInf : x.get<double>(); };
  for (const auto& e : reference().at("tree2")) {
    const std::array<double, 2> x1{{val(e.at("x1").at(0)), val(e.at("x1").at(1))}};
    std::array<std::array<double, 2>, 3> x2;
    const bool two = !e.at("x2").is_null();
    if (two)
      for (int c = 0; c < 3; ++c)
        x2[c] = {{val(e.at("x2").at(c).at(0)), val(e.at("x2").at(c).at(1))}};
    const auto sc = Engine::tree2Scores(x1, two ? &x2 : nullptr, e.at("P").get<double>());
    for (int c = 0; c < 3; ++c) {
      const double w = val(e.at("scores").at(c));
      REQUIRE((std::isinf(w) ? std::isinf(sc[c]) : same(sc[c], w, 1e-15)));
    }
    REQUIRE(static_cast<int>(std::min_element(sc.begin(), sc.end()) - sc.begin()) == e.at("choice").get<int>());
  }
}
