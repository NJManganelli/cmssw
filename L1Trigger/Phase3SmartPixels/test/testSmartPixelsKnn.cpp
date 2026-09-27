// Exactness of smartpixels::KnnTree: same neighbours, same (d2, index) order as brute force,
// including exact ties (duplicated points) and a low-dimensional lattice.
#include "catch2/catch_all.hpp"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsKnn.h"

#include <algorithm>
#include <random>
#include <utility>
#include <vector>

namespace {
  std::vector<std::pair<double, uint32_t>> brute(const std::vector<double>& pts, unsigned dim, const double* q, unsigned k) {
    std::vector<std::pair<double, uint32_t>> all;
    for (uint32_t i = 0; i < pts.size() / dim; ++i)
      all.emplace_back(smartpixels::KnnTree::dist2(q, &pts[i * dim], dim), i);
    std::sort(all.begin(), all.end());
    all.resize(std::min<size_t>(k, all.size()));
    return all;
  }

  void check(const std::vector<double>& pts, unsigned dim, unsigned k, std::mt19937& rng, int nq, double lo, double hi) {
    const smartpixels::KnnTree tree(pts, dim);
    std::uniform_real_distribution<double> u(lo, hi);
    std::vector<std::pair<double, uint32_t>> got;
    std::vector<double> q(dim);
    for (int t = 0; t < nq; ++t) {
      for (auto& x : q)
        x = u(rng);
      tree.query(q.data(), k, got);
      REQUIRE(got == brute(pts, dim, q.data(), k));
    }
  }
}  // namespace

TEST_CASE("KnnTree matches brute force", "[SmartPixelsKnn]") {
  std::mt19937 rng(12345);
  std::normal_distribution<double> g(0., 1.);
  SECTION("11-D gaussian cloud, k = 10") {
    std::vector<double> pts(5000 * 11);
    for (auto& x : pts)
      x = g(rng);
    check(pts, 11, 10, rng, 300, -2., 2.);
  }
  SECTION("2-D with exact duplicates, k = 32") {
    std::vector<double> pts;
    for (int i = 0; i < 3000; ++i) {
      const double a = std::round(g(rng) * 20.), b = std::round(g(rng) * 20.);
      pts.push_back(a);
      pts.push_back(b);
    }
    check(pts, 2, 32, rng, 500, -40., 40.);
  }
  SECTION("fewer points than k") {
    std::vector<double> pts = {0., 0., 1., 1., 2., 2.};
    check(pts, 2, 10, rng, 10, -3., 3.);
  }
}
