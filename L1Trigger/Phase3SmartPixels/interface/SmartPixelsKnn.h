#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsKnn_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsKnn_h

// Exact k-nearest-neighbour search on a static point set (Euclidean, any dimension),
// for the SmartPixels-output angle emulator.
//
// Why not a CMSSW facility: CommonTools/RecoAlgos KDTreeLinkerAlgo and FKDTree answer
// box (range) queries only, and ROOT TKDTree's search is not const, so it cannot be
// shared read-only between streams. This tree is built once (GlobalCache) and every
// query is const and allocation-free apart from the caller's result buffer.
//
// DETERMINISTIC ORDER: neighbours are ordered by (squared distance, point index), and
// the squared distance is the sequential per-dimension sum of (q - p)^2, contraction-free
// (SPX_NO_FP_CONTRACT), exactly as the python reference (emu_common.knn_exact) computes
// it, so both select and order the same points. Pruning is exact: a subtree is skipped
// only when its splitting-plane distance is STRICTLY larger than the current k-th
// distance, so equal-distance points are still visited and the index tie-break holds.

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#if defined(__clang__)
#define SPX_NO_FP_CONTRACT
#define SPX_FP_CONTRACT_OFF_BODY _Pragma("clang fp contract(off)")
#elif defined(__GNUC__)
#define SPX_NO_FP_CONTRACT __attribute__((optimize("fp-contract=off")))
#define SPX_FP_CONTRACT_OFF_BODY
#else
#define SPX_NO_FP_CONTRACT
#define SPX_FP_CONTRACT_OFF_BODY
#endif

namespace smartpixels {

  class KnnTree {
  public:
    KnnTree() = default;
    // points: n x dim, row-major
    KnnTree(std::vector<double> points, unsigned dim) : pts_(std::move(points)), dim_(dim) {
      const uint32_t n = dim_ ? pts_.size() / dim_ : 0;
      perm_.resize(n);
      for (uint32_t i = 0; i < n; ++i)
        perm_[i] = i;
      nodes_.reserve(2 * (n / kLeaf + 1));
      if (n > 0)
        build(0, n);
    }

    uint32_t size() const { return perm_.size(); }
    unsigned dim() const { return dim_; }
    const double* point(uint32_t i) const { return &pts_[size_t(i) * dim_]; }

    SPX_NO_FP_CONTRACT static double dist2(const double* a, const double* b, unsigned dim) {
      SPX_FP_CONTRACT_OFF_BODY
      double d2 = 0.;
      for (unsigned j = 0; j < dim; ++j) {
        const double t = a[j] - b[j];
        d2 += t * t;
      }
      return d2;
    }

    // out: the k nearest (d2, index), ascending by (d2, index); fewer if the set is smaller.
    SPX_NO_FP_CONTRACT void query(const double* q, unsigned k, std::vector<std::pair<double, uint32_t>>& out) const {
      SPX_FP_CONTRACT_OFF_BODY
      out.clear();
      if (nodes_.empty() || k == 0)
        return;
      // `out` is kept as a max-heap on (d2, index) while searching
      const auto worse = [](const std::pair<double, uint32_t>& a, const std::pair<double, uint32_t>& b) {
        return a < b;
      };
      // (node, lower bound on the squared distance of any point below it)
      std::pair<uint32_t, double> stack[128];
      int sp = 0;
      stack[sp++] = {0, 0.};
      while (sp > 0) {
        const auto [ni, bound] = stack[--sp];
        if (out.size() == k && bound > out.front().first)
          continue;
        const Node& nd = nodes_[ni];
        if (nd.dim < 0) {
          for (uint32_t i = nd.lo; i < nd.hi; ++i) {
            const uint32_t id = perm_[i];
            const std::pair<double, uint32_t> c(dist2(q, point(id), dim_), id);
            if (out.size() < k) {
              out.push_back(c);
              std::push_heap(out.begin(), out.end(), worse);
            } else if (c < out.front()) {
              std::pop_heap(out.begin(), out.end(), worse);
              out.back() = c;
              std::push_heap(out.begin(), out.end(), worse);
            }
          }
          continue;
        }
        const double diff = q[nd.dim] - nd.split;
        stack[sp++] = {diff < 0 ? nd.right : nd.left, std::max(bound, diff * diff)};
        stack[sp++] = {diff < 0 ? nd.left : nd.right, bound};  // near side, popped first
      }
      std::sort_heap(out.begin(), out.end(), worse);
    }

  private:
    static constexpr uint32_t kLeaf = 16;
    struct Node {
      int dim;  // -1 = leaf
      double split;
      uint32_t lo, hi, left, right;
    };

    // median split on the widest dimension; left = coordinates <= split, right >= split
    uint32_t build(uint32_t lo, uint32_t hi) {
      const uint32_t id = nodes_.size();
      nodes_.push_back(Node{-1, 0., lo, hi, 0, 0});
      if (hi - lo <= kLeaf)
        return id;
      int best = 0;
      double spread = -1.;
      for (unsigned j = 0; j < dim_; ++j) {
        double mn = point(perm_[lo])[j], mx = mn;
        for (uint32_t i = lo + 1; i < hi; ++i) {
          const double v = point(perm_[i])[j];
          mn = std::min(mn, v);
          mx = std::max(mx, v);
        }
        if (mx - mn > spread) {
          spread = mx - mn;
          best = j;
        }
      }
      if (spread <= 0.)
        return id;  // all points identical: one leaf
      const uint32_t mid = lo + (hi - lo) / 2;
      std::nth_element(perm_.begin() + lo, perm_.begin() + mid, perm_.begin() + hi, [&](uint32_t a, uint32_t b) {
        return point(a)[best] < point(b)[best];
      });
      const double split = point(perm_[mid])[best];
      const uint32_t l = build(lo, mid);
      const uint32_t r = build(mid, hi);
      nodes_[id] = Node{best, split, lo, hi, l, r};
      return id;
    }

    std::vector<double> pts_;
    unsigned dim_ = 0;
    std::vector<uint32_t> perm_;
    std::vector<Node> nodes_;
  };

}  // namespace smartpixels

#endif
