#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "sketches/hash.hpp"

namespace sketches {

// Cormode-Muthukrishnan Count-Min sketch over non-negative counts.
// Guarantee: f(x) <= estimate(x) <= f(x) + eps*N with probability >= 1 - delta,
// for width >= e/eps and depth >= ln(1/delta).
// Widths are rounded up to a power of two so indexing is a mask; this does not
// weaken the 2-universal bound, and it is what makes WeakHash fail loudly.
template <class Hash = TwoUniversal>
class CountMin {
 public:
  CountMin(size_t width, size_t depth, uint64_t seed = 1)
      : width_(RoundUpPow2(width)), depth_(depth), table_(width_ * depth, 0) {
    std::mt19937_64 rng(seed);
    for (size_t i = 0; i < depth_; ++i) hashes_.push_back(Hash::Random(rng));
  }

  static CountMin FromError(double eps, double delta, uint64_t seed = 1) {
    return CountMin(static_cast<size_t>(std::ceil(std::exp(1.0) / eps)),
                    static_cast<size_t>(std::ceil(std::log(1.0 / delta))), seed);
  }

  void Add(uint64_t key, uint64_t count = 1) {
    for (size_t i = 0; i < depth_; ++i) table_[i * width_ + Index(i, key)] += count;
  }

  uint64_t Estimate(uint64_t key) const {
    uint64_t best = UINT64_MAX;
    for (size_t i = 0; i < depth_; ++i)
      best = std::min(best, table_[i * width_ + Index(i, key)]);
    return best;
  }

  size_t width() const { return width_; }
  size_t depth() const { return depth_; }
  size_t bytes() const { return table_.size() * sizeof(uint64_t); }

  // Distinct cells a row actually uses. A 2-universal family drives this toward
  // the width; a family whose rows are correlated leaves it tiny.
  size_t Occupancy(size_t row) const {
    size_t used = 0;
    for (size_t j = 0; j < width_; ++j) used += table_[row * width_ + j] != 0;
    return used;
  }

 private:
  size_t Index(size_t row, uint64_t key) const {
    return hashes_[row](key) & (width_ - 1);
  }

  static size_t RoundUpPow2(size_t v) {
    size_t p = 1;
    while (p < v) p <<= 1;
    return p;
  }

  size_t width_, depth_;
  std::vector<uint64_t> table_;
  std::vector<Hash> hashes_;
};

}  // namespace sketches
