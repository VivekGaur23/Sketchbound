#pragma once
#include <bit>
#include <cmath>
#include <cstdint>
#include <vector>

#include "sketches/hash.hpp"

namespace sketches {

// Flajolet et al. HyperLogLog with a 64-bit hash. Standard error of the
// corrected estimate is 1.04/sqrt(m) for m = 2^precision registers.
// RawEstimate() is exposed without the small-range correction so the bias it
// removes can be measured; see REPORT.md, experiment 3.
class HyperLogLog {
 public:
  explicit HyperLogLog(unsigned precision, uint64_t seed = 1)
      : p_(precision), m_(1ULL << precision), seed_(seed), registers_(m_, 0) {}

  void Add(uint64_t key) {
    uint64_t h = Mix64(key ^ seed_);
    uint64_t index = h >> (64 - p_);
    uint64_t tail = (h << p_) | (1ULL << (p_ - 1));
    uint8_t rho = static_cast<uint8_t>(std::countl_zero(tail) + 1);
    registers_[index] = std::max(registers_[index], rho);
  }

  double RawEstimate() const {
    double inv = 0.0;
    for (uint8_t r : registers_) inv += std::ldexp(1.0, -r);
    return Alpha() * m_ * m_ / inv;
  }

  double Estimate() const {
    double e = RawEstimate();
    size_t zeros = 0;
    for (uint8_t r : registers_) zeros += (r == 0);
    if (e <= 2.5 * m_ && zeros > 0)
      return m_ * std::log(static_cast<double>(m_) / zeros);  // linear counting
    return e;
  }

  double StandardError() const { return 1.04 / std::sqrt(static_cast<double>(m_)); }
  size_t m() const { return m_; }
  size_t bytes() const { return registers_.size(); }

 private:
  double Alpha() const {
    switch (m_) {
      case 16: return 0.673;
      case 32: return 0.697;
      case 64: return 0.709;
      default: return 0.7213 / (1.0 + 1.079 / m_);
    }
  }

  unsigned p_;
  size_t m_;
  uint64_t seed_;
  std::vector<uint8_t> registers_;
};

}  // namespace sketches
