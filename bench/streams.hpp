#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <unordered_map>
#include <vector>

namespace streams {

using Stream = std::vector<uint64_t>;
using Exact = std::unordered_map<uint64_t, uint64_t>;

inline Stream Uniform(size_t n, uint64_t universe, uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<uint64_t> d(0, universe - 1);
  Stream s(n);
  for (auto& x : s) x = d(rng);
  return s;
}

inline Stream Zipf(size_t n, uint64_t universe, double skew, uint64_t seed) {
  std::vector<double> cdf(universe);
  double acc = 0.0;
  for (uint64_t i = 0; i < universe; ++i)
    cdf[i] = (acc += 1.0 / std::pow(i + 1.0, skew));
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> d(0.0, acc);
  Stream s(n);
  for (auto& x : s)
    x = std::lower_bound(cdf.begin(), cdf.end(), d(rng)) - cdf.begin();
  return s;
}

// Keys aligned to a stride: aligned pointers, IPs inside a fixed prefix, or
// timestamps on a fixed tick all look like this. Uniform over the aligned keys,
// so it is only adversarial with respect to the hash, not the frequencies.
inline Stream Strided(size_t n, uint64_t universe, uint64_t stride, uint64_t seed) {
  Stream s = Uniform(n, universe, seed);
  for (auto& x : s) x *= stride;
  return s;
}

inline Exact Frequencies(const Stream& s) {
  Exact f;
  f.reserve(s.size() / 4);
  for (uint64_t x : s) ++f[x];
  return f;
}

}  // namespace streams
