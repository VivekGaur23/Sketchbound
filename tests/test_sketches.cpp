// Correctness checks for the four sketches and the hash families.
// Build: g++ -std=c++20 -O2 -Iinclude -Ibench -o tests/test tests/test_sketches.cpp
#include <cassert>
#include <cmath>
#include <cstdio>
#include <set>
#include <unordered_map>

#include "sketches/count_min.hpp"
#include "sketches/hyperloglog.hpp"
#include "sketches/misra_gries.hpp"
#include "sketches/reservoir.hpp"
#include "streams.hpp"

using namespace sketches;

void TestTwoUniversalCollisionRate() {
  // Pr[h(x) == h(y)] <= 1/w for x != y. Estimated over random pairs.
  std::mt19937_64 rng(1);
  TwoUniversal h = TwoUniversal::Random(rng);
  const uint64_t w = 1024;
  size_t collisions = 0, trials = 200000;
  for (size_t i = 0; i < trials; ++i) {
    uint64_t x = rng(), y = rng();
    if (x != y) collisions += (h(x) & (w - 1)) == (h(y) & (w - 1));
  }
  double rate = double(collisions) / trials;
  assert(rate < 2.0 / w);
  std::printf("  2-universal collision rate %.5f vs 1/w = %.5f\n", rate, 1.0 / w);
}

void TestWeakHashRowsAreCorrelated() {
  // The failure this library is built to demonstrate: with a power-of-two width
  // every row of the weak family collides on exactly the same pairs.
  std::mt19937_64 rng(1);
  WeakHash a = WeakHash::Random(rng), b = WeakHash::Random(rng);
  const uint64_t w = 4096;
  for (uint64_t x = 0; x < 500; ++x) {
    uint64_t y = x + w;  // congruent mod w
    assert((a(x) & (w - 1)) == (a(y) & (w - 1)));
    assert((b(x) & (w - 1)) == (b(y) & (w - 1)));
  }
  std::puts("  weak hash: rows collide identically, as predicted");
}

void TestCountMinNeverUnderestimates() {
  auto stream = streams::Zipf(200000, 5000, 1.1, 3);
  auto freq = streams::Frequencies(stream);
  CountMin<> cm(1024, 4, 7);
  for (uint64_t x : stream) cm.Add(x);
  for (const auto& kv : freq) assert(cm.Estimate(kv.first) >= kv.second);
  std::puts("  count-min never underestimates");
}

void TestCountMinExactWhenUncrowded() {
  CountMin<> cm(4096, 4, 7);
  for (uint64_t x = 0; x < 20; ++x) cm.Add(x, x + 1);
  for (uint64_t x = 0; x < 20; ++x) assert(cm.Estimate(x) == x + 1);
  std::puts("  count-min exact on a stream far below capacity");
}

void TestCountMinBoundHolds() {
  const double eps = 0.001, delta = 0.01;
  auto stream = streams::Uniform(500000, 50000, 5);
  auto freq = streams::Frequencies(stream);
  auto cm = CountMin<>::FromError(eps, delta, 7);
  for (uint64_t x : stream) cm.Add(x);
  size_t violations = 0;
  for (const auto& kv : freq)
    violations += (cm.Estimate(kv.first) - kv.second) > eps * stream.size();
  assert(double(violations) / freq.size() <= delta);
  std::printf("  count-min bound: %zu/%zu keys over eps*N (delta = %.2f)\n",
              violations, freq.size(), delta);
}

void TestHyperLogLogWithinFiveSigma() {
  for (unsigned p : {8u, 12u, 14u}) {
    for (uint64_t n : {1000ULL, 100000ULL, 1000000ULL}) {
      HyperLogLog hll(p, 42);
      for (uint64_t i = 0; i < n; ++i) hll.Add(i);
      double rel = std::abs(hll.Estimate() - double(n)) / n;
      assert(rel < 5.0 * hll.StandardError());
    }
  }
  std::puts("  hyperloglog within 5 standard errors across p and n");
}

void TestHyperLogLogExactOnEmptyAndTiny() {
  HyperLogLog hll(14, 42);
  assert(hll.Estimate() == 0.0);
  for (uint64_t i = 0; i < 10; ++i) hll.Add(i);
  assert(std::abs(hll.Estimate() - 10.0) < 1.0);  // linear counting range
  std::puts("  hyperloglog exact on empty and tiny inputs");
}

void TestReservoirShapeAndMembership() {
  const size_t n = 1000, k = 50;
  Reservoir<uint64_t> r(k, 5);
  for (uint64_t i = 0; i < n; ++i) r.Offer(i);
  assert(r.sample().size() == k);
  assert(r.seen() == n);
  std::set<uint64_t> unique(r.sample().begin(), r.sample().end());
  assert(unique.size() == k);  // no item sampled twice
  for (uint64_t x : r.sample()) assert(x < n);

  Reservoir<uint64_t> small(k, 5);
  for (uint64_t i = 0; i < 10; ++i) small.Offer(i);
  assert(small.sample().size() == 10);  // shorter than capacity
  std::puts("  reservoir: size, uniqueness and short-stream behaviour");
}

void TestReservoirUniformity() {
  const size_t n = 500, k = 25;
  const int trials = 20000;
  std::vector<uint64_t> hits(n, 0);
  for (int t = 0; t < trials; ++t) {
    Reservoir<uint64_t> r(k, 100 + t);
    for (uint64_t i = 0; i < n; ++i) r.Offer(i);
    for (uint64_t x : r.sample()) ++hits[x];
  }
  double p = double(k) / n, sd = std::sqrt(p * (1 - p) / trials), maxz = 0;
  for (uint64_t h : hits) maxz = std::max(maxz, std::abs(h / double(trials) - p) / sd);
  assert(maxz < 5.0);
  std::printf("  reservoir inclusion: max |z| = %.2f over %zu positions\n", maxz, n);
}

void TestMisraGriesBound() {
  for (double skew : {1.1, 2.0}) {
    auto stream = streams::Zipf(200000, 5000, skew, 3);
    auto freq = streams::Frequencies(stream);
    MisraGries mg(64);
    for (uint64_t x : stream) mg.Add(x);
    for (const auto& kv : freq) {
      assert(mg.Estimate(kv.first) <= kv.second);  // never overestimates
      assert(double(kv.second) - mg.Estimate(kv.first) <= mg.ErrorBound());
      if (kv.second > mg.ErrorBound()) assert(mg.Estimate(kv.first) > 0);
    }
    assert(mg.table().size() <= 64);
  }
  std::puts("  misra-gries: bound holds and every heavy hitter is retained");
}

int main() {
  std::puts("hash");
  TestTwoUniversalCollisionRate();
  TestWeakHashRowsAreCorrelated();
  std::puts("count-min");
  TestCountMinNeverUnderestimates();
  TestCountMinExactWhenUncrowded();
  TestCountMinBoundHolds();
  std::puts("hyperloglog");
  TestHyperLogLogWithinFiveSigma();
  TestHyperLogLogExactOnEmptyAndTiny();
  std::puts("reservoir");
  TestReservoirShapeAndMembership();
  TestReservoirUniformity();
  std::puts("misra-gries");
  TestMisraGriesBound();
  std::puts("\nall tests passed");
  return 0;
}
