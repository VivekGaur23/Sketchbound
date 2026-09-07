#pragma once
#include <cstdint>
#include <random>

namespace sketches {

// Carter-Wegman: h(x) = ((a*x + b) mod (2^61 - 1)), a 2-universal family.
// Pr[h(x) == h(y)] <= 1/w after reduction mod w, which is the only property the
// Count-Min analysis needs.
class TwoUniversal {
 public:
  TwoUniversal() = default;
  TwoUniversal(uint64_t a, uint64_t b) : a_(a), b_(b) {}

  static TwoUniversal Random(std::mt19937_64& rng) {
    std::uniform_int_distribution<uint64_t> d(1, kPrime - 1);
    return TwoUniversal(d(rng), d(rng));
  }

  uint64_t operator()(uint64_t x) const {
    return Reduce(MulMod(a_, x % kPrime) + b_);
  }

  static const char* name() { return "two_universal"; }

 private:
  static constexpr uint64_t kPrime = (1ULL << 61) - 1;

  static uint64_t Reduce(uint64_t v) { return v >= kPrime ? v - kPrime : v; }

  static uint64_t MulMod(uint64_t x, uint64_t y) {
    __uint128_t z = static_cast<__uint128_t>(x) * y;
    return Reduce(static_cast<uint64_t>(z & kPrime) +
                  static_cast<uint64_t>(z >> 61));
  }

  uint64_t a_ = 1, b_ = 0;
};

// Knuth multiplicative hashing with a per-row additive seed. This is what an
// implementer reaching for "a fast hash" often writes, and it is not 2-universal:
// with a power-of-two table, C is odd and therefore invertible mod w, so
// h(x) == h(y) iff x == y (mod w) -- independently of the seed. Every row of a
// sketch collides on exactly the same pairs, so taking a minimum over rows buys
// nothing. See REPORT.md, experiment 2.
class WeakHash {
 public:
  WeakHash() = default;
  explicit WeakHash(uint64_t seed) : seed_(seed) {}

  static WeakHash Random(std::mt19937_64& rng) { return WeakHash(rng()); }

  uint64_t operator()(uint64_t x) const { return x * kKnuth + seed_; }

  static const char* name() { return "weak_multiplicative"; }

 private:
  static constexpr uint64_t kKnuth = 2654435761ULL;
  uint64_t seed_ = 0;
};

// splitmix64 finalizer: strong avalanche, used where the analysis assumes a
// random oracle rather than mere 2-universality (HyperLogLog).
inline uint64_t Mix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

}  // namespace sketches
