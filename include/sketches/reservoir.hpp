#pragma once
#include <cstdint>
#include <random>
#include <vector>

namespace sketches {

// Vitter's Algorithm R. After any prefix of n items, every item seen is in the
// sample with probability exactly min(1, k/n).
template <class T>
class Reservoir {
 public:
  Reservoir(size_t capacity, uint64_t seed = 1) : capacity_(capacity), rng_(seed) {
    sample_.reserve(capacity);
  }

  void Offer(const T& item) {
    ++seen_;
    if (sample_.size() < capacity_) {
      sample_.push_back(item);
      return;
    }
    std::uniform_int_distribution<uint64_t> d(0, seen_ - 1);
    uint64_t j = d(rng_);
    if (j < capacity_) sample_[j] = item;
  }

  const std::vector<T>& sample() const { return sample_; }
  uint64_t seen() const { return seen_; }

 private:
  size_t capacity_;
  uint64_t seen_ = 0;
  std::vector<T> sample_;
  std::mt19937_64 rng_;
};

}  // namespace sketches
