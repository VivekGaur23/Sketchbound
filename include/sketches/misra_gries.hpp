#pragma once
#include <cstdint>
#include <unordered_map>

namespace sketches {

// Misra-Gries frequent-items summary with k counters. Deterministic guarantee:
// f(x) - N/(k+1) <= estimate(x) <= f(x). Any item with f(x) > N/(k+1) is
// therefore retained, which is the heavy-hitter property.
class MisraGries {
 public:
  explicit MisraGries(size_t counters) : counters_(counters) {}

  void Add(uint64_t key, uint64_t count = 1) {
    total_ += count;
    auto it = table_.find(key);
    if (it != table_.end()) {
      it->second += count;
      return;
    }
    if (table_.size() < counters_) {
      table_.emplace(key, count);
      return;
    }
    // Decrement all counters by the smallest of {count, all current values}.
    uint64_t d = count;
    for (const auto& kv : table_) d = std::min(d, kv.second);
    for (auto it2 = table_.begin(); it2 != table_.end();)
      if ((it2->second -= d) == 0) it2 = table_.erase(it2); else ++it2;
    if (count > d) table_.emplace(key, count - d);
  }

  uint64_t Estimate(uint64_t key) const {
    auto it = table_.find(key);
    return it == table_.end() ? 0 : it->second;
  }

  double ErrorBound() const {
    return static_cast<double>(total_) / (counters_ + 1);
  }

  const std::unordered_map<uint64_t, uint64_t>& table() const { return table_; }
  uint64_t total() const { return total_; }

 private:
  size_t counters_;
  uint64_t total_ = 0;
  std::unordered_map<uint64_t, uint64_t> table_;
};

}  // namespace sketches
