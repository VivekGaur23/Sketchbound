#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "sketches/count_min.hpp"
#include "sketches/hyperloglog.hpp"
#include "sketches/misra_gries.hpp"
#include "sketches/reservoir.hpp"
#include "streams.hpp"

using namespace sketches;

namespace {

struct Stats {
  double mean = 0, p50 = 0, p95 = 0, p99 = 0, max = 0;
  size_t violations = 0;
};

Stats Summarize(std::vector<double> errs, double bound) {
  Stats s;
  if (errs.empty()) return s;
  std::sort(errs.begin(), errs.end());
  for (double e : errs) {
    s.mean += e;
    s.violations += (e > bound);
  }
  s.mean /= errs.size();
  auto at = [&](double q) {
    return errs[std::min(errs.size() - 1, size_t(q * errs.size()))];
  };
  s.p50 = at(0.50);
  s.p95 = at(0.95);
  s.p99 = at(0.99);
  s.max = errs.back();
  return s;
}

struct Case {
  std::string name;
  streams::Stream stream;
};

std::vector<double> Overestimates(const streams::Exact& freq,
                                  const auto& sketch) {
  std::vector<double> errs;
  errs.reserve(freq.size());
  for (const auto& kv : freq)
    errs.push_back(double(sketch.Estimate(kv.first)) - double(kv.second));
  return errs;
}

// Experiment 1: does the eps*N guarantee hold, and how much slack is in it?
void CountMinBound(const std::vector<Case>& cases) {
  std::ofstream out("results/cm_bound.csv");
  out << "stream,eps,delta,width,depth,bytes,n,distinct,mean_err,p50,p95,p99,"
         "max_err,bound,violations,viol_frac,max_over_bound\n";
  const double delta = 0.01;
  for (const Case& c : cases) {
    auto freq = streams::Frequencies(c.stream);
    for (double eps : {0.01, 0.005, 0.001, 0.0005}) {
      auto cm = CountMin<>::FromError(eps, delta, 12345);
      for (uint64_t x : c.stream) cm.Add(x);
      double bound = eps * c.stream.size();
      Stats s = Summarize(Overestimates(freq, cm), bound);
      out << c.name << "," << eps << "," << delta << "," << cm.width() << ","
          << cm.depth() << "," << cm.bytes() << "," << c.stream.size() << ","
          << freq.size() << "," << s.mean << "," << s.p50 << "," << s.p95 << ","
          << s.p99 << "," << s.max << "," << bound << "," << s.violations << ","
          << double(s.violations) / freq.size() << "," << s.max / bound << "\n";
    }
  }
  std::puts("results/cm_bound.csv");
}

// Experiment 2: the bound assumes a 2-universal family. What does violating that
// assumption cost, and on which streams is the violation even visible?
template <class Hash>
void HashRow(std::ofstream& out, const Case& c, const streams::Exact& freq,
             size_t width, size_t depth) {
  CountMin<Hash> cm(width, depth, 12345);
  for (uint64_t x : c.stream) cm.Add(x);
  double bound = std::exp(1.0) / cm.width() * c.stream.size();
  Stats s = Summarize(Overestimates(freq, cm), bound);
  out << c.name << "," << Hash::name() << "," << cm.width() << "," << cm.depth()
      << "," << c.stream.size() << "," << freq.size() << "," << cm.Occupancy(0)
      << "," << s.mean << "," << s.max << "," << bound << "," << s.violations
      << "," << double(s.violations) / freq.size() << "," << s.max / bound
      << "\n";
}

void HashIndependence() {
  std::ofstream out("results/cm_hash.csv");
  out << "stream,hash,width,depth,n,distinct,cells_used_row0,mean_err,max_err,"
         "bound,violations,viol_frac,max_over_bound\n";
  const size_t n = 1000000, universe = 4096, width = 4096, depth = 5;
  std::vector<Case> cases;
  cases.push_back({"uniform4k", streams::Uniform(n, universe, 7)});
  cases.push_back({"uniform100k", streams::Uniform(n, 100000, 7)});
  cases.push_back({"zipf1.1", streams::Zipf(n, universe, 1.1, 7)});
  for (uint64_t stride : {8ULL, 64ULL, 4096ULL})
    cases.push_back({"stride" + std::to_string(stride),
                     streams::Strided(n, universe, stride, 7)});
  for (const Case& c : cases) {
    auto freq = streams::Frequencies(c.stream);
    HashRow<TwoUniversal>(out, c, freq, width, depth);
    HashRow<WeakHash>(out, c, freq, width, depth);
  }
  std::puts("results/cm_hash.csv");
}

// Experiment 3: measured relative error against the 1.04/sqrt(m) prediction,
// with and without the small-range correction, over two key sets. The HLL
// analysis assumes a random oracle, which is a strictly stronger assumption than
// the 2-universality Count-Min needs, so it is worth testing whether one round
// of splitmix64 delivers it on structured input.
void HyperLogLogError() {
  std::ofstream out("results/hll.csv");
  out << "p,m,bytes,n_true,keys,trials,rmse_rel,rmse_rel_raw,mean_rel_bias,"
         "predicted_se\n";
  const int trials = 100;
  for (unsigned p : {4u, 6u, 8u, 10u, 12u, 14u, 16u}) {
    for (uint64_t n : {100ULL, 1000ULL, 10000ULL, 100000ULL, 1000000ULL}) {
      for (int random_keys = 0; random_keys < 2; ++random_keys) {
        double sq = 0, sq_raw = 0, bias = 0, predicted = 0;
        size_t bytes = 0;
        for (int t = 0; t < trials; ++t) {
          HyperLogLog hll(p, 1000 + t);
          std::mt19937_64 rng(7777 + t);
          for (uint64_t i = 0; i < n; ++i) hll.Add(random_keys ? rng() : i);
          double rel = (hll.Estimate() - double(n)) / n;
          double rel_raw = (hll.RawEstimate() - double(n)) / n;
          sq += rel * rel;
          sq_raw += rel_raw * rel_raw;
          bias += rel;
          predicted = hll.StandardError();
          bytes = hll.bytes();
        }
        out << p << "," << (1ULL << p) << "," << bytes << "," << n << ","
            << (random_keys ? "random" : "sequential") << "," << trials << ","
            << std::sqrt(sq / trials) << "," << std::sqrt(sq_raw / trials) << ","
            << bias / trials << "," << predicted << "\n";
      }
    }
  }
  std::puts("results/hll.csv");
}

// Experiment 4: is every stream position really retained with probability k/n?
void ReservoirUniformity() {
  std::ofstream summary("results/reservoir.csv");
  summary << "n,k,trials,expected_p,min_freq,max_freq,max_abs_z\n";
  std::ofstream positions("results/reservoir_positions.csv");
  positions << "position,inclusion_freq,expected_p\n";
  const std::vector<std::pair<size_t, size_t>> configs = {
      {1000, 50}, {1000, 200}, {10000, 100}};
  for (const auto& cfg : configs) {
    size_t n = cfg.first, k = cfg.second;
    const int trials = 20000;
    std::vector<uint64_t> hits(n, 0);
    for (int t = 0; t < trials; ++t) {
      Reservoir<uint64_t> r(k, 999 + t);
      for (uint64_t i = 0; i < n; ++i) r.Offer(i);
      for (uint64_t x : r.sample()) ++hits[x];
    }
    double p = double(k) / n;
    double sd = std::sqrt(p * (1 - p) / trials);
    double lo = 1.0, hi = 0.0, maxz = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double f = double(hits[i]) / trials;
      lo = std::min(lo, f);
      hi = std::max(hi, f);
      maxz = std::max(maxz, std::abs(f - p) / sd);
      if (n == 1000 && k == 50) positions << i << "," << f << "," << p << "\n";
    }
    summary << n << "," << k << "," << trials << "," << p << "," << lo << ","
            << hi << "," << maxz << "\n";
  }
  std::puts("results/reservoir.csv");
}

// Experiment 5: the only deterministic bound in the set. It should never break.
void MisraGriesBound(const std::vector<Case>& cases) {
  std::ofstream out("results/misra_gries.csv");
  out << "stream,k,n,distinct,max_underest,bound,violations,heavy_true,"
         "heavy_found\n";
  for (const Case& c : cases) {
    auto freq = streams::Frequencies(c.stream);
    for (size_t k : {8, 32, 128, 512}) {
      MisraGries mg(k);
      for (uint64_t x : c.stream) mg.Add(x);
      double bound = mg.ErrorBound();
      double worst = 0;
      size_t violations = 0, heavy_true = 0, heavy_found = 0;
      for (const auto& kv : freq) {
        double under = double(kv.second) - double(mg.Estimate(kv.first));
        worst = std::max(worst, under);
        violations += (under > bound);
        if (double(kv.second) > bound) {
          ++heavy_true;
          heavy_found += (mg.Estimate(kv.first) > 0);
        }
      }
      out << c.name << "," << k << "," << c.stream.size() << "," << freq.size()
          << "," << worst << "," << bound << "," << violations << ","
          << heavy_true << "," << heavy_found << "\n";
    }
  }
  std::puts("results/misra_gries.csv");
}

}  // namespace

int main() {
  const size_t n = 1000000, universe = 100000;
  std::vector<Case> cases;
  cases.push_back({"uniform", streams::Uniform(n, universe, 7)});
  cases.push_back({"zipf1.1", streams::Zipf(n, universe, 1.1, 7)});
  cases.push_back({"zipf2.0", streams::Zipf(n, universe, 2.0, 7)});

  CountMinBound(cases);
  HashIndependence();
  HyperLogLogError();
  ReservoirUniformity();
  MisraGriesBound(cases);
  return 0;
}
