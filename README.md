# sketchbound

Header-only C++20 streaming sketches, with an experimental check of whether each
one meets the error bound it is proved to have.

Bounded-memory estimation over a stream you only see once: Count-Min for
frequencies, HyperLogLog for cardinality, reservoir sampling for a uniform
sample, Misra-Gries for heavy hitters. This is how network telemetry counts flows
and how a query planner estimates cardinality without materialising anything;
Redis ships HyperLogLog for exactly this reason.

The library is the small half. The point of the project is
**[REPORT.md](REPORT.md)** — measuring observed error against the predicted bound
across uniform, Zipfian and adversarially structured streams, and showing what
breaks when the hash family stops satisfying the hypothesis the proof depends on.

Two findings, both from the hash rather than the algorithm, and neither visible on
random input:

- Swapping Count-Min's 2-universal family for Knuth multiplicative hashing changes
  nothing on random or Zipfian keys — it is *exactly* correct on some of them —
  and then exceeds the guarantee by **1506x** on stride-aligned keys, with 100% of
  keys violating a bound meant to fail for 1%.
- HyperLogLog matches 1.04/sqrt(m) to three digits on random keys, and carries a
  systematic **+2.3% bias at p = 14** on sequential keys through a splitmix64
  finalizer, because a bijection is not a random oracle.

## Layout

```
include/sketches/   hash.hpp  count_min.hpp  hyperloglog.hpp
                    reservoir.hpp  misra_gries.hpp
bench/              streams.hpp (uniform / Zipf / stride), experiments.cpp
tests/              test_sketches.cpp
analysis/           plot.py
results/            CSV output
figures/            PNG output
```

## Build and reproduce

```sh
./build.sh
```

Compiles the tests and the harness, runs both, and regenerates every figure. The
whole sweep takes about ten seconds. Requires a C++20 compiler; the plots
additionally need Python with matplotlib. Nothing else — the library itself has no
dependencies beyond the standard library.

By hand:

```sh
g++ -std=c++20 -O2 -Iinclude -Ibench -o tests/test_sketches tests/test_sketches.cpp && ./tests/test_sketches
g++ -std=c++20 -O2 -Iinclude -Ibench -o bench/experiments bench/experiments.cpp && ./bench/experiments
python analysis/plot.py
```

## Use

```cpp
#include "sketches/count_min.hpp"

auto cm = sketches::CountMin<>::FromError(/*eps=*/0.001, /*delta=*/0.01);
cm.Add(key);
uint64_t f = cm.Estimate(key);   // f(key) <= f <= f(key) + eps*N w.p. 1 - delta
```

`CountMin` is templated on its hash family so the two can be compared under
identical conditions; `CountMin<>` is the 2-universal default and is the one to
use. `sketches::WeakHash` exists to be measured, not to be deployed.

## Tests

`tests/test_sketches.cpp` covers the properties, not just the API: Count-Min never
underestimates and is exact below capacity, its violation rate stays under delta,
the 2-universal family's measured collision rate matches 1/w, the weak family's
rows are shown to collide identically, HyperLogLog stays within five standard
errors across precisions and cardinalities, reservoir inclusion is uniform across
stream positions, and the Misra-Gries bound holds with every heavy hitter
retained. Plain asserts, no framework.
