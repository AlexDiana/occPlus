#ifndef OCCJSDM_RNG_H
#define OCCJSDM_RNG_H

// Thread-safe, R-seeded random number generation.
//
// Background. The Polya-Gamma sampler runs inside `#pragma omp parallel for`
// regions, and R's RNG (unif_rand/norm_rand) is a single global that is not
// safe to call from more than one thread. The original code called it anyway;
// the first fix replaced it with per-thread C++ engines, which closed the race
// but left the package unreproducible, because neither engine ever consulted
// R's RNG:
//
//   * get_rng() was seeded from the literal `12345 + omp_get_thread_num()`.
//     Deterministic per process, but set.seed() had no effect on it, and the
//     stream was never reset, so a second fit in the same session continued
//     where the first left off and could not be reproduced.
//   * mvrnormArmaQuick_TS() was seeded from std::random_device{}(), i.e. OS
//     entropy, so it was not reproducible even across processes. (It is also
//     not guaranteed to be non-deterministic: on some Windows toolchains
//     std::random_device is a fixed PRNG, which would have handed every thread
//     an identical seed and perfectly correlated the draws across the parallel
//     loop.)
//
// The scheme here. R draws one base seed per fit (see setOccJSDMSeed(), called
// from runOccJSDM()), which makes set.seed() control the whole sampler and
// makes consecutive fits in one session independent. Each thread then derives
// its own stream from that base seed. Threads never touch R's RNG.
//
// Note that `arma::randn()` elsewhere in this package is *not* covered by this
// header: RcppArmadillo routes Armadillo's RNG to R's, so those draws are
// already reproducible under set.seed(). They are also, for the same reason,
// only safe outside parallel regions -- which is where they currently sit.

#include <random>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

// Base seed and a generation counter, held as statics inside inline functions
// so that every translation unit shares one copy rather than getting its own.
//
// Both are written only by setOccJSDMSeed(), on the master thread, outside any
// parallel region; worker threads only ever read them. That is what makes the
// plain (non-atomic) reads below safe.
inline unsigned int& occjsdm_rng_base_seed() {
  static unsigned int seed = 12345u;
  return seed;
}

inline unsigned int& occjsdm_rng_generation() {
  static unsigned int generation = 1u;   // threads start at 0, so first use seeds
  return generation;
}

// One mt19937 per thread, re-seeded whenever R has handed us a new base seed.
//
// The generation counter is what makes re-seeding work at all: a thread_local
// engine is constructed once per thread and would otherwise keep its original
// seed for the life of the R session, so setting a new base seed between fits
// would have no effect.
inline std::mt19937& get_rng() {
  thread_local std::mt19937 rng;
  thread_local unsigned int seeded_at = 0u;   // 0 == never seeded

  const unsigned int current = occjsdm_rng_generation();
  if (seeded_at != current) {
    unsigned int tid = 0u;
#ifdef _OPENMP
    tid = static_cast<unsigned int>(omp_get_thread_num());
#endif
    // seed_seq rather than `base + tid`: mt19937 states from adjacent seeds
    // are correlated, and thread ids are about as adjacent as it gets.
    //
    // The generation is deliberately NOT seed material -- it only decides
    // *whether* to re-seed. Mixing it in would make the same base seed yield a
    // different stream on the second fit of a session, which defeats the whole
    // point. Consecutive fits already differ because runOccJSDM() draws each
    // base seed from R's RNG, which advances.
    std::seed_seq seq{ occjsdm_rng_base_seed(), tid };
    rng.seed(seq);
    seeded_at = current;
  }
  return rng;
}

inline double runif() {
  // std::uniform_real_distribution is stateless, so unlike rnorm() below it
  // needs no reset when the engine is re-seeded.
  static thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(get_rng());
}

inline double rnorm() {
  static thread_local std::normal_distribution<double> dist(0.0, 1.0);
  static thread_local unsigned int dist_reset_at = 0u;

  std::mt19937& rng = get_rng();   // re-seeds if R has set a new base seed

  // Re-seeding the engine is not enough. std::normal_distribution generates
  // deviates in pairs (Box-Muller) and caches the second one; that cache
  // survives a call to rng.seed(). Without the reset below, the first rnorm()
  // of a fit can return a value left over from the *previous* fit, so whether
  // a fit reproduces depends on the parity of the normal draws that happened
  // to precede it in the session. That is exactly what made this reproducible
  // when run alone but not when run after other tests in the same process.
  const unsigned int current = occjsdm_rng_generation();
  if (dist_reset_at != current) {
    dist.reset();
    dist_reset_at = current;
  }
  return dist(rng);
}

#endif // OCCJSDM_RNG_H
