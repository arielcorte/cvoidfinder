/*
 * rng.h - small self-contained pseudo-random number generator and the
 * shuffle / next-permutation helpers used by the void finder.
 *
 * CosmoBolognaLib relies on std::mt19937 / GSL plus std::shuffle and
 * std::next_permutation. Since the LaZeVo finder is intrinsically random
 * (its seeds come from the wall clock), bit-for-bit reproduction of CBL
 * output is neither possible nor required; we use a good portable PRNG
 * (SplitMix64 -> xoshiro256**) and standard-library-equivalent algorithms.
 */

#ifndef VF_RNG_H
#define VF_RNG_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint64_t s[4];
} vf_rng;

/* Seed an RNG state from a 64-bit seed. */
void vf_rng_seed(vf_rng *r, uint64_t seed);

/* Uniform 64-bit integer. */
uint64_t vf_rng_next_u64(vf_rng *r);

/* Uniform double in [0, 1). */
double vf_rng_uniform(vf_rng *r);

/* Uniform double in [lo, hi). */
double vf_rng_uniform_range(vf_rng *r, double lo, double hi);

/* Uniform integer in [0, n) (unbiased). */
uint64_t vf_rng_below(vf_rng *r, uint64_t n);

/* Uniform integer in [lo, hi] inclusive. */
unsigned vf_rng_int_inclusive(vf_rng *r, unsigned lo, unsigned hi);

/* Fisher-Yates in-place shuffle of an unsigned array of length n. */
void vf_shuffle_u(vf_rng *r, unsigned *a, size_t n);

/*
 * std::next_permutation equivalent for an array of unsigned of length n.
 * Rearranges into the next lexicographically greater permutation and
 * returns 1; if already the largest, rearranges to the smallest and
 * returns 0.
 */
int vf_next_permutation_u(unsigned *a, size_t n);

#endif /* VF_RNG_H */
