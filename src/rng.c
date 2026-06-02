#include "rng.h"

static inline uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline uint64_t rotl(uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

void vf_rng_seed(vf_rng *r, uint64_t seed)
{
    uint64_t sm = seed ? seed : 0x123456789ABCDEFULL;
    for (int i = 0; i < 4; i++)
        r->s[i] = splitmix64(&sm);
}

/* xoshiro256** */
uint64_t vf_rng_next_u64(vf_rng *r)
{
    const uint64_t result = rotl(r->s[1] * 5, 7) * 9;
    const uint64_t t = r->s[1] << 17;
    r->s[2] ^= r->s[0];
    r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2];
    r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3] = rotl(r->s[3], 45);
    return result;
}

double vf_rng_uniform(vf_rng *r)
{
    /* 53-bit mantissa -> [0,1) */
    return (vf_rng_next_u64(r) >> 11) * (1.0 / 9007199254740992.0);
}

double vf_rng_uniform_range(vf_rng *r, double lo, double hi)
{
    return lo + (hi - lo) * vf_rng_uniform(r);
}

uint64_t vf_rng_below(vf_rng *r, uint64_t n)
{
    if (n == 0) return 0;
    /* Lemire-style rejection for unbiased [0,n). */
    uint64_t t = (-n) % n; /* 2^64 mod n */
    uint64_t x;
    do {
        x = vf_rng_next_u64(r);
    } while (x < t);
    return x % n;
}

unsigned vf_rng_int_inclusive(vf_rng *r, unsigned lo, unsigned hi)
{
    if (hi <= lo) return lo;
    return lo + (unsigned)vf_rng_below(r, (uint64_t)(hi - lo) + 1);
}

void vf_shuffle_u(vf_rng *r, unsigned *a, size_t n)
{
    for (size_t i = n; i > 1; i--) {
        size_t j = (size_t)vf_rng_below(r, (uint64_t)i);
        unsigned tmp = a[i - 1];
        a[i - 1] = a[j];
        a[j] = tmp;
    }
}

int vf_next_permutation_u(unsigned *a, size_t n)
{
    if (n < 2) return 0;
    size_t i = n - 1;
    for (;;) {
        size_t i1 = i;
        if (a[--i] < a[i1]) {
            size_t i2 = n;
            while (!(a[i] < a[--i2]))
                ;
            unsigned tmp = a[i];
            a[i] = a[i2];
            a[i2] = tmp;
            /* reverse [i1, n) */
            for (size_t lo = i1, hi = n - 1; lo < hi; lo++, hi--) {
                unsigned t = a[lo];
                a[lo] = a[hi];
                a[hi] = t;
            }
            return 1;
        }
        if (i == 0) {
            /* reverse whole range -> smallest permutation */
            for (size_t lo = 0, hi = n - 1; lo < hi; lo++, hi--) {
                unsigned t = a[lo];
                a[lo] = a[hi];
                a[hi] = t;
            }
            return 0;
        }
    }
}
