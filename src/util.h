/*
 * util.h - small math helpers shared across the void finder, ported from
 * the corresponding CosmoBolognaLib functions.
 */

#ifndef VF_UTIL_H
#define VF_UTIL_H

#include <math.h>

#ifndef VF_PI
#define VF_PI 3.1415926535897932384626433832795
#endif

/* cbl::Euclidean_distance (Func.cpp:250) */
static inline double vf_euclidean_distance(double x1, double x2,
                                           double y1, double y2,
                                           double z1, double z2)
{
    const double dx = x1 - x2;
    const double dy = y1 - y2;
    const double dz = z1 - z2;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/* cbl::volume_sphere (Func.h:1231): 4/3 * pi * R^3 */
static inline double vf_volume_sphere(double r)
{
    return 4.0 / 3.0 * VF_PI * r * r * r;
}

/*
 * Linear interpolation/extrapolation between two points, matching how
 * cbl::interpolated(value, {x0,x1}, {y0,y1}, "Linear") behaves for the
 * 2-point inputs used by the void finder.
 */
static inline double vf_interp_linear2(double x, double x0, double x1,
                                       double y0, double y1)
{
    if (x1 == x0) return y0;
    return y0 + (x - x0) / (x1 - x0) * (y1 - y0);
}

/* round-half-to-even-ish nint used by the chain mesh (cbl::nint). */
static inline int vf_nint(double x)
{
    return (int)floor(x + 0.5);
}

static inline int vf_imax(int a, int b) { return a > b ? a : b; }
static inline int vf_imin(int a, int b) { return a < b ? a : b; }

#endif /* VF_UTIL_H */
