/*
 * voidfinder.h - Public API of the C "back-in-time" void finder.
 *
 * This library is a pure-C port of the LaZeVo back-in-time void finder
 * implemented in CosmoBolognaLib (Catalogue/VoidCatalogue.cpp, the
 * `Catalogue(VoidAlgorithm::_LaZeVo_, ...)` constructor).
 *
 * The algorithm identifies and counts cosmic voids in a catalogue of
 * tracers distributed in a cubic comoving box. It works in three steps:
 *
 *   1. Reconstruction (LaZeVo): a uniform random catalogue is generated
 *      inside the tracer box and tracers are optimally paired with random
 *      points by iteratively minimising pair distances, until the fraction
 *      of changed pairings drops below a threshold. The tracer->random
 *      offsets form a displacement field.
 *   2. Divergence field: the displacement field is deposited on a regular
 *      grid, Gaussian-smoothed, and the negative-divergence local minima
 *      are taken as void centres.
 *   3. Rescaling: each void centre is assigned a radius from the tracer
 *      density-contrast profile around it.
 *
 * Only the cubic-box geometry is supported (comoving X/Y/Z), matching the
 * original LaZeVo constructor; the light-cone variant is not ported.
 */

#ifndef VOIDFINDER_H
#define VOIDFINDER_H

#ifdef __cplusplus
extern "C" {
#endif

/* Return codes. */
typedef enum {
    VF_OK = 0,
    VF_ERR_IO = 1,           /* could not read/write a file               */
    VF_ERR_EMPTY = 2,        /* tracer catalogue has no objects           */
    VF_ERR_PARAM = 3,        /* invalid parameter                         */
    VF_ERR_NOMEM = 4,        /* memory allocation failed                  */
    VF_ERR_NO_VOIDS = 5      /* no local minima / voids found             */
} vf_status;

/*
 * Configuration of a void-finder run. Use vf_default_params() to obtain a
 * struct pre-filled with the values used in the CBL cosmicVoids example,
 * then override what you need.
 */
typedef struct {
    int    n_rec;        /* number of LaZeVo realisations (e.g. 5)            */
    double cellsize;     /* chain-mesh cell size; <=0 -> auto = 4 * mps      */
    double step_size;    /* divergence-grid step, in units of mps (e.g. 2.5/3)*/
    double threshold;    /* LaZeVo convergence threshold (e.g. 0.0)          */
    unsigned long seed;  /* RNG seed; 0 -> derive a seed internally          */

    int col_x;           /* 1-indexed column of X in the tracer file (e.g. 1)*/
    int col_y;           /* 1-indexed column of Y (e.g. 2)                   */
    int col_z;           /* 1-indexed column of Z (e.g. 3)                   */

    int verbose;         /* 0 = silent, 1 = progress messages               */
} vf_params;

/* Fill *p with sensible defaults (mirrors the CBL voidFinder example). */
void vf_default_params(vf_params *p);

/*
 * Run the back-in-time void finder.
 *
 *   tracer_file        path to a whitespace-separated tracer catalogue;
 *                      columns col_x/col_y/col_z give comoving X/Y/Z.
 *   output_voids_file  path to write the resulting void catalogue; each line
 *                      is "x y z radius". May be NULL to skip writing.
 *   params             configuration (NULL -> defaults).
 *   n_voids_out        if non-NULL, receives the number of voids found.
 *
 * Returns VF_OK on success or a vf_status error code otherwise.
 */
vf_status vf_run_void_finder(const char *tracer_file,
                             const char *output_voids_file,
                             const vf_params *params,
                             int *n_voids_out);

/* Human-readable message for a vf_status code. */
const char *vf_strerror(vf_status status);

#ifdef __cplusplus
}
#endif

#endif /* VOIDFINDER_H */
