/*
 * catalogue.h - a minimal point catalogue (comoving X/Y/Z) and the box
 * statistics / random-catalogue generation used by the void finder,
 * ported from cbl::catalogue::Catalogue (box methods only).
 */

#ifndef VF_CATALOGUE_H
#define VF_CATALOGUE_H

#include <stddef.h>
#include "rng.h"

typedef struct {
    size_t n;
    double *x;
    double *y;
    double *z;
} vf_catalogue;

/* Allocate a catalogue of n points (coordinates uninitialised). */
int vf_cat_alloc(vf_catalogue *c, size_t n);

/* Free the coordinate arrays. */
void vf_cat_free(vf_catalogue *c);

/*
 * Read a whitespace-separated tracer file. col_x/col_y/col_z are 1-indexed
 * columns. Blank lines and lines starting with '#' are skipped.
 * Returns 0 on success, non-zero on I/O error.
 */
int vf_cat_read(vf_catalogue *c, const char *path,
                int col_x, int col_y, int col_z);

/* Per-axis minimum/maximum over the catalogue (axis: 0=x,1=y,2=z). */
double vf_cat_min(const vf_catalogue *c, int axis);
double vf_cat_max(const vf_catalogue *c, int axis);

/* Box volume from the bounding box (matches Catalogue::volume default). */
double vf_cat_volume(const vf_catalogue *c);

/* Number density n/V. */
double vf_cat_numdensity(const vf_catalogue *c);

/* Mean particle separation numdensity^(-1/3). */
double vf_cat_mps(const vf_catalogue *c);

/*
 * Generate a uniform random catalogue with the same number of objects as
 * `tracers`, drawn inside the tracer bounding box (cbl _createRandom_box_,
 * N_R = 1). Returns 0 on success.
 */
int vf_cat_random_box(vf_catalogue *out, const vf_catalogue *tracers,
                      vf_rng *rng);

#endif /* VF_CATALOGUE_H */
