#include "catalogue.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int vf_cat_alloc(vf_catalogue *c, size_t n)
{
    c->n = n;
    c->x = (double *)malloc(n * sizeof(double));
    c->y = (double *)malloc(n * sizeof(double));
    c->z = (double *)malloc(n * sizeof(double));
    if (n > 0 && (!c->x || !c->y || !c->z)) {
        vf_cat_free(c);
        return -1;
    }
    return 0;
}

void vf_cat_free(vf_catalogue *c)
{
    free(c->x);
    free(c->y);
    free(c->z);
    c->x = c->y = c->z = NULL;
    c->n = 0;
}

/* Parse the value at 1-indexed column `col` of a whitespace-separated line.
   Returns 0 on success. */
static int parse_column(const char *line, int col, double *out)
{
    const char *p = line;
    int cur = 0;
    while (*p) {
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (!*p) break;
        cur++;
        char *end = NULL;
        double v = strtod(p, &end);
        if (end == p) return -1; /* not a number */
        if (cur == col) {
            *out = v;
            return 0;
        }
        p = end;
    }
    return -1; /* column not found */
}

static int is_blank_or_comment(const char *line)
{
    const char *p = line;
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    return (*p == '\0' || *p == '#');
}

int vf_cat_read(vf_catalogue *c, const char *path,
                int col_x, int col_y, int col_z)
{
    FILE *f = fopen(path, "r");
    if (!f) return -1;

    size_t cap = 1024, n = 0;
    double *x = (double *)malloc(cap * sizeof(double));
    double *y = (double *)malloc(cap * sizeof(double));
    double *z = (double *)malloc(cap * sizeof(double));
    if (!x || !y || !z) {
        free(x); free(y); free(z); fclose(f);
        return -1;
    }

    char *line = NULL;
    size_t linecap = 0;
    ssize_t len;
    int rc = 0;
    while ((len = getline(&line, &linecap, f)) != -1) {
        if (is_blank_or_comment(line)) continue;
        double vx, vy, vz;
        if (parse_column(line, col_x, &vx) != 0 ||
            parse_column(line, col_y, &vy) != 0 ||
            parse_column(line, col_z, &vz) != 0) {
            continue; /* skip malformed lines, like CBL's tolerant readers */
        }
        if (n == cap) {
            cap *= 2;
            double *nx = (double *)realloc(x, cap * sizeof(double));
            double *ny = (double *)realloc(y, cap * sizeof(double));
            double *nz = (double *)realloc(z, cap * sizeof(double));
            if (!nx || !ny || !nz) { rc = -1; x = nx ? nx : x; y = ny ? ny : y; z = nz ? nz : z; break; }
            x = nx; y = ny; z = nz;
        }
        x[n] = vx; y[n] = vy; z[n] = vz;
        n++;
    }
    free(line);
    fclose(f);

    if (rc != 0) { free(x); free(y); free(z); return rc; }

    c->n = n;
    c->x = x;
    c->y = y;
    c->z = z;
    return 0;
}

static const double *axis_ptr(const vf_catalogue *c, int axis)
{
    return axis == 0 ? c->x : (axis == 1 ? c->y : c->z);
}

double vf_cat_min(const vf_catalogue *c, int axis)
{
    const double *a = axis_ptr(c, axis);
    double m = a[0];
    for (size_t i = 1; i < c->n; i++)
        if (a[i] < m) m = a[i];
    return m;
}

double vf_cat_max(const vf_catalogue *c, int axis)
{
    const double *a = axis_ptr(c, axis);
    double m = a[0];
    for (size_t i = 1; i < c->n; i++)
        if (a[i] > m) m = a[i];
    return m;
}

double vf_cat_volume(const vf_catalogue *c)
{
    return (vf_cat_max(c, 0) - vf_cat_min(c, 0)) *
           (vf_cat_max(c, 1) - vf_cat_min(c, 1)) *
           (vf_cat_max(c, 2) - vf_cat_min(c, 2));
}

double vf_cat_numdensity(const vf_catalogue *c)
{
    return (double)c->n / vf_cat_volume(c);
}

double vf_cat_mps(const vf_catalogue *c)
{
    return pow(vf_cat_numdensity(c), -1.0 / 3.0);
}

int vf_cat_random_box(vf_catalogue *out, const vf_catalogue *tracers,
                      vf_rng *rng)
{
    if (vf_cat_alloc(out, tracers->n) != 0) return -1;

    const double xmin = vf_cat_min(tracers, 0), xmax = vf_cat_max(tracers, 0);
    const double ymin = vf_cat_min(tracers, 1), ymax = vf_cat_max(tracers, 1);
    const double zmin = vf_cat_min(tracers, 2), zmax = vf_cat_max(tracers, 2);

    for (size_t i = 0; i < out->n; i++) {
        out->x[i] = vf_rng_uniform_range(rng, xmin, xmax);
        out->y[i] = vf_rng_uniform_range(rng, ymin, ymax);
        out->z[i] = vf_rng_uniform_range(rng, zmin, zmax);
    }
    return 0;
}
