#include "chainmesh.h"
#include "util.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ----------------------------- query scratch ---------------------------- */

void vf_query_init(vf_query *q)
{
    q->idx = NULL;
    q->dist = NULL;
    q->count = 0;
    q->cap = 0;
}

void vf_query_free(vf_query *q)
{
    free(q->idx);
    free(q->dist);
    vf_query_init(q);
}

static int query_reserve(vf_query *q, int need)
{
    if (need <= q->cap) return 0;
    int cap = q->cap ? q->cap : 64;
    while (cap < need) cap *= 2;
    unsigned *ni = (unsigned *)realloc(q->idx, (size_t)cap * sizeof(unsigned));
    double   *nd = (double   *)realloc(q->dist, (size_t)cap * sizeof(double));
    if (!ni || !nd) { if (ni) q->idx = ni; if (nd) q->dist = nd; return -1; }
    q->idx = ni;
    q->dist = nd;
    q->cap = cap;
    return 0;
}

static int query_push(vf_query *q, unsigned i, double d)
{
    if (query_reserve(q, q->count + 1) != 0) return -1;
    q->idx[q->count] = i;
    q->dist[q->count] = d;
    q->count++;
    return 0;
}

/* ------------------------------- cells ---------------------------------- */

static int cell_push(vf_cell *c, unsigned v)
{
    if (c->count == c->cap) {
        int cap = c->cap ? c->cap * 2 : 4;
        unsigned *p = (unsigned *)realloc(c->idx, (size_t)cap * sizeof(unsigned));
        if (!p) return -1;
        c->idx = p;
        c->cap = cap;
    }
    c->idx[c->count++] = v;
    return 0;
}

/* ------------------------------- build ---------------------------------- */

static inline int clampi(int v, int lo, int hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

static inline int cell_coord(const vf_chainmesh *cm, double v, int axis)
{
    return (int)((v - cm->lim[axis][0]) / cm->cellsize);
}

static inline size_t cell_index(const vf_chainmesh *cm, int ix, int iy, int iz)
{
    return (size_t)ix * (size_t)cm->ncells * (size_t)cm->ncells +
           (size_t)iy * (size_t)cm->ncells + (size_t)iz;
}

int vf_chainmesh_build(vf_chainmesh *cm, const vf_catalogue *cat,
                       double cellsize)
{
    memset(cm, 0, sizeof(*cm));
    cm->cat = cat;

    const double start_cellsize = cellsize;
    double delta_max = 0.0;
    for (int a = 0; a < 3; a++) {
        double lo = vf_cat_min(cat, a) - 0.05 * start_cellsize;
        double hi = vf_cat_max(cat, a) + 0.05 * start_cellsize;
        cm->lim[a][0] = lo;
        cm->lim[a][1] = hi;
        double d = hi - lo;
        if (d > delta_max) delta_max = d;
    }

    int ncells = 0;
    double cs = delta_max;
    while (cs > start_cellsize) {
        ncells++;
        cs = delta_max / ncells;
    }
    if (ncells < 1) { ncells = 1; cs = delta_max; }

    cm->ncells = ncells;
    cm->cellsize = cs;
    cm->ntot = (size_t)ncells * ncells * ncells;
    cm->cells = (vf_cell *)calloc(cm->ntot, sizeof(vf_cell));
    if (!cm->cells) return -1;

    for (size_t i = 0; i < cat->n; i++) {
        int ix = clampi(cell_coord(cm, cat->x[i], 0), 0, ncells - 1);
        int iy = clampi(cell_coord(cm, cat->y[i], 1), 0, ncells - 1);
        int iz = clampi(cell_coord(cm, cat->z[i], 2), 0, ncells - 1);
        if (cell_push(&cm->cells[cell_index(cm, ix, iy, iz)], (unsigned)i) != 0)
            return -1;
    }
    return 0;
}

int vf_chainmesh_clone(vf_chainmesh *dst, const vf_chainmesh *src)
{
    *dst = *src;
    dst->cells = (vf_cell *)calloc(src->ntot, sizeof(vf_cell));
    if (!dst->cells) return -1;
    for (size_t i = 0; i < src->ntot; i++) {
        int n = src->cells[i].count;
        if (n == 0) continue;
        dst->cells[i].idx = (unsigned *)malloc((size_t)n * sizeof(unsigned));
        if (!dst->cells[i].idx) return -1;
        memcpy(dst->cells[i].idx, src->cells[i].idx, (size_t)n * sizeof(unsigned));
        dst->cells[i].count = n;
        dst->cells[i].cap = n;
    }
    return 0;
}

void vf_chainmesh_free(vf_chainmesh *cm)
{
    if (cm->cells) {
        for (size_t i = 0; i < cm->ntot; i++) free(cm->cells[i].idx);
        free(cm->cells);
    }
    cm->cells = NULL;
    cm->ntot = 0;
}

void vf_chainmesh_delete(vf_chainmesh *cm, unsigned index)
{
    const vf_catalogue *c = cm->cat;
    int ix = clampi(cell_coord(cm, c->x[index], 0), 0, cm->ncells - 1);
    int iy = clampi(cell_coord(cm, c->y[index], 1), 0, cm->ncells - 1);
    int iz = clampi(cell_coord(cm, c->z[index], 2), 0, cm->ncells - 1);
    vf_cell *cell = &cm->cells[cell_index(cm, ix, iy, iz)];
    for (int i = 0; i < cell->count; i++) {
        if (cell->idx[i] == index) {
            cell->idx[i] = cell->idx[cell->count - 1];
            cell->count--;
            return;
        }
    }
}

/* ------------------------------ queries --------------------------------- */

void vf_chainmesh_close_objects(const vf_chainmesh *cm, const double pos[3],
                                double rMax, vf_query *q)
{
    q->count = 0;
    const vf_catalogue *c = cm->cat;
    const int nc = cm->ncells;

    int cx = clampi(cell_coord(cm, pos[0], 0), 0, nc - 1);
    int cy = clampi(cell_coord(cm, pos[1], 1), 0, nc - 1);
    int cz = clampi(cell_coord(cm, pos[2], 2), 0, nc - 1);

    int h = (int)(rMax / cm->cellsize) + 2; /* cells to span, with margin */

    int x0 = clampi(cx - h, 0, nc - 1), x1 = clampi(cx + h, 0, nc - 1);
    int y0 = clampi(cy - h, 0, nc - 1), y1 = clampi(cy + h, 0, nc - 1);
    int z0 = clampi(cz - h, 0, nc - 1), z1 = clampi(cz + h, 0, nc - 1);

    for (int ix = x0; ix <= x1; ix++)
        for (int iy = y0; iy <= y1; iy++)
            for (int iz = z0; iz <= z1; iz++) {
                const vf_cell *cell = &cm->cells[cell_index(cm, ix, iy, iz)];
                for (int t = 0; t < cell->count; t++) {
                    unsigned k = cell->idx[t];
                    double d = vf_euclidean_distance(pos[0], c->x[k],
                                                     pos[1], c->y[k],
                                                     pos[2], c->z[k]);
                    if (d < rMax) query_push(q, k, d);
                }
            }
}

/* Build an index permutation of q sorted ascending by distance (shell sort). */
static void sort_query_indices(const vf_query *q, unsigned *order)
{
    for (int i = 0; i < q->count; i++) order[i] = (unsigned)i;
    int n = q->count;
    for (int gap = n / 2; gap > 0; gap /= 2) {
        for (int i = gap; i < n; i++) {
            unsigned tmp = order[i];
            int j = i;
            while (j >= gap && q->dist[order[j - gap]] > q->dist[tmp]) {
                order[j] = order[j - gap];
                j -= gap;
            }
            order[j] = tmp;
        }
    }
}

/* shared core for the two N-nearest variants, collecting in expanding cubes */
static int nnearest_core(const vf_chainmesh *cm, const double pos[3],
                         unsigned N, unsigned *out, vf_query *q)
{
    q->count = 0;
    const vf_catalogue *c = cm->cat;
    const int nc = cm->ncells;

    int cx = clampi(cell_coord(cm, pos[0], 0), 0, nc - 1);
    int cy = clampi(cell_coord(cm, pos[1], 1), 0, nc - 1);
    int cz = clampi(cell_coord(cm, pos[2], 2), 0, nc - 1);

    int h = 0;
    int prev_x0 = 1, prev_x1 = 0, prev_y0 = 1, prev_y1 = 0, prev_z0 = 1, prev_z1 = 0;

    for (;;) {
        int x0 = clampi(cx - h, 0, nc - 1), x1 = clampi(cx + h, 0, nc - 1);
        int y0 = clampi(cy - h, 0, nc - 1), y1 = clampi(cy + h, 0, nc - 1);
        int z0 = clampi(cz - h, 0, nc - 1), z1 = clampi(cz + h, 0, nc - 1);

        /* scan only the newly added shell of cells (avoid recounting) */
        for (int ix = x0; ix <= x1; ix++)
            for (int iy = y0; iy <= y1; iy++)
                for (int iz = z0; iz <= z1; iz++) {
                    int in_prev = (ix >= prev_x0 && ix <= prev_x1 &&
                                   iy >= prev_y0 && iy <= prev_y1 &&
                                   iz >= prev_z0 && iz <= prev_z1);
                    if (in_prev) continue;
                    const vf_cell *cell = &cm->cells[cell_index(cm, ix, iy, iz)];
                    for (int t = 0; t < cell->count; t++) {
                        unsigned k = cell->idx[t];
                        double d = vf_euclidean_distance(pos[0], c->x[k],
                                                         pos[1], c->y[k],
                                                         pos[2], c->z[k]);
                        query_push(q, k, d);
                    }
                }

        prev_x0 = x0; prev_x1 = x1;
        prev_y0 = y0; prev_y1 = y1;
        prev_z0 = z0; prev_z1 = z1;

        int covered_all = (x0 == 0 && y0 == 0 && z0 == 0 &&
                           x1 == nc - 1 && y1 == nc - 1 && z1 == nc - 1);

        if ((unsigned)q->count >= N) {
            /* guaranteed radius fully searched is h*cellsize */
            double safe = (double)h * cm->cellsize;
            /* find current Nth-nearest distance via a partial pass */
            /* (cheap because we only need to know if it is <= safe) */
            /* Count how many candidates are within safe. */
            int within = 0;
            for (int i = 0; i < q->count; i++)
                if (q->dist[i] <= safe) within++;
            if ((unsigned)within >= N || covered_all) break;
        } else if (covered_all) {
            break;
        }
        h++;
    }

    /* sort all collected candidates by distance, output the first N */
    unsigned *order = (unsigned *)malloc((size_t)q->count * sizeof(unsigned));
    if (!order) return 0;
    sort_query_indices(q, order);

    int outn = (int)N <= q->count ? (int)N : q->count;
    for (int i = 0; i < outn; i++) out[i] = q->idx[order[i]];
    /* pad if not enough (rare) with the farthest found */
    for (int i = outn; i < (int)N; i++)
        out[i] = q->count ? q->idx[order[q->count - 1]] : 0;

    free(order);
    return outn;
}

int vf_chainmesh_nnearest_pos(const vf_chainmesh *cm, const double pos[3],
                              unsigned N, unsigned *out, vf_query *q)
{
    return nnearest_core(cm, pos, N, out, q);
}

int vf_chainmesh_nnearest_obj(const vf_chainmesh *cm, unsigned obj,
                              unsigned N, unsigned *out, vf_query *q)
{
    const vf_catalogue *c = cm->cat;
    double pos[3] = { c->x[obj], c->y[obj], c->z[obj] };
    return nnearest_core(cm, pos, N, out, q);
}
