/*
 * cvoidfinder - back-in-time (LaZeVo) cosmic void finder in C.
 * Copyright (C) 2026 Ariel Corte
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <https://www.gnu.org/licenses/>.
 *
 * This file is a derivative work ported from CosmoBolognaLib
 * (https://gitlab.com/federicomarulli/CosmoBolognaLib), GPL-2.0-or-later,
 * Copyright (C) by Federico Marulli, Simone Sartori, Sofia Contarini,
 * Carlo Cannarozzo and Tommaso Ronconi.
 */
/*
 * voidfinder.c - C port of the LaZeVo "back-in-time" void finder from
 * CosmoBolognaLib (Catalogue/VoidCatalogue.cpp:56, the
 * Catalogue(VoidAlgorithm::_LaZeVo_, ...) constructor), box geometry only.
 *
 * Step 1: LaZeVo reconstruction -> displacement field
 * Step 2: divergence field -> void centres (negative-divergence minima)
 * Step 3: rescaling -> void radii
 */

#include "voidfinder.h"

#include "catalogue.h"
#include "chainmesh.h"
#include "rng.h"
#include "util.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define VF_LOG(p, ...)                                                          \
    do {                                                                       \
        if ((p)->verbose) {                                                    \
            fprintf(stderr, __VA_ARGS__);                                      \
        }                                                                      \
    } while (0)

void vf_default_params(vf_params *p)
{
    p->n_rec = 5;
    p->cellsize = -1.0;     /* auto: 4 * mps */
    p->step_size = 2.5 / 3.0;
    p->threshold = 0.0;
    p->seed = 0;            /* derive internally */
    p->col_x = 1;
    p->col_y = 2;
    p->col_z = 3;
    p->verbose = 1;
}

const char *vf_strerror(vf_status s)
{
    switch (s) {
    case VF_OK:          return "success";
    case VF_ERR_IO:      return "file I/O error";
    case VF_ERR_EMPTY:   return "empty tracer catalogue";
    case VF_ERR_PARAM:   return "invalid parameter";
    case VF_ERR_NOMEM:   return "out of memory";
    case VF_ERR_NO_VOIDS:return "no voids found";
    default:             return "unknown error";
    }
}

/* ----------------------------- helpers --------------------------------- */

/* sum of pair distances for the 4-point assignment H[j] <-> R[j] */
static double assign_cost(const vf_catalogue *tr, const vf_catalogue *rnd,
                          const unsigned H[4], const unsigned R[4])
{
    double d = 0.0;
    for (int j = 0; j < 4; j++)
        d += vf_euclidean_distance(tr->x[H[j]], rnd->x[R[j]],
                                   tr->y[H[j]], rnd->y[R[j]],
                                   tr->z[H[j]], rnd->z[R[j]]);
    return d;
}

/* ----------------------- step 1: LaZeVo ------------------------------- */

/* randoms_out: n_rec arrays of length n, the tracer->random pairing. */
static vf_status lazevo_reconstruction(const vf_params *p,
                                       const vf_catalogue *tracers,
                                       const vf_catalogue *randoms_cat,
                                       const vf_chainmesh *cm_tr,
                                       const vf_chainmesh *cm_rnd,
                                       double mps,
                                       vf_rng *rng,
                                       unsigned **randoms_out)
{
    const size_t n = tracers->n;
    const unsigned N_near = (unsigned)(3.0 * 3.0 * 4.0 * VF_PI); /* = 113 */

    /* near_part[i][0..N_near-1]: N_near nearest tracers to tracer i */
    VF_LOG(p, "  * searching close particles (N=%u per tracer)...\n", N_near);
    unsigned *near_part = (unsigned *)malloc(n * (size_t)N_near * sizeof(unsigned));
    if (!near_part) return VF_ERR_NOMEM;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        vf_query q;
        vf_query_init(&q);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (long long i = 0; i < (long long)n; i++)
            vf_chainmesh_nnearest_obj(cm_tr, (unsigned)i, N_near,
                                      &near_part[(size_t)i * N_near], &q);
        vf_query_free(&q);
    }

    /* initialise pairings as shuffled identities */
    for (int rec = 0; rec < p->n_rec; rec++) {
        for (size_t i = 0; i < n; i++) randoms_out[rec][i] = (unsigned)i;
        vf_shuffle_u(rng, randoms_out[rec], n);
    }

    /* ---- starting configuration: spatially coherent initial pairing ---- */
    VF_LOG(p, "  * setting starting configuration...\n");
    const double dist = 4.0 * mps;
    const double minX = cm_tr->lim[0][0] + dist / 10, maxX = cm_tr->lim[0][1] - dist / 10;
    const double minY = cm_tr->lim[1][0] + dist / 10, maxY = cm_tr->lim[1][1] - dist / 10;
    const double minZ = cm_tr->lim[2][0] + dist / 10, maxZ = cm_tr->lim[2][1] - dist / 10;

    vf_query qc, qr;
    vf_query_init(&qc);
    vf_query_init(&qr);
    unsigned *close_rnd = NULL;
    size_t close_rnd_cap = 0;

    for (int rec = 0; rec < p->n_rec; rec++) {
        vf_chainmesh tr_copy, rnd_copy;
        if (vf_chainmesh_clone(&tr_copy, cm_tr) != 0 ||
            vf_chainmesh_clone(&rnd_copy, cm_rnd) != 0) {
            vf_query_free(&qc); vf_query_free(&qr); free(close_rnd);
            free(near_part);
            return VF_ERR_NOMEM;
        }

        size_t removed = 0;
        /* safety cap: avoid an unbounded loop if coverage stalls */
        unsigned long long guard = 0, guard_max = (unsigned long long)n * 1000 + 100000;
        while (removed < n && guard++ < guard_max) {
            double pos[3] = {
                vf_rng_uniform_range(rng, minX, maxX),
                vf_rng_uniform_range(rng, minY, maxY),
                vf_rng_uniform_range(rng, minZ, maxZ)
            };
            vf_chainmesh_close_objects(&tr_copy, pos, dist, &qc);
            int S = qc.count;
            int tr_to_rmv = S < 100 ? S : 100;
            if (tr_to_rmv > 0) {
                if ((size_t)S > close_rnd_cap) {
                    close_rnd_cap = (size_t)S;
                    unsigned *nc = (unsigned *)realloc(close_rnd, close_rnd_cap * sizeof(unsigned));
                    if (!nc) { vf_chainmesh_free(&tr_copy); vf_chainmesh_free(&rnd_copy);
                               vf_query_free(&qc); vf_query_free(&qr); free(close_rnd);
                               free(near_part); return VF_ERR_NOMEM; }
                    close_rnd = nc;
                }
                vf_chainmesh_nnearest_pos(&rnd_copy, pos, (unsigned)S, close_rnd, &qr);

                /* working copies of the candidate index lists we erase from */
                /* qc.idx holds the close tracers; close_rnd the near randoms */
                int csize = S, rsize = S;
                removed += (size_t)tr_to_rmv;
                while (tr_to_rmv > 0) {
                    unsigned r1 = vf_rng_int_inclusive(rng, 0, (unsigned)tr_to_rmv - 1);
                    unsigned r2 = vf_rng_int_inclusive(rng, 0, (unsigned)tr_to_rmv - 1);
                    unsigned tr_idx = qc.idx[r1];
                    unsigned rn_idx = close_rnd[r2];
                    randoms_out[rec][tr_idx] = rn_idx;
                    vf_chainmesh_delete(&tr_copy, tr_idx);
                    vf_chainmesh_delete(&rnd_copy, rn_idx);
                    /* erase r1 from qc.idx, r2 from close_rnd (order-preserving) */
                    memmove(&qc.idx[r1], &qc.idx[r1 + 1], (size_t)(csize - r1 - 1) * sizeof(unsigned));
                    memmove(&close_rnd[r2], &close_rnd[r2 + 1], (size_t)(rsize - r2 - 1) * sizeof(unsigned));
                    csize--; rsize--;
                    tr_to_rmv--;
                }
            }
        }
        vf_chainmesh_free(&tr_copy);
        vf_chainmesh_free(&rnd_copy);
        if (guard >= guard_max)
            VF_LOG(p, "  ! starting configuration hit safety cap (rec %d)\n", rec + 1);
    }
    vf_query_free(&qc);
    vf_query_free(&qr);
    free(close_rnd);

    /* --------------------- convergence iterations ---------------------- */
    VF_LOG(p, "  * performing iterations...\n");
    unsigned *index = (unsigned *)malloc(n * sizeof(unsigned));
    char *index_bool = (char *)malloc(n * sizeof(char));
    if (!index || !index_bool) { free(index); free(index_bool); free(near_part); return VF_ERR_NOMEM; }

    for (int rec = 0; rec < p->n_rec; rec++) {
        for (size_t i = 0; i < n; i++) index[i] = (unsigned)i;
        double ratio = 1.0;
        unsigned n_iter = 0;
        while (ratio > p->threshold) {
            vf_shuffle_u(rng, index, n);
            memset(index_bool, 0, n);

            for (size_t i = 0; i < n; i++) {
                unsigned base = index[i];
                unsigned r1 = vf_rng_int_inclusive(rng, 1, N_near - 1);
                unsigned r2 = vf_rng_int_inclusive(rng, 1, N_near - 1);
                unsigned r3 = vf_rng_int_inclusive(rng, 1, N_near - 1);
                while (r1 == r2 || r1 == r3 || r2 == r3) {
                    r2 = vf_rng_int_inclusive(rng, 1, N_near - 1);
                    r3 = vf_rng_int_inclusive(rng, 1, N_near - 1);
                }
                const unsigned *np = &near_part[(size_t)base * N_near];
                unsigned H[4] = { np[0], np[r1], np[r2], np[r3] };
                unsigned R[4] = { randoms_out[rec][H[0]], randoms_out[rec][H[1]],
                                  randoms_out[rec][H[2]], randoms_out[rec][H[3]] };
                unsigned R_copy[4]; memcpy(R_copy, R, sizeof(R));
                unsigned R_def[4]; memcpy(R_def, R, sizeof(R));

                /* minimise total pair distance over permutations of R, scanning
                   from the current order up to the lexicographic maximum, as in
                   the original CBL do/while(next_permutation) loop. */
                double dmin = 1e30;
                do {
                    double d = assign_cost(tracers, randoms_cat, H, R);
                    if (d < dmin) { dmin = d; memcpy(R_def, R, sizeof(R)); }
                } while (vf_next_permutation_u(R, 4));

                if (memcmp(R_def, R_copy, sizeof(R)) != 0) index_bool[i] = 1;
                for (int j = 0; j < 4; j++) randoms_out[rec][H[j]] = R_def[j];
            }

            size_t changed = 0;
            for (size_t i = 0; i < n; i++) changed += (size_t)index_bool[i];
            ratio = (double)changed / (double)n;
            VF_LOG(p, "    realisation %d/%d | iter %u | ratio %.6f (thr %.6f)\r",
                   rec + 1, p->n_rec, n_iter, ratio, p->threshold);
            n_iter++;
        }
        VF_LOG(p, "\n");
    }

    free(index);
    free(index_bool);
    free(near_part);
    return VF_OK;
}

/* ----------------- steps 2 & 3: divergence + rescaling ---------------- */

static vf_status divergence_and_voids(const vf_params *p,
                                      const vf_catalogue *tracers,
                                      const vf_catalogue *randoms_cat,
                                      const vf_chainmesh *cm_tr,
                                      double mps,
                                      unsigned **randoms,
                                      const char *output_file,
                                      int *n_voids_out)
{
    const size_t n = tracers->n;
    const double density = vf_cat_numdensity(tracers);
    const double step = p->step_size * mps;

    double mn[3], mx[3];
    for (int a = 0; a < 3; a++) {
        double rmin = vf_cat_min(randoms_cat, a), tmin = vf_cat_min(tracers, a);
        double rmax = vf_cat_max(randoms_cat, a), tmax = vf_cat_max(tracers, a);
        mn[a] = (rmin < tmin ? rmin : tmin) - 0.5 * step;
        mx[a] = (rmax > tmax ? rmax : tmax) + 0.5 * step;
    }
    int nCells = 0;
    for (int a = 0; a < 3; a++) {
        int c = (int)((mx[a] - mn[a]) / step);
        if (c > nCells) nCells = c;
    }
    if (nCells < 3) return VF_ERR_PARAM;

    /* re-centre the grid so each axis spans exactly step*nCells */
    for (int a = 0; a < 3; a++) {
        double span = mx[a] - mn[a];
        double pad = (step * nCells - span) / 2.0;
        mn[a] -= pad;
        mx[a] += pad;
    }

    const size_t NC = (size_t)nCells;
    const size_t ncube = NC * NC * NC;
    double *Div = (double *)calloc(ncube, sizeof(double));
    double *Dcopy = (double *)malloc(ncube * sizeof(double));
    if (!Div || !Dcopy) { free(Div); free(Dcopy); return VF_ERR_NOMEM; }
#define DIDX(i, j, k) (((size_t)(i) * NC + (size_t)(j)) * NC + (size_t)(k))

    VF_LOG(p, "  * estimating the divergence field...\n");
    for (int rec = 0; rec < p->n_rec; rec++) {
        for (size_t i = 0; i < n; i++) {
            unsigned ri = randoms[rec][i];
            double P0[3] = { randoms_cat->x[ri], randoms_cat->y[ri], randoms_cat->z[ri] };
            double P1[3] = { tracers->x[i], tracers->y[i], tracers->z[i] };
            double displ[3] = { P0[0] - P1[0], P0[1] - P1[1], P0[2] - P1[2] };
            int inds[2][3];
            for (int j = 0; j < 3; j++) {
                inds[0][j] = (int)((P0[j] - mn[j]) / step);
                inds[1][j] = (int)((P1[j] - mn[j]) / step);
            }
            if (inds[0][0] == inds[1][0] && inds[0][1] == inds[1][1] &&
                inds[0][2] == inds[1][2])
                continue;

            for (int part = 0; part < 2; part++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 3; k++) {
                        if (displ[k] == 0.0) continue;
                        double face = (inds[part][k] + j) * step + mn[k];
                        double t = (P0[k] - face) / displ[k];
                        if (t < 0.0 || t > 1.0) continue;
                        double ic[3];
                        for (int a = 0; a < 3; a++)
                            ic[a] = P0[a] - (P0[a] - P1[a]) * t;
                        int a1 = (k == 0) ? 1 : 0;
                        int a2 = (k == 2) ? 1 : 2;
                        double lo1 = inds[part][a1] * step + mn[a1];
                        double hi1 = (inds[part][a1] + 1) * step + mn[a1];
                        double lo2 = inds[part][a2] * step + mn[a2];
                        double hi2 = (inds[part][a2] + 1) * step + mn[a2];
                        if (ic[a1] > lo1 && ic[a1] <= hi1 &&
                            ic[a2] > lo2 && ic[a2] <= hi2) {
                            int ii = inds[part][0], jj = inds[part][1], kk = inds[part][2];
                            if (ii < 0 || jj < 0 || kk < 0 ||
                                ii >= nCells || jj >= nCells || kk >= nCells)
                                continue;
                            if (part == 0) Div[DIDX(ii, jj, kk)] -= fabs(displ[k]) / step;
                            else           Div[DIDX(ii, jj, kk)] += fabs(displ[k]) / step;
                        }
                    }
                }
            }
        }
    }

    VF_LOG(p, "  * smoothing...\n");
    for (size_t i = 0; i < ncube; i++) Dcopy[i] = Div[i] / (double)p->n_rec;

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int k = 0; k < nCells; k++) {
                double sum = 0.0, sum_div = 0.0;
                int i0 = vf_imax(0, i - 1), i1 = vf_imin(nCells, i + 2);
                int j0 = vf_imax(0, j - 1), j1 = vf_imin(nCells, j + 2);
                int k0 = vf_imax(0, k - 1), k1 = vf_imin(nCells, k + 2);
                for (int ii = i0; ii < i1; ii++)
                    for (int jj = j0; jj < j1; jj++)
                        for (int kk = k0; kk < k1; kk++) {
                            double dd = vf_euclidean_distance(
                                mn[0] + step * (i + 0.5), mn[0] + step * (ii + 0.5),
                                mn[1] + step * (j + 0.5), mn[1] + step * (jj + 0.5),
                                mn[2] + step * (k + 0.5), mn[2] + step * (kk + 0.5));
                            double ex = exp(-dd * dd / (2.0 * step * step));
                            sum += ex;
                            sum_div += Dcopy[DIDX(ii, jj, kk)] * ex;
                        }
                Div[DIDX(i, j, k)] = sum_div / sum;
            }

    /* ---------------------- identify void centres --------------------- */
    VF_LOG(p, "  * identifying voids...\n");
    size_t vcap = 1024, vcount = 0;
    double *vx = (double *)malloc(vcap * sizeof(double));
    double *vy = (double *)malloc(vcap * sizeof(double));
    double *vz = (double *)malloc(vcap * sizeof(double));
    if (!vx || !vy || !vz) { free(Div); free(Dcopy); free(vx); free(vy); free(vz); return VF_ERR_NOMEM; }

    size_t neg_cells = 0;
    for (int i = 1; i < nCells - 1; i++)
        for (int j = 1; j < nCells - 1; j++)
            for (int k = 1; k < nCells - 1; k++) {
                double center = Div[DIDX(i, j, k)];
                if (center >= 0.0) continue;
                neg_cells++;
                int control = 0;
                double num[3] = {0, 0, 0}, den[3] = {0, 0, 0};
                for (int ii = i - 1; ii <= i + 1 && !control; ii++)
                    for (int jj = j - 1; jj <= j + 1 && !control; jj++)
                        for (int kk = k - 1; kk <= k + 1 && !control; kk++) {
                            if (ii == i && jj == j && kk == k) continue;
                            double nb = Div[DIDX(ii, jj, kk)];
                            double cd = vf_euclidean_distance((double)i, (double)ii,
                                                              (double)j, (double)jj,
                                                              (double)k, (double)kk);
                            int dc[3] = { i - ii, j - jj, k - kk };
                            for (int t = 0; t < 3; t++) {
                                num[t] += nb * (step / cd) * dc[t];
                                den[t] += fabs(nb);
                            }
                            if (nb <= center || nb > 0.0) control = 1;
                        }
                if (!control) {
                    if (vcount == vcap) {
                        vcap *= 2;
                        double *nx = realloc(vx, vcap * sizeof(double));
                        double *ny = realloc(vy, vcap * sizeof(double));
                        double *nz = realloc(vz, vcap * sizeof(double));
                        if (!nx || !ny || !nz) { free(Div); free(Dcopy);
                            free(nx ? nx : vx); free(ny ? ny : vy); free(nz ? nz : vz);
                            return VF_ERR_NOMEM; }
                        vx = nx; vy = ny; vz = nz;
                    }
                    vx[vcount] = mn[0] + step * (i + 0.5) + num[0] / den[0];
                    vy[vcount] = mn[1] + step * (j + 0.5) + num[1] / den[1];
                    vz[vcount] = mn[2] + step * (k + 0.5) + num[2] / den[2];
                    vcount++;
                }
            }

    free(Div);
    free(Dcopy);
#undef DIDX

    VF_LOG(p, "  cells with negative divergence: %zu\n", neg_cells);
    if (vcount == 0) { free(vx); free(vy); free(vz); return VF_ERR_NO_VOIDS; }
    VF_LOG(p, "  number of voids: %zu\n", vcount);

    /* ------------------------- rescale radii -------------------------- */
    VF_LOG(p, "  * rescaling voids...\n");
    double *radius = (double *)malloc(vcount * sizeof(double));
    if (!radius) { free(vx); free(vy); free(vz); return VF_ERR_NOMEM; }

    vf_query qc;
    vf_query_init(&qc);
    double *dd0 = NULL, *dd1 = NULL, *dsort = NULL;
    size_t prof_cap = 0;

    for (size_t v = 0; v < vcount; v++) {
        double pos[3] = { vx[v], vy[v], vz[v] };
        vf_chainmesh_close_objects(cm_tr, pos, 5.0 * mps, &qc);
        int nclose = qc.count;
        if (nclose < 3) { radius[v] = 0.0; continue; }

        if ((size_t)nclose > prof_cap) {
            prof_cap = (size_t)nclose;
            dd0 = realloc(dd0, prof_cap * sizeof(double));
            dd1 = realloc(dd1, prof_cap * sizeof(double));
            dsort = realloc(dsort, prof_cap * sizeof(double));
            if (!dd0 || !dd1 || !dsort) { free(vx); free(vy); free(vz); free(radius);
                free(dd0); free(dd1); free(dsort); vf_query_free(&qc); return VF_ERR_NOMEM; }
        }

        /* NOTE: the original CBL code sorts the tracer index list but uses the
           distance array in its (unsorted) collection order when building the
           cumulative radial density profile, which makes the profile ill-defined.
           Here we sort the distances ascending, which is what a cumulative
           density-contrast profile requires. This is the one intentional
           deviation from the literal CBL source (see README). */
        for (int j = 0; j < nclose; j++) dsort[j] = qc.dist[j];
        /* insertion-ish shell sort ascending */
        for (int gap = nclose / 2; gap > 0; gap /= 2)
            for (int a = gap; a < nclose; a++) {
                double tmp = dsort[a];
                int b = a;
                while (b >= gap && dsort[b - gap] > tmp) { dsort[b] = dsort[b - gap]; b -= gap; }
                dsort[b] = tmp;
            }

        for (int j = 0; j < nclose - 1; j++) {
            dd0[j] = (dsort[j] + dsort[j + 1]) / 2.0;
            dd1[j] = (double)(j + 1) / (vf_volume_sphere(dd0[j]) * density);
        }
        dd0[nclose - 1] = 2.0 * dd0[nclose - 2] - dd0[nclose - 3];
        dd1[nclose - 1] = (double)nclose / (vf_volume_sphere(dd0[nclose - 1]) * density);

        int N = 0;
        while (!(dd1[N] < p->threshold && dd1[N + 1] > p->threshold) && N < nclose - 1)
            N++;

        if (N == nclose - 1) {
            int idx = nclose - 1;
            double maxv = dd1[nclose - 1];
            for (int kk = 0; kk < nclose - 1; kk++) {
                if (dd1[kk] < dd1[kk + 1] && dd1[kk + 1] < 0.5 && dd1[kk + 1] > maxv) {
                    idx = kk;
                    maxv = dd1[kk + 1];
                }
            }
            radius[v] = vf_interp_linear2(0.5, 0.0, dd1[idx], 0.0, dd0[idx]);
        } else {
            radius[v] = vf_interp_linear2(0.5, dd1[N], dd1[N + 1], dd0[N], dd0[N + 1]);
        }
    }
    vf_query_free(&qc);
    free(dd0); free(dd1); free(dsort);

    /* ----------------------------- output ----------------------------- */
    if (output_file) {
        FILE *f = fopen(output_file, "w");
        if (!f) { free(vx); free(vy); free(vz); free(radius); return VF_ERR_IO; }
        for (size_t v = 0; v < vcount; v++)
            fprintf(f, "%.7g %.7g %.7g %.7g\n", vx[v], vy[v], vz[v], radius[v]);
        fclose(f);
        VF_LOG(p, "  wrote %zu voids to %s\n", vcount, output_file);
    }

    if (n_voids_out) *n_voids_out = (int)vcount;

    free(vx); free(vy); free(vz); free(radius);
    return VF_OK;
}

/* ------------------------------- driver -------------------------------- */

vf_status vf_run_void_finder(const char *tracer_file,
                             const char *output_voids_file,
                             const vf_params *params,
                             int *n_voids_out)
{
    vf_params p;
    if (params) p = *params; else vf_default_params(&p);
    if (n_voids_out) *n_voids_out = 0;
    if (p.n_rec < 1) return VF_ERR_PARAM;
    if (p.step_size <= 0.0) return VF_ERR_PARAM;

    VF_LOG(&p, "> > > BACK-IN-TIME VOID FINDER (LaZeVo, C port) < < <\n");

    vf_catalogue tracers = {0};
    if (vf_cat_read(&tracers, tracer_file, p.col_x, p.col_y, p.col_z) != 0)
        return VF_ERR_IO;
    if (tracers.n == 0) { vf_cat_free(&tracers); return VF_ERR_EMPTY; }
    VF_LOG(&p, "  read %zu tracers\n", tracers.n);

    const double mps = vf_cat_mps(&tracers);
    const double cellsize = (p.cellsize > 0.0) ? p.cellsize : 4.0 * mps;
    VF_LOG(&p, "  mps = %.6g, cellsize = %.6g, n_rec = %d\n", mps, cellsize, p.n_rec);

    uint64_t seed = p.seed;
    if (seed == 0) seed = (uint64_t)time(NULL) ^ ((uint64_t)tracers.n << 32);
    vf_rng rng;
    vf_rng_seed(&rng, seed);

    /* random catalogue (box, same N as tracers) */
    vf_catalogue randoms_cat = {0};
    if (vf_cat_random_box(&randoms_cat, &tracers, &rng) != 0) {
        vf_cat_free(&tracers); return VF_ERR_NOMEM;
    }

    /* chain meshes */
    vf_chainmesh cm_tr, cm_rnd;
    if (vf_chainmesh_build(&cm_tr, &tracers, cellsize) != 0) {
        vf_cat_free(&tracers); vf_cat_free(&randoms_cat); return VF_ERR_NOMEM;
    }
    if (vf_chainmesh_build(&cm_rnd, &randoms_cat, cellsize) != 0) {
        vf_chainmesh_free(&cm_tr); vf_cat_free(&tracers); vf_cat_free(&randoms_cat);
        return VF_ERR_NOMEM;
    }

    /* pairings */
    unsigned **randoms = (unsigned **)malloc((size_t)p.n_rec * sizeof(unsigned *));
    if (!randoms) { vf_chainmesh_free(&cm_tr); vf_chainmesh_free(&cm_rnd);
        vf_cat_free(&tracers); vf_cat_free(&randoms_cat); return VF_ERR_NOMEM; }
    int alloc_ok = 1;
    for (int r = 0; r < p.n_rec; r++) {
        randoms[r] = (unsigned *)malloc(tracers.n * sizeof(unsigned));
        if (!randoms[r]) alloc_ok = 0;
    }

    vf_status st = VF_OK;
    if (!alloc_ok) st = VF_ERR_NOMEM;

    if (st == VF_OK)
        st = lazevo_reconstruction(&p, &tracers, &randoms_cat, &cm_tr, &cm_rnd,
                                   mps, &rng, randoms);

    if (st == VF_OK)
        st = divergence_and_voids(&p, &tracers, &randoms_cat, &cm_tr, mps,
                                  randoms, output_voids_file, n_voids_out);

    for (int r = 0; r < p.n_rec; r++) free(randoms[r]);
    free(randoms);
    vf_chainmesh_free(&cm_tr);
    vf_chainmesh_free(&cm_rnd);
    vf_cat_free(&tracers);
    vf_cat_free(&randoms_cat);
    return st;
}
