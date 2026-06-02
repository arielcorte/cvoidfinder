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
 * chainmesh.h - cubic-grid spatial index, a C port of
 * cbl::catalogue::CatalogueChainMesh restricted to the operations used by
 * the LaZeVo void finder: nearest-neighbour and radius queries with
 * particle deletion.
 *
 * The grid references coordinates stored in a vf_catalogue (it does not own
 * them); a clone shares the same coordinates but keeps an independent set
 * of per-cell particle lists, so deletions on a clone are isolated.
 */

#ifndef VF_CHAINMESH_H
#define VF_CHAINMESH_H

#include <stddef.h>
#include "catalogue.h"

typedef struct {
    unsigned *idx;   /* particle indices in this cell */
    int count;
    int cap;
} vf_cell;

typedef struct {
    const vf_catalogue *cat; /* not owned */
    int ncells;              /* cells per dimension (m_dimension)        */
    double cellsize;         /* m_cellsize                               */
    double lim[3][2];        /* padded per-axis bounds                   */
    vf_cell *cells;          /* ncells^3 cells, flattened                */
    size_t ntot;             /* ncells^3                                 */
} vf_chainmesh;

/* Reusable scratch for queries (so hot loops avoid malloc churn). */
typedef struct {
    unsigned *idx;
    double   *dist;
    int       count;
    int       cap;
} vf_query;

void vf_query_init(vf_query *q);
void vf_query_free(vf_query *q);

/*
 * Build a chain mesh over catalogue `cat` with the requested cell size,
 * following CBL's cell-count derivation. Returns 0 on success.
 */
int vf_chainmesh_build(vf_chainmesh *cm, const vf_catalogue *cat,
                       double cellsize);

/* Deep-copy the per-cell lists (shares the catalogue). Returns 0 on success. */
int vf_chainmesh_clone(vf_chainmesh *dst, const vf_chainmesh *src);

void vf_chainmesh_free(vf_chainmesh *cm);

/* Remove particle `index` from its cell. */
void vf_chainmesh_delete(vf_chainmesh *cm, unsigned index);

/*
 * Collect all particles within [0, rMax) of pos into q (q->idx / q->count).
 * Equivalent to CBL Close_objects(pos, rMax).
 */
void vf_chainmesh_close_objects(const vf_chainmesh *cm, const double pos[3],
                                double rMax, vf_query *q);

/*
 * The N nearest particles to position pos, written (sorted by distance)
 * into out[0..N-1]. If fewer than N exist, the remaining slots are filled
 * with the farthest available index. q is scratch. Returns the count found.
 */
int vf_chainmesh_nnearest_pos(const vf_chainmesh *cm, const double pos[3],
                              unsigned N, unsigned *out, vf_query *q);

/*
 * The N nearest particles to particle `obj` (obj itself is included and is
 * element 0). Returns the count found.
 */
int vf_chainmesh_nnearest_obj(const vf_chainmesh *cm, unsigned obj,
                              unsigned N, unsigned *out, vf_query *q);

#endif /* VF_CHAINMESH_H */
