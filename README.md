# cvoidfinder — Back-in-time cosmic void finder (C)

A self-contained **C** port of the **LaZeVo "back-in-time" void finder** from
[CosmoBolognaLib](https://gitlab.com/federicomarulli/CosmoBolognaLib)
(`Catalogue/VoidCatalogue.cpp`, the
`Catalogue(VoidAlgorithm::_LaZeVo_, ...)` constructor).

The algorithm identifies and counts cosmic voids in a catalogue of tracers
(galaxies / haloes) distributed in a **cubic comoving box**.

## What it does

The finder runs in three steps, faithfully ported from the original C++:

1. **LaZeVo reconstruction.** A uniform random catalogue (same number of
   objects as the tracers) is generated inside the tracer box. Tracers are
   paired with random points; the pairing is first seeded by a spatially
   coherent pass and then refined by repeatedly taking quadruples of nearby
   tracers and re-assigning their random partners to minimise the total pair
   distance. Iteration continues until the fraction of pairings that change in
   a sweep drops to/below `threshold`. The tracer→random offsets form a
   **displacement field**. Several independent realisations (`n_rec`) are
   averaged.

2. **Divergence field.** The displacement field is deposited onto a regular
   grid using a line-crossing estimator, divided by `n_rec`, and
   Gaussian-smoothed over the 3×3×3 neighbourhood. Cells that are strict local
   minima of negative divergence (with all 26 neighbours non-positive) become
   **void centres**, with a sub-cell position from a divergence-weighted
   centroid.

3. **Rescaling.** Each void centre is given a **radius** from the tracer
   density-contrast profile around it, taken where the cumulative profile
   reaches half the mean density.

Only the **cubic box** geometry is supported (comoving X/Y/Z), matching the
original LaZeVo constructor. The light-cone variant (which needs cosmological
distance integration) is intentionally **not** ported.

## Build

Pure C11, no third-party libraries. OpenMP is optional.

```sh
make              # build libvoidfinder.a + example (serial)
make OPENMP=1     # build with OpenMP parallelism
make clean
```

This produces `libvoidfinder.a` and the `void_finder_example` program.

## Library API

Public header: [`include/voidfinder.h`](include/voidfinder.h).

```c
#include "voidfinder.h"

vf_params p;
vf_default_params(&p);   /* CBL cosmicVoids-example defaults */
p.n_rec     = 5;         /* number of LaZeVo realisations    */
p.step_size = 2.5/3.0;   /* divergence-grid step, in mps     */
p.threshold = 0.0;       /* convergence threshold            */
p.cellsize  = -1.0;      /* <=0 -> auto = 4 * mps            */
p.col_x = 1; p.col_y = 2; p.col_z = 3;  /* 1-indexed columns */
p.seed  = 0;             /* 0 -> derive from time + size      */

int n_voids = 0;
vf_status st = vf_run_void_finder("tracers.txt", "voids.txt", &p, &n_voids);
if (st != VF_OK) fprintf(stderr, "%s\n", vf_strerror(st));
```

`vf_run_void_finder` is the single entry point ("the method to call the
algorithm"):

| argument            | meaning                                                       |
|---------------------|--------------------------------------------------------------|
| `tracer_file`       | whitespace-separated catalogue; `col_x/y/z` give comoving X/Y/Z |
| `output_voids_file` | output path, one void per line: `x y z radius` (NULL to skip) |
| `params`            | configuration (NULL → defaults)                              |
| `n_voids_out`       | receives the number of voids found (may be NULL)             |

Returns `VF_OK` (0) on success, otherwise a `vf_status` error code
(`vf_strerror` gives a message).

## Input / output format

* **Input** — one tracer per line; X, Y, Z read from the configured columns.
  Blank lines and `#` comments are skipped, as are unparseable lines. This
  matches `CosmoBolognaLib/Examples/cosmicVoids/input/halo_catalogue.txt`.
* **Output** — one void per line: `x y z radius` (Mpc/h), the same layout as
  the CBL `Voids_*` output.

## Example

```sh
make OPENMP=1
OMP_NUM_THREADS=$(nproc) ./void_finder_example tracers.txt voids.txt 5 42
#                                               ^tracers  ^out  ^n_rec ^seed
```

## Notes on fidelity

* The LaZeVo finder is **intrinsically random** (the original seeds itself from
  the wall clock). Bit-for-bit reproduction of CBL output is therefore neither
  possible nor meaningful. This port uses a deterministic, seedable PRNG
  (SplitMix64 → xoshiro256**); **a fixed `seed` gives reproducible results**,
  and the serial and OpenMP builds agree.
* The convergence iteration is run **serially**. In CBL it carries an OpenMP
  pragma, but the loop body sets and immediately clears its `used[]` guards
  within each element, so under serialisation every element is always
  processed — i.e. the serial path *is* the intended per-element greedy update,
  without the data races the parallel version risks. The genuinely parallel,
  race-free step (per-tracer nearest-neighbour search) is parallelised when
  built with `OPENMP=1`.
* **One intentional correctness fix.** In the rescaling step the original CBL
  code sorts the *tracer-index* list but then builds the cumulative radial
  density profile from the *distance* array left in its unsorted collection
  order. A cumulative profile requires the distances in ascending order, so
  this port sorts the **distances**. This is the only deliberate deviation from
  the literal CBL source; see the comment in `src/voidfinder.c`
  (`divergence_and_voids`, rescaling loop).

## Layout

```
cvoidfinder/
├── include/voidfinder.h        public API
├── src/
│   ├── voidfinder.c            3-step algorithm + driver
│   ├── chainmesh.{c,h}         cubic-grid spatial index (CatalogueChainMesh)
│   ├── catalogue.{c,h}         point catalogue I/O, box stats, random box
│   ├── rng.{c,h}               PRNG, shuffle, next_permutation
│   └── util.h                  Euclidean distance, sphere volume, interpolation
├── examples/void_finder_example.c
└── Makefile
```

## Source provenance

Ported from CosmoBolognaLib:
* `Catalogue/VoidCatalogue.cpp:56` — the LaZeVo void-finder constructor (steps 1–3)
* `Catalogue/CatalogueChainMesh.cpp` — the spatial index
* `Catalogue/RandomCatalogue.cpp` — the cubic random-box generator
* `Headers/Catalogue.h` — `volume` / `numdensity` / `mps`
* `Func/Func.cpp`, `Headers/Func.h` — `Euclidean_distance`, `volume_sphere`,
  `interpolated`
