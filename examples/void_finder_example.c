/*
 * void_finder_example.c - minimal driver for libvoidfinder.
 *
 * Usage:
 *   void_finder_example <tracers_file> <voids_output_file> [n_rec] [seed]
 *
 * The tracer file is whitespace-separated with X,Y,Z in columns 1,2,3
 * (the format of CosmoBolognaLib's Examples/cosmicVoids/input/halo_catalogue.txt).
 */

#include <stdio.h>
#include <stdlib.h>

#include "voidfinder.h"

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr,
                "usage: %s <tracers_file> <voids_output_file> [n_rec] [seed]\n",
                argv[0]);
        return 2;
    }

    vf_params p;
    vf_default_params(&p);
    if (argc >= 4) p.n_rec = atoi(argv[3]);
    if (argc >= 5) p.seed = strtoul(argv[4], NULL, 10);

    int n_voids = 0;
    vf_status st = vf_run_void_finder(argv[1], argv[2], &p, &n_voids);
    if (st != VF_OK) {
        fprintf(stderr, "void finder failed: %s\n", vf_strerror(st));
        return 1;
    }

    printf("Found %d voids -> %s\n", n_voids, argv[2]);
    return 0;
}
