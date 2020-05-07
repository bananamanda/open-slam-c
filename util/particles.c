#include <stdlib.h>

#include "particles.h"

particle_t *init_particles(int n, int nlm) {
  particle_t *p = malloc(sizeof(particle_t) * n);

  int i, j;
  for (i = 0; i < n; i++) {
    p[i].w = 1 / n;

    p[i].xv[0] = 0.0;
    p[i].xv[1] = 0.0;
    p[i].xv[2] = 0.0;

    p[i].Nf = 0;
    p[i].xf = malloc(sizeof(double *) * nlm);
    p[i].Pf = malloc(sizeof(double *) * nlm);
    for (j = 0; j < nlm; j++) {
      p[i].xf[j] = malloc(sizeof(double) * 2);
      p[i].Pf[j] = malloc(sizeof(double) * 4);
    }
  }

  return p;
}

void free_particles(particle_t *p, int n, int nlm) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < nlm; j++) {
      free(p[i].xf[j]);
      free(p[i].Pf[j]);
    }
    free(p[i].xf);
    free(p[i].Pf);
  }
  free(p);
}
