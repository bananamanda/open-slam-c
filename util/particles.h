typedef struct {
  double w;
  double xv[3];
  
  int Nf;
  double **xf;
  double **Pf;
} particle_t;

particle_t *init_particles(int n, int nlm);
void free_particles(particle_t *p, int n, int nlm);
