#ifndef OMP
#define OMP 0
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if OMP
#include <omp.h>
#else
#include "../util/fake_omp.h"
#endif

#include "../util/particles.h"
#include "../util/util.h"
#include "configfile.h"

#include "../util/cycletimer.h"
#include "../util/instrument.h"

extern int errno;
extern char *optarg;
extern int opterr, optind;

typedef struct {
  int print;
  int predict_noise;
  int seed_random;
  int nparts;
  FILE *landmarks;
  FILE *waypoints;
  FILE *path;
} options_t;

void compute_steering(double x[3], double **wp, int nwp, int *iwp, double minD,
                      double *G, double rateG, double maxG, double dt);
void predict_true(double (*x)[3], double V, double G, int WB, double dt);
void add_control_noise(double V, double G, double Q[2][2], double *Vn,
                       double *Gn);
void predict(particle_t *particle, double V, double G, double Q[2][2], int WB,
             double dt, int addrandom);
void get_observations(double x[3], double **landmarks, int nlm, int *idf,
                      int rmax, double (*observations)[][2],
                      int (*observationftag)[], int *nlmo);
void add_observation_noise(double (*observations)[][2], int nlmo,
                           double R[2][2]);
double compute_weight(particle_t *particle, double observations[][2], int nlmo,
                      double R[2][2]);
void resample_particles(particle_t *particles, int nparts, int Nmin,
                        int do_resample);

/**
 *
 *
 **/
int main(int argc, char *argv[]) {
  bool instrument = false;
  int thread_count = 0;

  int i, j;
  char *helpstr =
      "FastSLAM 1.0 \n a fast SLAM simulation written "
      "in C. to use, create a set of landmarks and vehicle waypoints (i.e. "
      "waypoints for the desired vehicle path.) the program 'frontend' may "
      "be "
      "used to create this simulated environment - type help frontend for "
      "more "
      "information \n the configuration of the simulator is managed by the "
      "file 'configfile.h'. to alter the parameters of the vehicle, sensors, "
      "etc. adjust this file. there are also several switches that control "
      "certain filter options. \n usage: simulator [-o] [-p] [-r] [-i] [-l "
      "landmarks] "
      "[-w waypoints] [-n particles] [-s pathfile] \n\t -p \t\t predict noise "
      "switch \n\t "
      "-r "
      "\t\t seed random switch \n\t -i \t\t instrumentation switch \n\t -l "
      "landmarks \t the landmarks file \n\t "
      "-w "
      "waypoints \t the waypoints file \n\t -n particles \t the number of "
      "particles \n\t -s path file \t the file to save path information to \n";

  int opt;
  options_t options = {0, 0, 0, NPARTICLES, stdin, stdin, NULL};
  opterr = 0;

  while ((opt = getopt(argc, argv, "opril:w:n:t:s:")) != EOF) {
    switch (opt) {
    case 'o':
      options.print += 1;
      break;
    case 'p':
      options.predict_noise += 1;
      break;
    case 'r':
      options.seed_random += 1;
      break;
    case 'i':
      instrument = true;
      break;
    case 'l':
      if (!(options.landmarks = fopen(optarg, "r"))) {
        perror("fopen failed for landmark file");
        fprintf(stderr, helpstr);
        return -1;
      }
      break;
    case 'w':
      if (!(options.waypoints = fopen(optarg, "r"))) {
        perror("fopen failed for waypoint file");
        fprintf(stderr, helpstr);
        return -1;
      }
      break;
    case 'n':
      if (!(options.nparts = atoi(optarg))) {
        perror("failed to parse int");
        fprintf(stderr, helpstr);
        return -1;
      }
      break;
    case 't':
      thread_count = atoi(optarg);
      break;
    case 's':
      if (!(options.path = fopen(optarg, "w"))) {
        perror("fopen failed for path file");
        fprintf(stderr, helpstr);
        return -1;
      }
      break;
    default:
      fprintf(stderr, helpstr);
      return -1;
    }
  }

  if (!options.landmarks || !options.waypoints) {
    fprintf(stderr, helpstr);
    return -1;
  }

#if DEBUG
  fprintf(stderr, "Running with %d threads. Max possible is %d\n", thread_count,
          omp_get_max_threads());
#endif

  if (!options.predict_noise)
    fprintf(stderr, "Sampling from predict noise is necessary for FastSLAM 1.0 "
                    "particle diversity \n");

  track_activity(instrument);
  START_ACTIVITY(ACTIVITY_STARTUP);
  omp_set_num_threads(thread_count);

  // parse landmarks
  int nlm = 0;
  char buf[1024];
  if (fgets(buf, 1024, options.landmarks) != NULL) {
    nlm = atoi(buf);
  }

  double **landmarks = malloc(sizeof(double *) * nlm);

  i = 0;
  while (i < nlm && fgets(buf, 1024, options.landmarks) != NULL) {
    landmarks[i] = malloc(sizeof(double) * 2);

    char *tmp;
    tmp = strtok(buf, " \t\n");
    landmarks[i][0] = atof(tmp);
    tmp = strtok(NULL, " \t\n");
    landmarks[i][1] = atof(tmp);

    i++;
  }

  // parse waypoints
  int nwp = 0;
  if (fgets(buf, 1024, options.waypoints) != NULL) {
    nwp = atoi(buf);
  }

  double **waypoints = malloc(sizeof(double *) * nwp);

  i = 0;
  while (i < nwp && fgets(buf, 1024, options.waypoints) != NULL) {
    waypoints[i] = malloc(sizeof(double) * 2);

    char *tmp;
    tmp = strtok(buf, " \t\n");
    waypoints[i][0] = atof(tmp);
    tmp = strtok(NULL, " \t\n");
    waypoints[i][1] = atof(tmp);

    i++;
  }

  // initializations
  particle_t *particles = init_particles(options.nparts, nlm);
  double xtrue[3] = {0.0, 0.0, 0.0};

  int nlmo = 0;
  double observations[nlm][2];
  int observationftag[nlm];

  double dt = DT_CONTROLS;
  double dtsum = 0.0;
  int ftag[nlm];
  // double da_table[nlm];
  int iwp = 0;    // index of the initial waypoint
  double G = 0.0; // initial steer angle
  double V = 1.0;
  int nloops = NUMBER_LOOPS;
  if (options.seed_random) {
    time_t t;
    srand((unsigned)time(&t));
  }

  int m = (SWITCH_INFLATE_NOISE ? 2 : 1);
  double Qe[2][2]; //, Re[2][2];
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      Qe[i][j] = m * Q[i][j];
      //   Re[i][j] = m * R[i][j];
    }
  }

  for (i = 0; i < nlm; i++) {
    ftag[i] = i;
    // da_table[i] = 0.0;
  }
  FINISH_ACTIVITY(ACTIVITY_STARTUP);

  // main loop
  while (iwp < nwp) {
    // compute true data
    int oldwp = iwp;
    compute_steering(xtrue, waypoints, nwp, &iwp, AT_WAYPOINT, &G, RATEG, MAXG,
                     dt);
    if (oldwp != iwp && options.print)
      printf("Reached waypoint %d \t [%5.2f %5.2f]\n", oldwp,
             waypoints[oldwp][0], waypoints[oldwp][1]);

    if (iwp > nwp && nloops > 1) {
      iwp = 0;
      nloops -= 1;
    }
    predict_true(&xtrue, V, G, WHEELBASE, dt);

    // add process noise
    double Vn, Gn;
    if (SWITCH_CONTROL_NOISE) {
      add_control_noise(V, G, Q, &Vn, &Gn);
    } else {
      Vn = V;
      Gn = G;
    }

    // predict step
    START_ACTIVITY(ACTIVITY_PREDICT);
#if OMP
#pragma omp parallel
#endif
    {
#if OMP
#pragma omp for schedule(static)
#endif
      for (i = 0; i < options.nparts; i++) {
        predict(&(particles[i]), Vn, Gn, Qe, WHEELBASE, dt,
                options.predict_noise);
        if (SWITCH_HEADING_KNOWN)
          particles[i].xv[2] = xtrue[2];
#if DEBUG
        fprintf(stderr, "%d : %2.2f %2.2f %2.2f [%2.2f]\n", i,
                particles[i].xv[0], particles[i].xv[1], particles[i].xv[2],
                particles[i].w);
#endif
      }
    }
    FINISH_ACTIVITY(ACTIVITY_PREDICT);

    // observe step
    dtsum += dt;
    if (dtsum >= DT_OBSERVE) {
      dtsum = 0;
      // compute true data, then add noise
      nlmo = 0;
      get_observations(xtrue, landmarks, nlm, ftag, MAX_RANGE, &observations,
                       &observationftag, &nlmo);

      if (SWITCH_CONTROL_NOISE)
        add_observation_noise(&observations, nlmo, R);

      // data associate known

      // update
      for (i = 0; i < options.nparts; i++) {
        double w = compute_weight(&(particles[i]), observations, nlmo, R);
        particles[i].w *= w;
        // feature_update(&(particles[i]), observations, idf, R);
      }

      resample_particles(particles, options.nparts, NEFFECTIVE,
                         SWITCH_RESAMPLE);
    }
  }
  if (options.print)
    printf("Program end.\n");

  // free resources
  for (i = 0; i < nlm; i++)
    free(landmarks[i]);
  free(landmarks);

  for (i = 0; i < nwp; i++)
    free(waypoints[i]);
  free(waypoints);

  free_particles(particles, options.nparts, nlm);

  if (instrument)
    show_activity(stdout, nwp, nlm);
  return 0;
}

/**
 *
 **/
void compute_steering(double x[3], double **wp, int nwp, int *iwp, double minD,
                      double *G, double rateG, double maxG, double dt) {
  double *cwp = wp[*iwp];
  double d2 =
      (cwp[0] - x[0]) * (cwp[0] - x[0]) + (cwp[1] - x[1]) * (cwp[1] - x[1]);

  if (d2 < minD * minD) {
    *iwp += 1; // switch to next
    if (*iwp >= nwp)
      return;
    else
      cwp = wp[*iwp];
  }

  // compute change in G to point towards current waypoint
  double deltaG = atan2(cwp[1] - x[1], cwp[0] - x[0]) - x[2] - *G;
  pi_to_pi(&deltaG, 1, 0);

  // limit rate
  double maxDelta = rateG * dt;
  if (abs(deltaG) > maxDelta)
    deltaG = (deltaG > 0.0 ? 1.0 : -1.0) * maxDelta;

  // limit angle
  *G = *G + deltaG;
  if (abs(*G) > maxG)
    *G = (*G > 0.0 ? 1.0 : -1.0) * maxG;
}

/**
 *
 **/
void predict_true(double (*x)[3], double V, double G, int WB, double dt) {
  (*x)[0] += V * dt * cos(G + (*x)[2]);
  (*x)[1] += V * dt * sin(G + (*x)[2]);
  (*x)[2] += (V * dt * sin(G)) / WB;
  pi_to_pi((*x), 3, 2);
}

void add_control_noise(double V, double G, double Q[2][2], double *Vn,
                       double *Gn) {
  *Vn = V + (rand() / (double)RAND_MAX) * sqrt(Q[0][0]);
  *Gn = G + (rand() / (double)RAND_MAX) * sqrt(Q[1][1]);
}

void get_observations(double x[3], double **landmarks, int nlm, int *idf,
                      int range, double (*observations)[][2],
                      int (*observationidf)[], int *nlmo) {
  START_ACTIVITY(ACTIVITY_OBS);
  int dx, dy, i;
  for (i = 0; i < nlm; i++) {
    dx = landmarks[i][0] - x[0];
    dy = landmarks[i][1] - x[1];
    if ((abs(dx) < range) && (abs(dy) < range) &&
        (dx * cos(x[2]) + dy * sin(x[2]) > 0) &&
        (dx * dx + dy * dy < range * range)) {
      (*observations)[*nlmo][0] = dx * dx + dy * dy;
      (*observations)[*nlmo][1] = atan2(dy, dx) - x[2];
      (*observationidf)[*nlmo] = idf[i];
      *nlmo += 1;
    }
  }
  FINISH_ACTIVITY(ACTIVITY_OBS);
}

void add_observation_noise(double (*observations)[][2], int nlmo,
                           double R[2][2]) {
  START_ACTIVITY(ACTIVITY_OBS);
  int i;
  for (i = 0; i < nlmo; i++) {
    (*observations)[i][0] += rand() * sqrt(R[0][0]);
    (*observations)[i][1] += rand() * sqrt(R[1][1]);
  }
  FINISH_ACTIVITY(ACTIVITY_OBS);
}

void resample_particles(particle_t *particles, int nparts, int Nmin,
                        int doresample) {
  START_ACTIVITY(ACTIVITY_RESAMPLE);
  int i;
  double wsum = 0.0, wsumsq = 0.0;
  double *w = malloc(sizeof(double) * nparts);

  /*#if OMP
  #pragma omp parallel for reduction(+ : wsum) reduction(+ : wsumsq)
  #endif*/
  for (i = 0; i < nparts; i++) {
    w[i] = particles[i].w;
    wsum += w[i];
    wsumsq += w[i] * w[i];
  }

  /*#if OMP
  #pragma omp parallel for shared(w)
  #endif*/
  for (i = 0; i < nparts; i++) {
    w[i] = w[i] / wsum;
    particles[i].w = particles[i].w / wsum;
  }
  FINISH_ACTIVITY(ACTIVITY_RESAMPLE);
}

void compute_jacobians(particle_t *particle, int idf, double R[2][2],
                       double (*zp)[][2], double (*Hv)[][2][3],
                       double (*Hf)[][2][2], double (*Sf)[][2][2]) {
  double *xv, **xf, **Pf;
  xv = particle->xv;
  xf = particle->xf;
  Pf = particle->Pf;

  int i;
  for (i = 0; i < idf; i++) {
    double dx = xf[i][0] - xv[0];
    double dy = xf[i][1] - xv[1];
    double d2 = dx * dx + dy * dy;
    double d = sqrt(d2);

    (*zp)[i][0] = d;
    (*zp)[i][1] = atan2(dy, dx) - particle->xv[2];
    pi_to_pi((*zp)[i], 2, 1);

    // (*Hv)[i] = {{-dx / d, -dy / d, 0}, {dy / d2, -dx / d2, -1}};
    (*Hv)[i][0][0] = -dx / d;
    (*Hv)[i][0][1] = -dy / d;
    (*Hv)[i][0][2] = 0;
    (*Hv)[i][1][0] = dy / d2;
    (*Hv)[i][1][1] = -dx / d2;
    (*Hv)[i][1][2] = -1;

    // (*Hf)[i] = {{dx / d, dy / d}, {-dy / d2, dx / d2}};
    (*Hf)[i][0][0] = dx / d;
    (*Hf)[i][0][1] = dy / d;
    (*Hf)[i][1][0] = -dy / d2;
    (*Hf)[i][1][1] = dx / d2;

    // Sf[i] = Hf[i] * Pf[i] * Hf[i]' + R;
    (*Sf)[i][0][0] =
        (*Hf)[i][0][0] *
            (Pf[0][0] * (*Hf)[i][0][0] + Pf[0][1] * (*Hf)[i][0][1]) +
        (*Hf)[i][0][1] *
            (Pf[1][0] * (*Hf)[i][0][0] + Pf[1][1] * (*Hf)[i][0][1]) +
        R[0][0];
    (*Sf)[i][0][1] =
        (*Hf)[i][0][0] *
            (Pf[0][0] * (*Hf)[i][1][0] + Pf[0][1] * (*Hf)[i][1][1]) +
        (*Hf)[i][0][1] *
            (Pf[1][0] * (*Hf)[i][1][0] + Pf[1][1] * (*Hf)[i][1][1]) +
        R[0][1];
    (*Sf)[i][1][0] =
        (*Hf)[i][1][0] *
            (Pf[0][0] * (*Hf)[i][0][0] + Pf[0][1] * (*Hf)[i][0][1]) +
        (*Hf)[i][1][1] *
            (Pf[1][0] * (*Hf)[i][0][0] + Pf[1][1] * (*Hf)[i][0][1]) +
        R[1][0];
    (*Sf)[i][1][1] =
        (*Hf)[i][1][0] *
            (Pf[0][0] * (*Hf)[i][1][0] + Pf[0][1] * (*Hf)[i][1][1]) +
        (*Hf)[i][1][1] *
            (Pf[1][0] * (*Hf)[i][1][0] + Pf[1][1] * (*Hf)[i][1][1]) +
        R[2][2];
  }
}

double compute_weight(particle_t *particle, double observations[][2], int nlmo,
                      double R[2][2]) {
  double zp[nlmo][2], Hv[nlmo][2][3], Hf[nlmo][2][2], Sf[nlmo][2][2];
  compute_jacobians(particle, nlmo, R, &zp, &Hv, &Hf, &Sf);
  START_ACTIVITY(ACTIVITY_WEIGHTS);

  double v[nlmo][2];
  int i;
  for (i = 0; i < nlmo; i++) {
    v[i][0] = observations[i][0] - zp[i][0];
    v[i][1] = observations[i][1] - zp[i][1];
    pi_to_pi(v[i], 2, 1);
  }

  double w = 1.0;
  for (i = 0; i < nlmo; i++) {
    double S[2][2] = {{Sf[i][0][0], Sf[i][0][1]}, {Sf[i][1][0], Sf[i][1][1]}};

    double dom = 6.28 * sqrt(det(S));

    double r1[2][2], r2[2];
    inv(S, &r1);
    multiply2by1(r1, v[i], &r2);

    double num = exp(-0.5 * dot(v[i], r2));
    w *= num / dom;
  }
  FINISH_ACTIVITY(ACTIVITY_WEIGHTS);

  return w;
}

void predict(particle_t *particle, double V, double G, double Q[2][2], int WB,
             double dt, int addrandom) {
  if (addrandom)
    add_control_noise(V, G, Q, &V, &G);

  double *xv;
  xv = particle->xv;
  particle->xv[0] = xv[0] + V * dt * cos(G + xv[2]);
  particle->xv[1] = xv[1] + V * dt * sin(G + xv[2]);
  particle->xv[2] = xv[2] + V * dt * sin(G) / WB;
  pi_to_pi(particle->xv, 3, 2);
}
