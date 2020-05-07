/* Simulator for EKF SLAM
 */

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

#include "configfile.h"
//#include "particles.h"
#include "../util/util.h"

#include "../util/cycletimer.h"
#include "../util/instrument.h"

#if OMP
#include <omp.h>
#else
#include "../util/fake_omp.h"
#endif

extern int errno;
extern char *optarg;
extern int opterr, optind;

typedef struct {
  int print;
  int sensor_noise;
  int control_noise;
  int seed_random;
  FILE *landmarks;
  FILE *waypoints;
} options_t;

void compute_steering(double x[3], double **wp, int nwp, int *iwp, double minD,
                      double *G, double rateG, double maxG, double dt);
void vehicle_model(double (*x)[3], double V, double G, int WB, double dt);
void add_control_noise(double *V, double *G, double Q[2][2]);
void predict(double (*x)[5], double (*P)[3][3], double V, double G,
             double Q[2][2], int WB, double dt);
void get_observations(double x[5], double **landmarks, int *idf, int nlm,
                      int range, double (*observations)[][2],
                      int (*observationidf)[], int *nlmo);
void add_observation_noise(double (*observations)[][2], int nlmo,
                           double R[2][2]);
void update(double (*x)[5], double (*P)[3][3], double observations[][2],
            int nlmo, double R[2][2], int *idf);
void KF_update(double (*x)[5], double (*P)[3][3], double v[2], double R[2][2],
               double H[2][5]);
void make_symmetric(double P[3][3], double (*res)[3][3]);

/**
 *
 *
 **/
int main(int argc, char *argv[]) {
  bool instrument = false;
  int thread_count = 0;

  char *helpstr =
      "EKFSLAM 1.0 \n a fast SLAM simulation written "
      "in C. to use, create a set of landmarks and vehicle waypoints (i.e. "
      "waypoints for the desired vehicle path.) the program 'frontend' may "
      "be "
      "used to create this simulated environment - type help frontend for "
      "more "
      "information \n the configuration of the simulator is managed by the "
      "file 'configfile.h'. to alter the parameters of the vehicle, sensors, "
      "etc. adjust this file. there are also several switches that control "
      "certain filter options. \n usage: simulator [-o] [-p] [-c] [-r] [-i] "
      "[-l "
      "landmarks] "
      "[-w waypoints] [-t threads]\n\t -p \t\t sensor noise switch \n\t -c "
      "\t\t control "
      "noise "
      "switch \n\t -r "
      "\t\t seed random switch \n\t -i \t\t include instrumentation \n\t -l "
      "landmarks \t the landmarks file \n\t "
      "-w "
      "waypoints \t the waypoints file \n\t -t threads \t the number of "
      "threads\n";

  int opt;
  options_t options = {0, 0, 0, 0, stdin, stdin};
  opterr = 0;

  while ((opt = getopt(argc, argv, "opril:w:n:t:")) != EOF) {
    switch (opt) {
    case 'o':
      options.print += 1;
      break;
    case 'p':
      options.sensor_noise += 1;
      break;
    case 'c':
      options.control_noise += 1;
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
    case 't':
      thread_count = atoi(optarg);
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

  int i = 0;
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
  double xtrue[3] = {0.0, 0.0, 0.0};
  double x[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double P[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  double observations[nlm][2];
  int observationftag[nlm];
  double dt = DT_CONTROLS;
  double dtsum = 0.0;
  int ftag[nlm];
  double da_table[nlm];
  int iwp = 0;    // index of the initial waypoint
  double G = 0.0; // initial steer angle
  double V = 1.0;
  int nlmo = 0;

  if (options.seed_random) {
    time_t t;
    srand((unsigned)time(&t));
  }

  int m = (SWITCH_INFLATE_NOISE ? 2 : 1);
  double Qe[2][2], Re[2][2];
  int j;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      Qe[i][j] = m * Q[i][j];
      Re[i][j] = m * R[i][j];
    }
  }

  for (i = 0; i < nlm; i++) {
    observations[i][0] = 0.0;
    observations[i][1] = 0.0;
    observationftag[i] = 0;
    ftag[i] = i;
    da_table[i] = 0.0;
  }

  int nloops = NUMBER_LOOPS;
  FINISH_ACTIVITY(ACTIVITY_STARTUP);

  if (options.print)
    printf("Startup finished. Going to main loop\n");

  // main loop
  while (iwp < nwp) {
    int oldwp = iwp;
    // compute true data
    compute_steering(xtrue, waypoints, nwp, &iwp, AT_WAYPOINT, &G, RATEG, MAXG,
                     dt);
    if (iwp != oldwp && options.print)
      printf("Reached waypoint %f %f!\n", waypoints[oldwp][0],
             waypoints[oldwp][1]);

    if (iwp > nwp && nloops > 1) {
      iwp = 0;
      nloops--;
    }

    vehicle_model(&xtrue, V, G, WHEELBASE, dt); // vehicle_model

    // add process noise
    if (SWITCH_CONTROL_NOISE)
      add_control_noise(&V, &G, Q);

    // ekf predict step
    predict(&x, &P, V, G, Q, WHEELBASE, dt);

    // observe heading if we want to

    dtsum = dtsum + dt;
    if (dtsum >= DT_OBSERVE) {
      // get observations
      nlmo = 0;
      get_observations(x, landmarks, ftag, nlm, MAX_RANGE, &observations,
                       &observationftag, &nlmo);

      // add observation noise
      if (SWITCH_CONTROL_NOISE) {
        add_observation_noise(&observations, nlmo, R);
      }

      // data associate known

      // update
      update(&x, &P, observations, nlmo, Re, observationftag);

      // augment
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

  if (instrument)
    show_activity(stdout, nwp, nlm);

  return 0;
}

/**
 * Functions
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

void vehicle_model(double (*x)[3], double V, double G, int WB, double dt) {
  (*x)[0] += V * dt * cos(G + (*x)[2]);
  (*x)[1] += V * dt * sin(G + (*x)[2]);
  (*x)[2] += (V * dt * sin(G)) / WB;
  pi_to_pi(*x, 3, 2);
}

// multivariate_gauss(double *x, int lx, double **P, int lP, int n) {

void add_control_noise(double *V, double *G, double Q[2][2]) {
  *V = *V + rand() * sqrt(Q[0][0]);
  *G = *G + rand() * sqrt(Q[1][1]);
}

void predict(double (*x)[5], double (*P)[3][3], double V, double G,
             double Q[2][2], int WB, double dt) {
  START_ACTIVITY(ACTIVITY_PREDICT);
  double s = sin(G + (*x)[2]);
  double c = cos(G + (*x)[2]);
  double vts = V * dt * s;
  double vtc = V * dt * c;
  double Gv[3][3] = {{1, 0, -vts}, {0, 1, vtc}, {0, 0, 1}};
  double Gvt[3][3] = {{1, 0, 0}, {0, 1, 0}, {-vts, vtc, 1}};
  double Gu[3][2] = {
      {dt * c, -vts}, {dt * s, vtc}, {dt * sin(G) / WB, V * dt * cos(G) / WB}};
  double Gut[2][3] = {{dt * c, dt * s, dt * sin(G) / WB},
                      {-vts, vtc, V * dt * cos(G) / WB}};
  double temp1[3][3];
  double temp2[3][3];
  double temp3[3][2];
  double temp4[3][3];
  int i, j;
  multiply1(Gv, *P, &temp1);
  multiply1(temp1, Gvt, &temp2);
  multiply2(Gu, Q, &temp3);
  multiply3(temp3, Gut, &temp4);

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      (*P)[i][j] = temp2[i][j] + temp4[i][j];
    }
  }
  (*x)[0] += vtc;
  (*x)[1] += vts;
  (*x)[2] += V * dt * sin(G) / WB;
  pi_to_pi(*x, 3, 2);
  FINISH_ACTIVITY(ACTIVITY_PREDICT);
}

void get_observations(double x[5], double **landmarks, int *idf, int nlm,
                      int range, double (*observations)[][2],
                      int (*observationidf)[], int *nlmo) {
  START_ACTIVITY(ACTIVITY_OBS);

  int i;
  for (i = 0; i < nlm; i++) {
    double dx = landmarks[i][0] - x[0];
    double dy = landmarks[i][1] - x[1];
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

void update(double (*x)[5], double (*P)[3][3], double observations[][2],
            int nlmo, double R[2][2], int *idf) {
  START_ACTIVITY(ACTIVITY_UPDATE);
  int fpos = 3;
  double dx, dy, d;
  double z[2];
  double H[2][5];
  double v[2];

  /*#if OMP
  #pragma omp parallel
  #endif
    {*/
  int i;

#if OMP
//#pragma omp for schedule(dynamic)
#pragma omp parallel for
#endif
  for (i = 0; i < nlmo; i++) {
    // observe model
    dx = (*x)[fpos] - (*x)[0];
    dy = (*x)[fpos] - (*x)[1];
    d = sqrt(dx * dx + dy * dy);

    z[0] = d;
    z[1] = atan2(dy, dx) - (*x)[2];
    H[0][0] = -dx / d;
    H[0][1] = -dy / d;
    H[0][2] = 0;
    H[1][0] = dy / (d * d);
    H[1][1] = -dx / (d * d);
    H[1][2] = -1;
    H[0][fpos] = dx / d;
    H[1][fpos] = -dy / (d * d);
    H[0][fpos + 1] = dy / d;
    H[1][fpos + 1] = dx / (d * d);

    v[0] = observations[i][0] - z[0];
    v[1] = observations[i][1] - z[1];
    pi_to_pi(v, 2, 1);

    KF_update(x, P, v, R, H);
  }
  //}
  FINISH_ACTIVITY(ACTIVITY_UPDATE);
}

void KF_update(double (*x)[5], double (*P)[3][3], double v[2], double R[2][2],
               double H[2][5]) {
  double PHt[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  double S[2][2];
  double Si[2][2];
  double Sisym[2][2];
  double W[3][2];
  double temp1[3][2];
  double temp2[3][3];
  double temp3[3][3];
  int i, j, k;
  double temp;
  //#if OMP
  //#pragma omp parallel
  //#endif
  {
    //#if OMP
    //#pragma omp for
    //#endif
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 2; j++) {
        // temp = 0;
        for (k = 0; k < 3; k++) {
          PHt[i][j] += (*P)[i][k] * H[j][k];
        }
        // PHt[i][j] = temp;
      }
    }
    //#if OMP
    //#pragma omp for
    //#endif

    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        temp = 0;
        for (k = 0; k < 3; k++) {
          temp = temp + H[i][k] * PHt[k][j];
        }
        S[i][j] = temp + R[i][j];
      }
    }
    double detS = sqrt(S[0][0] * S[1][1] - S[0][1] * S[1][0]);
    Si[0][0] = S[1][1] / detS;
    Si[0][1] = -S[0][1] / detS;
    Si[1][0] = -S[1][0] / detS;
    Si[1][1] = S[0][0] / detS;

    // make_symmetric(Si, &Sisym);

    multiply2(PHt, Si, &W);

    (*x)[0] = (*x)[0] + W[0][0] * v[0] + W[1][0] * v[0] + W[2][0] * v[0];
    (*x)[1] = (*x)[1] + W[0][1] * v[1] + W[1][1] * v[1] + W[2][1] * v[1];

    multiply2(W, S, &temp1);

#if OMP
#pragma omp for
#endif
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        temp = 0;
        for (k = 0; k < 2; k++) {
          temp = temp + temp1[i][k] * W[j][k];
        }
        temp2[i][j] = temp;
      }
    }

    make_symmetric(temp2, &temp3);

#if OMP
#pragma omp for
#endif
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        (*P)[i][j] = (*P)[i][j] - temp3[i][j];
      }
    }
  }
}

void make_symmetric(double P[3][3], double (*res)[3][3]) {
  {
    int i, j;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        (*res)[i][j] = (P[i][j] + P[j][i]) * 0.5;
      }
    }
  }
}
