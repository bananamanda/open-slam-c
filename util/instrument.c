#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef OMP
#define OMP 0
#endif

#if OMP
#include <omp.h>
#else
#include "fake_omp.h"
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#include "cycletimer.h"
#include "instrument.h"

#define MAX_THREAD 64

/* Instrument different sections of program */
static char *activity_name[ACTIVITY_COUNT] = {
    "unknown", "startup",         "predict", "get_observations",
    "update",  "compute_weights", "resample"};

static bool initialized = false;

static bool tracking = false;
static double global_start_time = 0.0;

#define MAXDEPTH 100

static activity_t activity_stack[MAXDEPTH];
static int stack_level = 0;

static double current_start_time = 0.0;
static double accum[MAX_THREAD][ACTIVITY_COUNT];
static double global_accum[ACTIVITY_COUNT];

void track_activity(bool enable) { tracking = enable; }

static void init_instrument() {
  if (!tracking)
    return;
  if (initialized)
    return;
  initialized = true;
  global_start_time = currentSeconds();
  memset(accum, 0, ACTIVITY_COUNT * MAX_THREAD * sizeof(double));
  memset(global_accum, 0, ACTIVITY_COUNT * sizeof(double));
  stack_level = 0;
  activity_stack[stack_level] = ACTIVITY_NONE;
}

void start_activity(activity_t a) {
  if (!tracking)
    return;
  init_instrument();
  int olda = activity_stack[stack_level];
  double new_time = currentSeconds();
  global_accum[olda] += new_time - current_start_time;
  current_start_time = new_time;
  activity_stack[++stack_level] = a;
  if (stack_level >= MAXDEPTH) {
    fprintf(stderr, "Runaway instrumentation activity stack.  Disabling\n");
    tracking = false;
    return;
  }
}

void finish_activity(activity_t a) {
  if (!tracking)
    return;
  init_instrument();
  int olda = activity_stack[stack_level];
  if (a != olda) {
    fprintf(stderr,
            "Warning.  Started activity %s, but now finishing activity %s.  "
            "Disabling\n",
            activity_name[olda], activity_name[a]);
    tracking = false;
    return;
  }
  double new_time = currentSeconds();
  global_accum[olda] += (new_time - current_start_time);
  current_start_time = new_time;
  stack_level--;
  if (stack_level < 0) {
    fprintf(stderr, "Warning, popped off bottom of instrumentation activity "
                    "stack.  Disabling\n");
    tracking = false;
    return;
  }
}

void show_activity(FILE *f, int num_waypoints, int num_landmarks) {
  if (!tracking)
    return;
  init_instrument();
  int a;
  double elapsed = currentSeconds() - global_start_time;
  double unknown = elapsed;
  for (a = 1; a < ACTIVITY_COUNT; a++)
    unknown -= global_accum[a];
  global_accum[0] = unknown;
  fprintf(f, "    %8d waypoints %8d landmarks\n", num_waypoints, num_landmarks);
  for (a = 0; a < ACTIVITY_COUNT; a++) {
    if (global_accum[a] == 0.0)
      continue;
    double ms = global_accum[a] * 1000.0;
    double pct = global_accum[a] / elapsed * 100.0;
    fprintf(f, "    %8d ms    %5.1f %%    %s\n", (int)ms, pct,
            activity_name[a]);
  }
  fprintf(f, "    %8d ms    %5.1f %%    elapsed\n", (int)(elapsed * 1000.0),
          100.0);
}
