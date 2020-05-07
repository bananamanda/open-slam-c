#ifndef INSTRUMENT_H
#include <stdbool.h>

#ifndef TRACK
#define TRACK 1
#endif

/* Keep track of how time gets used */

/* Categories of activities */

typedef enum { ACTIVITY_NONE, ACTIVITY_STARTUP, ACTIVITY_PREDICT, ACTIVITY_OBS, ACTIVITY_UPDATE, ACTIVITY_WEIGHTS, ACTIVITY_RESAMPLE, ACTIVITY_COUNT} activity_t;

void track_activity(bool enable);

void start_activity(activity_t a);
void finish_activity(activity_t a);
void show_activity(FILE *f, int num_waypoints, int num_landmarks);

#if TRACK
#define TRACK_ACTIVITY(e) track_activity(e)
#define START_ACTIVITY(a) start_activity(a)
#define FINISH_ACTIVITY(a) finish_activity(a)
#define SHOW_ACTIVITY(f,nn,ne) show_activity(f,nn,ne)
#else
#define TRACK_ACTIVITY(e)  /* Optimized out */
#define START_ACTIVITY(a)   /* Optimized out */
#define FINISH_ACTIVITY(a)  /* Optimized out */
#define SHOW_ACTIVITY(f)  /* Optimized out */
#endif

#define INSTRUMENT_H
#endif
