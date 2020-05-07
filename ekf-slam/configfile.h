/**
 * Configuration file
 * Permits various adjustments to parameters of the FastSLAM algorithm.
 **/

/* control parameters */
//#define V 3.0 		// m/s
#define MAXG 		30*(3.14 / 180)	// radians
#define RATEG		20*(3.14 / 180)	// rad/s
#define WHEELBASE 	4		// meters
#define DT_CONTROLS 	0.025		// seconds
#define DT_OBSERVE 8*DT_CONTROLS
#define MAX_RANGE 30.0

/* control noises */
#define SIGMAV 		0.3		// m/s
#define SIGMAG 		3.0*(3.14/180)	// radians
double Q[2][2] = { { SIGMAV * SIGMAV, 0 }, { 0, SIGMAG * SIGMAG } };

/* observation parameters */
#define SIGMAR		0.1 		// meters
#define SIGMAB		1.0*(3.14/180)	// radians
double R[2][2] = { { SIGMAR * SIGMAR, 0 }, { 0, SIGMAB * SIGMAB } };

/* waypoint proximity */
#define AT_WAYPOINT	1.0 		// meters
#define NUMBER_LOOPS 	0			// num loops

/* resampling */
#define NPARTICLES 	100
#define NEFFECTIVE 	0.75 * NPARTICLES

#define SWITCH_CONTROL_NOISE 	0
#define SWITCH_SENSOR_NOISE	 0
#define SWITCH_INFLATE_NOISE	0
#define SWITCH_SAMPLE_PROPOSAL	1
#define	SWITCH_HEADING_KNOWN	0
#define SWITCH_RESAMPLE		1
#define SWITCH_PROFILE		1
