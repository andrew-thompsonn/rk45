#ifndef _EQUATIONS_H_
#define _EQUATIONS_H_

#include <stdint.h>
#include <string.h>

#include <math.h>
#include "definitions.h"

#define TRUE 	(1)
#define FALSE 	(0)
#define THREE_BODY_STATE_SIZE 	(12)

#define RESULT_COLLISION_EARTH 	(1)
#define RESULT_COLLISION_MOON   (2)
#define RESULT_ESCAPE 			(3)

/* A struct to represent a state */
typedef struct {

	/* Spacecraft position & velocity */
	double xs;
	double ys;
	double vxs;
	double vys;

	/* Earth position & velocity */
	double xe;
	double ye;
	double vxe;
	double vye;

	/* Moon position & velocity */
	double xm;
	double ym;
	double vxm;
	double vym;

} state_t;

/**
 * Takes a 12x1 array representing the state of the system and fills that array
 * with the derivative of the state at the specified time. Return value indicates
 * whether integration should terminate early.
 */
uint8_t equations(double time, double *stateIn);

/* Compute the force between two bodies */
double force(double mass1, double mass2, double X1, double X2, double Y1, double Y2, double dir);

/* Compute 2D distance between two objects */
double distance(double x1, double y1, double x2, double y2);

/* Check if a collision has occurred based on the current state */
uint8_t checkCollision(state_t state);

/* Differentiate the state */
void differentiate(state_t *state, double axSat, double aySat, double axMoon, double ayMoon);

/* Set the clearance variable */
void setClearance(double clearanceIn); 

/* Check if a collision has occurred from a state array */
uint8_t checkCollisionArray(double *stateIn);

#endif /* _EQUATIONS_H_ */

