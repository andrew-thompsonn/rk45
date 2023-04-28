#include "equations.h"
#include <stdio.h>

static double clearance;

uint8_t equations(double time, double *stateBuffer) {

	/* Copy the state into a struct (more readable) */
	state_t state;
	memcpy(&state, stateBuffer, sizeof(state_t));

	/* Check for a collision */	
	uint8_t status = checkCollision(state);

	/* Forces acting on the spacecraft */
	double fxMoonOnSat  = force(MASS_MOON, MASS_SAT, state.xs,  state.xm, state.ys, state.ym, 1);
	double fyMoonOnSat  = force(MASS_MOON, MASS_SAT, state.xs,  state.xm, state.ys, state.ym, 2);
	double fxEarthOnSat = force(MASS_EARTH, MASS_SAT, state.xs, state.xe, state.ys, state.ye, 1);
	double fyEarthOnSat = force(MASS_EARTH, MASS_SAT, state.xs,  state.xe, state.ys, state.ye, 2);

	/* Forces acting on the moon */
	double fxEarthOnMoon = force(MASS_EARTH, MASS_MOON, state.xm, state.xe, state.ym, state.ye, 1);
	double fyEarthOnMoon = force(MASS_EARTH, MASS_MOON, state.xm, state.xe, state.ym, state.ye, 2);
	double fxSatOnMoon   = -fxMoonOnSat;
	double fySatOnMoon   = -fyMoonOnSat;

	/* Accelerations */
	double axSat  = (fxEarthOnSat + fxMoonOnSat)/MASS_SAT;
	double aySat  = (fyEarthOnSat + fyMoonOnSat)/MASS_SAT;
	double axMoon = (fxSatOnMoon + fxEarthOnMoon)/MASS_MOON;
	double ayMoon = (fySatOnMoon + fyEarthOnMoon)/MASS_MOON;

	/* Differentiate the state and copy to the buffer */
	differentiate(&state, axSat, aySat, axMoon, ayMoon);
	memcpy(stateBuffer, &state, sizeof(state));

	if (status) return status;

	/* If we are here, then a collision has not occurred */
	return 0;
}

uint8_t checkCollisionArray(double *stateIn) {

    state_t state;
    memcpy(&state, stateIn, sizeof(state_t));
    return checkCollision(state);
}


uint8_t checkCollision(state_t state) {
	
	/* Compute distances */
	double distanceEarthMoon = distance(state.xe, state.ye, state.xm, state.ym);
	double distanceEarthSat  = distance(state.xe, state.ye, state.xs, state.ys);
	double distanceMoonSat 	 = distance(state.xm, state.ym, state.xs, state.ys);

	/* Check if a collision occurred with the moon, earth, or spacecraft has escaped */
	uint8_t collidedWithMoon  = (distanceMoonSat < RADIUS_MOON + clearance);
	uint8_t collidedWithEarth = (distanceEarthSat < RADIUS_EARTH);
	uint8_t escapedOrbit      = (distanceEarthSat > 2*distanceEarthMoon);

	/* Integration should terminate if any of these conditions occur */
	if (collidedWithMoon)
		return RESULT_COLLISION_MOON;
	if (collidedWithEarth)
		return RESULT_COLLISION_EARTH;
	if (escapedOrbit)
		return RESULT_ESCAPE;
	return 0;
}


double distance(double x1, double y1, double x2, double y2) {

	/* Return the scalar distance */
	return sqrt( powf((x2-x1), 2) + powf((y2-y1), 2) );
}


double force(double mass1, double mass2, double X1, double X2, double Y1, double Y2, double dir) {

	double d = sqrt( powf(X2-X1, 2) + powf(Y2-Y1, 2) );
	if (dir == 1){
		double result = G*mass1*mass2*(X2 - X1)/(powf(d, 3));
        return result;
	}
	else{
		double result = G*mass1*mass2*(Y2 - Y1)/(powf(d, 3));
        return result;
		
	}
	/* Return the scalar force acting between the two bodies */
}



void differentiate(state_t *state, double axSat, double aySat, 
				   			      double axMoon, double ayMoon) {
	/* Spacecraft */
	state->xs  = state->vxs;
	state->ys  = state->vys;
	state->vxs = axSat;
	state->vys = aySat;

	/* Earth */
	state->xe  = 0;
	state->ye  = 0;
	state->vxe = 0;
	state->vye = 0;

	/* Moon */
	state->xm  = state->vxm;
	state->ym  = state->vym;
	state->vxm = axMoon;
	state->vym = ayMoon;
}


void setClearance(double clearanceIn) {
	clearance = clearanceIn;
}

