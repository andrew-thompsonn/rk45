#include "util.h"


uint8_t parseArguments(int argc, char *argv[], configuration_t *configuration) {

	/* If we don't get the correct number of arguments */
	if (EXPECTED_ARGS != argc)
		return 0;

	/* Retrieve arguments */
	configuration->objective = strtod(argv[1], (char **)NULL);
	configuration->clearance = strtod(argv[2], (char **)NULL);
	configuration->accuracy  = strtod(argv[3], (char **)NULL);
	configuration->startTime = START_TIME;
	configuration->endTime   = END_TIME;
	configuration->timeStep  = TIME_STEP;
	configuration->stateSize = THREE_BODY_STATE_SIZE;
    configuration->loggingEnabled = 0;
	
	sprintf(configuration->fileName, "output/Optimum_%d_%.3f_%.3f", 
					configuration->objective,
					configuration->clearance,
					configuration->accuracy);

	removeDots(configuration->fileName);
	return 1;
}


void removeDots(char string[MAX_FILE_NAME_SIZE]) {
	for (uint8_t index = 0; index < MAX_FILE_NAME_SIZE; index++) {
		char character = string[index];
		if (character == '\0') return;
		if (character == '.') 
			string[index] = 'p';
	}	
}



configuration_t getConfiguration() {

	/* Configuration for the integration */
	configuration_t config = {
		.startTime = 0.0,
		.endTime = 10.0,
		.timeStep = 0.01,
		.stateSize = 12
	};
	return config;
}


void fillInitialConditions(double *stateBuffer, uint8_t size) {

	/* Separated these computation to avoid overflow */
	double numerator = sqrt(G)*MASS_EARTH; 	
	double velocityMoon = numerator/sqrt((MASS_EARTH+MASS_MOON)*DIST_EARTH_MOON);

	state_t state;
	state.xs  = DIST_EARTH_SAT*cosf(THETA0_RAD);
	state.ys  = DIST_EARTH_SAT*sinf(THETA0_RAD);
	state.vxs = VEL_SAT*cosf(THETA0_RAD);
	state.vys = VEL_SAT*sinf(THETA0_RAD);
	//state.xm  = DIST_EARTH_MOON*cosf(THETA_M0_RAD);
	//state.ym  = DIST_EARTH_MOON*sinf(THETA_M0_RAD);
	//state.vxm = -velocityMoon*sinf(THETA_M0_RAD);
	//state.vym = velocityMoon*cosf(THETA_M0_RAD);
	state.xe  = 0;
	state.ye  = 0;
	state.vxe = 0;
	state.vye = 0;
	state.xm  = DIST_EARTH_MOON*cosf(THETA_M0_RAD);
	state.ym  = DIST_EARTH_MOON*sinf(THETA_M0_RAD);
	state.vxm = -velocityMoon*sinf(THETA_M0_RAD);
	state.vym = velocityMoon*cosf(THETA_M0_RAD);

	memcpy(stateBuffer, &state, size*sizeof(double));
}

