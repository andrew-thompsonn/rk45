#include "optimizer.h"


void optimizeDeltaV(configuration_t configuration, double *optdvx, double *optdvy) {

	uint8_t (*diffEquation)(double time, double *stateVector) = &equations;
	double initialConditions[configuration.stateSize];
	fillInitialConditions(initialConditions, configuration.stateSize);

	double dv = 100000;
	*optdvx = 0;
	*optdvy = 0;
    double stopTime;

    printf("\nPerforming grid search for minimal delta V...\n");

	for (double dvx = -100; dvx<= 100; dvx+=configuration.accuracy){
		for (double dvy = -100; dvy<= 100; dvy+=configuration.accuracy){
            if (dvx == 0 || dvy == 0) continue;

            printf("\tTesting %.1f, %.1f\n", dvx, dvy); 
			initialConditions[2]+=dvx;
			initialConditions[3]+=dvy;

            //uint8_t result = rk45(diffEquation, initialConditions, configuration, &stopTime);
            uint8_t result = euler(diffEquation, initialConditions, configuration, &stopTime);
            if (RESULT_COLLISION_EARTH == result) {
				if (sqrt( powf(dvx, 2) + powf(dvy, 2) ) < dv){
					dv = sqrt( powf(dvx, 2) + powf(dvy, 2));
					*optdvx = dvx;
					*optdvy = dvy;
				}
            }
			initialConditions[2]-=dvx;
			initialConditions[3]-=dvy;
		}
	}
}


double optimizeReturnTime(configuration_t configuration, double *optdvx, double *optdvy) {

    /* Setup the function pointer to be integrated, fill the initial state */
	uint8_t (*diffEquation)(double time, double *stateVector) = &equations;
	double initialConditions[configuration.stateSize];
	fillInitialConditions(initialConditions, configuration.stateSize);

    /* Initialize values for the delta V, and stop/start time */
	*optdvx = 0;
	*optdvy = 0;
    double stopTime = configuration.endTime;
    double bestTime = configuration.endTime;

    for (double dvx = -100; dvx < 100; dvx+=configuration.accuracy) {
        for (double dvy = -100; dvy < 100; dvy+=configuration.accuracy) {
            if (dvx == 0 || dvy == 0) continue;
            printf("%.1f, %.1f\n", dvx, dvy);
			initialConditions[2]+=dvx;
			initialConditions[3]+=dvy;

            uint8_t result = rk45(diffEquation, initialConditions, configuration, &stopTime);
            //uint8_t result = euler(diffEquation, initialConditions, configuration, &stopTime);

            if (RESULT_COLLISION_EARTH == result) {
                if (stopTime < bestTime) {
                    bestTime = stopTime;
					*optdvx = dvx;
					*optdvy = dvy;
                }
            }
		    initialConditions[2]-=dvx;
			initialConditions[3]-=dvy;
        }
    }
    return bestTime;
}


