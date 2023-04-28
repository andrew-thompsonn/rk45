/* ASEN 4057 Final Project */

#include <stdlib.h>
#include <time.h>

#include "integrator.h"
#include "equations.h"
#include "util.h"
#include "optimizer.h"

//#define _DEBUG
#define DEBUG_DVX (-1)
#define DEBUG_DVY (49)

int main(int argc, char *argv[]) {

    clock_t start = clock();

    /* Parse command line arguments */
	configuration_t configuration;
	if (!parseArguments(argc, argv, &configuration))
		return EXIT_FAILURE;

    /* Set the clearance of the spacecraft and moon */
	setClearance((double)configuration.clearance);

    /* Initialize impulses and return time */
    double optdvx, optdvy, bestTime;
#ifndef _DEBUG 
    /* Optimized the desired objective */
    switch (configuration.objective) {
        case OBJECTIVE_1:
            optimizeDeltaV(configuration, &optdvx, &optdvy);
            break;
        case OBJECTIVE_2:
            bestTime = optimizeReturnTime(configuration, &optdvx, &optdvy);
            break;
    }
#else 
        optdvx = DEBUG_DVX;
        optdvy = DEBUG_DVY;
#endif 
    /* Create function pointer, initial conditions buffer */
	uint8_t (*diffEquation)(double time, double *stateVector) = &equations;
    double initialConditions[configuration.stateSize];  
    
    /* Get initial conditions for optimal configuration */
	fillInitialConditions(initialConditions, configuration.stateSize);
    initialConditions[2] += optdvx;
    initialConditions[3] += optdvy;

    /* Integrate optimal initial conditions with logging enabled */   
    double time;
    configuration.loggingEnabled = 1;
	rk45(diffEquation, initialConditions, configuration, &time);

    printf("\n\tSolution: (dvx, dvy) = (%.2f, %.2f)\n", optdvx, optdvy);
    printf("\n\t* Output written to: %s\n\n", configuration.fileName);

    clock_t end = clock();
    double runTime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("Run time: %.3f seconds\n\n", runTime);

	return EXIT_SUCCESS;
}
