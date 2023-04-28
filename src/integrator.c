#include "integrator.h"

/**
 * Log a single state to the file in the format: "time, state[0], state[1], ...state[n] \n"
 */
static void write(FILE *file, double *state, uint8_t stateSize, double time);

/**
 * Scalar multiplication
 */
static void multiplyState(double *state, double coeff, uint8_t length);

/**
 * Print a state
 */
static void printState(double *state);

/**
 * Add two states
 */
static void incrementState(double *state, double *increment, uint8_t length);

/**
 * Allocate memory for all buffers used in integration
 */
static void allocate(configuration_t config, double **a, double **b, double **c, 
                                             double **d, double **e, double **f);
/**
 * Free all memory previously allocated for integration
 */
static void freeAll();

/**
 * Construct the k1 state
 */
static void constructK1(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config);
/**
 * Construct the k2 state
 */
static void constructK2(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config);
/**
 * Construct the k3 state
 */
static void constructK3(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config);
/**
 * Construct the k4 state
 */
static void constructK4(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config);
/**
 * Construct the k5 state
 */
static void constructK5(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config);
/**
 * Construct the k6 state
 */
static void constructK6(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config);
/**
 * Construct the 4th order state 'A'
 */
static void constructA(double *A, double *state, configuration_t config);

/**
 * Construct the 5th order state 'B'
 */
static void constructB(double *B, double *state, configuration_t config);

/**
 * Global buffers private to this file (true values)
 */
static double *k1,  *k2,  *k3,  *k4,  *k5, *k6;

/**
 * Global buffers private to this file (copies of true values)
 */
static double *k1c, *k2c, *k3c, *k4c, *k5c, *k6c;


uint8_t rk45(uint8_t (*function)(double time, double *stateVector),
			   double *initialConditions, configuration_t config, double *stopTime) {

    /* Initialize k buffers and their copies */
    allocate(config, &k1, &k2, &k3, &k4, &k5, &k6);
    allocate(config, &k1c, &k2c, &k3c, &k4c, &k5c, &k6c);

    /* Initialize the current state from the initial conditions */
    double currentState[config.stateSize];
    memcpy(currentState, initialConditions, config.stateSize*sizeof(double));
	
    /* Other buffers, and return code for integration termination */
	double stateA[config.stateSize], stateB[config.stateSize], difference[config.stateSize];
	
    /* If logging is enabled, open the output file  */
    FILE *file;
	if (config.loggingEnabled) 
        file = fopen(config.fileName, "w");
	
    /* The integration will always start at time t = 0 */	
	double time = 0;
	uint8_t returnCode = 0;

	/* While the absolute return code does not indicate a collision */	
	while (returnCode == 0 && time <= config.endTime) {

        /* Construct states k1 through k6 */
        constructK1(function, time, currentState, config);
        constructK2(function, time, currentState, config);
        constructK3(function, time, currentState, config);
        constructK4(function, time, currentState, config);
        constructK5(function, time, currentState, config);
        constructK6(function, time, currentState, config);
		
	    /* Construct fourth and fifth order states */	
        constructA(stateA, currentState, config);    
        constructB(stateB, currentState, config);    

		/* Copy the 4th order solution, and initialize the difference buffer */
		double possibleSolution[config.stateSize];
		memcpy(possibleSolution, stateA, config.stateSize*sizeof(double));
		memcpy(difference, stateB, config.stateSize*sizeof(double));

        /* Compute the difference between the 4th and 5th order solution */
		multiplyState(stateA, -1.0f, config.stateSize);	
		incrementState(difference, stateA, config.stateSize);

		/* Get the magnitude of the difference */
		double sum = 0;
		for (uint8_t i = 0; i < config.stateSize; i++) 
			sum += powf(difference[i], 2.0f);
		double norm = sqrtf(sum);

		/* Compute delta */
		double delta = DELTA_COEF*powf((RK45_TOL/norm), 1.0f/4.0f);

		/* If the accuracy is acceptable, */
		if (norm/config.timeStep <= RK45_TOL) {

			/* Increment the current time, and copy the new state */
			time += config.timeStep;
			memcpy(currentState, possibleSolution, config.stateSize*sizeof(double));
		    
            /* Write the resulting state to the output file */
		    if (config.loggingEnabled) 
                write(file, currentState, config.stateSize, time);
        
            /* Check for a collision */ // TODO: This function should be a parameter to rk45
            returnCode = checkCollisionArray(currentState);
            if (returnCode != 0) break;
		} 
		/* Increment time step */
		config.timeStep *= delta;
	}
    /* Close file, free memory, etc. */
	if (config.loggingEnabled) fclose(file);
	(*stopTime) = time;
    freeAll();
	return returnCode;
}

void constructK1(uint8_t (*func)(double t, double *s), double t, double *s, 
                 configuration_t config) {

	/* Initialize first true buffer with the state */
	memcpy(k1, s, config.stateSize*sizeof(double));
	(func)(t, k1);	

	multiplyState(k1, config.timeStep, config.stateSize);

	/* Reset copies, initialize next buffer with k1 */
	memcpy(k1c, k1, config.stateSize*sizeof(double));
	memcpy(k2,  k1, config.stateSize*sizeof(double));
}


void constructK2(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config) {

	multiplyState(k2, K2_K1_COEF, config.stateSize);
	incrementState(k2, s, config.stateSize);

	(func)(t + config.timeStep*K2_H_COEF, k2);
	multiplyState(k2, config.timeStep, config.stateSize);
	
	/* Reset copies, initialize next buffer with k2 */
	memcpy(k2c, k2, config.stateSize*sizeof(double));
	memcpy(k3,  k2, config.stateSize*sizeof(double));
}

void constructK3(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config) {

	multiplyState(k1c, K3_K1_COEF, config.stateSize);
	multiplyState(k3,  K3_K2_COEF, config.stateSize);

	incrementState(k3, s, config.stateSize);
	incrementState(k3, k1c, config.stateSize);

	(func)(t + config.timeStep*K3_H_COEF, k3);
	multiplyState(k3, config.timeStep, config.stateSize);

	/* Reset copies, initialize next buffer with k3 */
	memcpy(k3c, k3, config.stateSize*sizeof(double));
	memcpy(k1c, k1, config.stateSize*sizeof(double));
	memcpy(k4,  k3, config.stateSize*sizeof(double));

}

void constructK4(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config) {

	multiplyState(k1c, K4_K1_COEF, config.stateSize);	
	multiplyState(k2c, K4_K2_COEF, config.stateSize);
	multiplyState(k4,  K4_K3_COEF, config.stateSize);

	incrementState(k4, s, config.stateSize);
	incrementState(k4, k1c, config.stateSize);
	incrementState(k4, k2c, config.stateSize);

	(func)(t + config.timeStep*K4_H_COEF, k4);
	multiplyState(k4, config.timeStep, config.stateSize);

	/* Reset copies, initialize next buffer with k4 */
	memcpy(k4c, k4, config.stateSize*sizeof(double));
	memcpy(k1c, k1, config.stateSize*sizeof(double));
	memcpy(k2c, k2, config.stateSize*sizeof(double));
	memcpy(k5,  k4, config.stateSize*sizeof(double));
}

void constructK5(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config) {

	multiplyState(k1c, K5_K1_COEF, config.stateSize);
	multiplyState(k2c, K5_K2_COEF, config.stateSize);
	multiplyState(k3c, K5_K3_COEF, config.stateSize);
	multiplyState(k5,  K5_K4_COEF, config.stateSize);

	incrementState(k5, s, config.stateSize);
	incrementState(k5, k1c, config.stateSize);
	incrementState(k5, k2c, config.stateSize);
	incrementState(k5, k3c, config.stateSize);

	(func)(t + config.timeStep*K5_H_COEF, k5);
	multiplyState(k5, config.timeStep, config.stateSize);
	
	/* Reset copies, initialize next buffer with k5 */	
	memcpy(k5c, k5, config.stateSize*sizeof(double));
	memcpy(k1c, k1, config.stateSize*sizeof(double));
	memcpy(k2c, k2, config.stateSize*sizeof(double));
	memcpy(k3c, k3, config.stateSize*sizeof(double));
	memcpy(k6,  k5, config.stateSize*sizeof(double));
}

void constructK6(uint8_t (*func)(double t, double *s), double t, double *s, 
                        configuration_t config) {
	multiplyState(k1c, K5_K1_COEF, config.stateSize);
	multiplyState(k2c, K5_K1_COEF, config.stateSize);
	multiplyState(k3c, K5_K1_COEF, config.stateSize);
	multiplyState(k4c, K5_K1_COEF, config.stateSize);
	multiplyState(k6,  K5_K1_COEF, config.stateSize);

	incrementState(k6, s, config.stateSize);
	incrementState(k6, k1c, config.stateSize);
	incrementState(k6, k2c, config.stateSize);
	incrementState(k6, k3c, config.stateSize);
	incrementState(k6, k4c, config.stateSize);

	(func)(t + config.timeStep*K6_H_COEF, k6);
	multiplyState(k6, config.timeStep, config.stateSize);

	/* Reset copies */
	memcpy(k6c, k6, config.stateSize*sizeof(double));
	memcpy(k1c, k1, config.stateSize*sizeof(double));
	memcpy(k2c, k2, config.stateSize*sizeof(double));
	memcpy(k3c, k3, config.stateSize*sizeof(double));
	memcpy(k4c, k4, config.stateSize*sizeof(double));
}

void constructB(double *B, double *state, configuration_t config) {

	memcpy(B, state, config.stateSize*sizeof(double));

    multiplyState(k1c, STATE_B_K1_COEF, config.stateSize);
    multiplyState(k3c, STATE_B_K3_COEF, config.stateSize);
    multiplyState(k4c, STATE_B_K4_COEF, config.stateSize);
    multiplyState(k5c, STATE_B_K5_COEF, config.stateSize);
    multiplyState(k6c, STATE_B_K6_COEF, config.stateSize);

    incrementState(B, k1c, config.stateSize);
    incrementState(B, k3c, config.stateSize);
    incrementState(B, k4c, config.stateSize);
    incrementState(B, k5c, config.stateSize);
    incrementState(B, k6c, config.stateSize);
}

void constructA(double *A, double *state, configuration_t config) {

    memcpy(A, state, config.stateSize*sizeof(double));
	
    multiplyState(k1c, STATE_A_K1_COEF, config.stateSize);
    multiplyState(k3c, STATE_A_K3_COEF, config.stateSize);
    multiplyState(k4c, STATE_A_K4_COEF, config.stateSize);
    multiplyState(k5c, STATE_A_K5_COEF, config.stateSize);

    incrementState(A, k1c, config.stateSize);
    incrementState(A, k3c, config.stateSize);
    incrementState(A, k4c, config.stateSize);
    incrementState(A, k5c, config.stateSize);

    /* Reset copies */
    memcpy(k1c, k1, config.stateSize*sizeof(double));
    memcpy(k2c, k2, config.stateSize*sizeof(double));
    memcpy(k3c, k3, config.stateSize*sizeof(double));
    memcpy(k3c, k4, config.stateSize*sizeof(double));
}

uint8_t euler(uint8_t (*function)(double time, double *stateVector),
			   double *initialConditions, configuration_t config, double *stopTime) {

	/* Declare buffers for state & state derivative */
	double currentState[config.stateSize];
	double stateDerivative[config.stateSize];

	/* Fill the buffers with the initial conditions */	
	memcpy(currentState, initialConditions, config.stateSize*sizeof(double));
	memcpy(stateDerivative, initialConditions, config.stateSize*sizeof(double));

	/* Open the output file and write the initial state */
    FILE *file;
	if (config.loggingEnabled) {
        file = fopen(config.fileName, "w");
	    write(file, initialConditions, config.stateSize, (double)0);
    }
	/* For every time step */
	for (double currentTime = config.startTime; currentTime < config.endTime; 
					currentTime += config.timeStep) {

		/* Fill derivative buffer with current state */
		memcpy(stateDerivative, currentState, config.stateSize*sizeof(double));

		/* Get the state, and check if a terminal condition occurred */
		uint8_t returnCode = (function)(currentTime, stateDerivative);
		if (returnCode != 0) {
			if (config.loggingEnabled)  fclose(file);
            *stopTime = currentTime;
			return returnCode;
		}
		/* Compute the derivative, multiply by time */
		multiplyState(stateDerivative, config.timeStep, config.stateSize);

		/* Increment the current state by the derivative multiplied by time */
		incrementState(currentState, stateDerivative, config.stateSize);

		/* Write the resulting state to the ouptut file */
		if (config.loggingEnabled) write(file, currentState, config.stateSize, currentTime);
	}
    /* Close file and return success */
	if (config.loggingEnabled) fclose(file);
    *stopTime = config.endTime;
	return 0;
}

void write(FILE *file, double *state, uint8_t stateSize, double time) {

	/* First column is time */
	fprintf(file, "%f ", time);

	/* Write each element of the state to the file, followed by new line */	
	for (uint8_t index = 0; index < stateSize; index++) 
		fprintf(file, "%f ", state[index]);
	fprintf(file, "\n");	
}

void multiplyState(double *state, double coeff, uint8_t length) {

	/* Multiply each element of the state by the coefficient */
	for (uint8_t index = 0; index < length; index++) 
		state[index] *= coeff;
}

void incrementState(double *state, double *increment, uint8_t length) {

	/* Elementwise addition for two states */
	for (uint8_t index = 0; index < length; index++)
		state[index] += increment[index];	
}	

void allocate(configuration_t config, double **a, double **b, double **c, 
                                             double **d, double **e, double **f) {
    /* Allocate memory for each pointer (a-f) */
    (*a) = (double *)calloc(config.stateSize, sizeof(double));
    (*b) = (double *)calloc(config.stateSize, sizeof(double));
    (*c) = (double *)calloc(config.stateSize, sizeof(double));
    (*d) = (double *)calloc(config.stateSize, sizeof(double));
    (*e) = (double *)calloc(config.stateSize, sizeof(double));
    (*f) = (double *)calloc(config.stateSize, sizeof(double));
}

void freeAll() {

    /* Free the memory of all static arrays */
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
    free(k1c);
    free(k2c);
    free(k3c);
    free(k4c);
    free(k5c);
    free(k6c);
}

void printState(double *state) {
	printf("\tSPACECRAFT POSITION X:\t\t%.6f\n", state[0]);
	printf("\tSPACECRAFT POSITION Y:\t\t%.6f\n", state[1]);
	printf("\tSPACECRAFT VELOCITY X:\t\t%.6f\n", state[2]);
	printf("\tSPACECRAFT VELOCITY Y:\t\t%.6f\n", state[3]); printf("\n");
	printf("\tEARTH POSITION X:     \t\t%.6f\n", state[4]);
	printf("\tEARTH POSITION Y:     \t\t%.6f\n", state[5]);
	printf("\tEARTH VELOCITY X:     \t\t%.6f\n", state[6]);
	printf("\tEARTH VELOCITY Y:     \t\t%.6f\n", state[7]); printf("\n");
	printf("\tMOON POSITION X:      \t\t%.6f\n", state[8]);
	printf("\tMOON POSITION Y:      \t\t%.6f\n", state[9]);
	printf("\tMOON VELOCITY X:      \t\t%.6f\n", state[10]);
	printf("\tMOON VELOCITY Y:      \t\t%.6f\n", state[11]);
	printf("\n");
}

