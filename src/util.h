#ifndef _UTIL_H_
#define _UTIL_H_

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "configuration.h"
#include "equations.h"

#define EXPECTED_ARGS 	(4)

/* Represents all the arguments to the program */
typedef struct {
	uint8_t objective;
	double clearance;
	double accuracy;
} args_t;

/* Parse command line arguments */
uint8_t parseArguments(int argc, char *argv[], configuration_t *configuration);

/* Retrieve an integration configuration */
configuration_t getConfiguration();

/* Fill the initial conditions for our specific scenario */
void fillInitialConditions(double *state, uint8_t size);

void removeDots(char string[MAX_FILE_NAME_SIZE]);

#endif /* _UTIL_H_ */

