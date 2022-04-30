#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "configuration.h"

/* Runge Kutta 45 integration */
uint8_t rk45(uint8_t (*function)(double time, double *stateVector),
		double *initialConditions, configuration_t config, double *stopTime);

#endif /* _INTEGRATOR_H_ */
