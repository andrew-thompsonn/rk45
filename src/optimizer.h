#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "util.h"
#include "integrator.h"


void optimizeDeltaV(configuration_t configuration, double *optdvx, double *optdvy);

double optimizeReturnTime(configuration_t configuration, double *optdvx, double *optdvy);


#endif /* _OPTIMIZER_H_ */
