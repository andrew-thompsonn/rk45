#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <stdint.h>

#define MAX_FILE_NAME_SIZE      (40)
#define START_TIME              (0)
#define END_TIME                (1E8)
#define TIME_STEP               (5)
#define RK45_TOL                (2E3)
#define RK45_MIN_STEP           (1)

/**
 * Parameters for integration
 */
typedef struct {

    /* Timing */
	double startTime;  
	double endTime;    
	double timeStep;   
    
    /* Arguments */
	uint8_t stateSize;
    
    /* Configuration */
	uint8_t objective;
	double clearance;
	double accuracy;

    /* File output */
    uint8_t loggingEnabled;
	char fileName[MAX_FILE_NAME_SIZE];

} configuration_t;

#endif /* _CONFIGURATION_H_ */
