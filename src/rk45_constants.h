#ifndef _RK45_CONSTANTS_H_
#define _RK45_CONSTANTS_H_

#define K2_H_COEF 		1.0f/4.0f
#define K2_K1_COEF 		1.0f/4.0f

#define K3_H_COEF 		3.0f/8.0f
#define K3_K1_COEF 		3.0f/32.0f
#define K3_K2_COEF 		9.0f/32.0f

#define K4_H_COEF 		12.0f/13.0f
#define K4_K1_COEF 		1932.0f/2197.0f
#define K4_K2_COEF 		-7200.0f/2197.0f
#define K4_K3_COEF 		7296.0f/2197.0f

#define K5_H_COEF 		1.0f
#define K5_K1_COEF 		439.0f/216.0f
#define K5_K2_COEF 		-8.0f
#define K5_K3_COEF 		8680.0f/513.0f
#define K5_K4_COEF 		-845.0f/4104.0f

#define K6_H_COEF 		1.0f/2.0f
#define K6_K1_COEF 		-8.0f/27.0f
#define K6_K2_COEF 		2.0f
#define K6_K3_COEF 		-3544.0f/2565.0f
#define K6_K4_COEF 		1859.0f/4104.0f
#define K6_K5_COEF 		-11.0f/40.0f

#define STATE_A_K1_COEF		25.0f/216.0f
#define STATE_A_K3_COEF 	1408.0f/2565.0f
#define STATE_A_K4_COEF 	2197.0f/4104.0f	
#define STATE_A_K5_COEF 	-1.0f/5.0f

#define STATE_B_K1_COEF 	16.0f/135.0f
#define STATE_B_K3_COEF 	6656.0f/12825.0f
#define STATE_B_K4_COEF 	28561.0f/56430.0f
#define STATE_B_K5_COEF 	-9.0f/50.0f
#define STATE_B_K6_COEF 	2.0f/55.0f

#define DELTA_COEF 		0.84f

#endif /* __RK45_CONSTANTS_H_ */


