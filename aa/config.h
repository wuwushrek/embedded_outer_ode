#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdint.h>
#include <cmath>
#include <cstdio>

// check if compilation is done for microcontroller -> chose your custom assert function
#ifndef USE_MCU
	#include <assert.h>
	#define assert_af(e) assert(e)
#else
	#include "cfassert.h"
	#define assert_af(e) ASSERT(e)
#endif

// Define the number of noise symbol if not already done during compilation process
#ifndef N_NOISE
	#define N_NOISE 1
#endif

// Use double type if it is specified else use float 
#ifdef USE_DOUBLE
	typedef double real;
#else
	typedef float real;
#endif

// Define FAST_MULT for using trivial multiplication of affine forms instead of a more precized one

// Define FAST_DIV for using trivial division of affine forms instead of a more precized one

#define MIN_RAD 1e-8

#define EPS_ZERO 1e-10

#ifndef M_PI
#define M_PI (4*atan(1.0))
#endif 

#endif // __CONFIG_H__