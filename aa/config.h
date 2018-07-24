#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdint.h>
#include <cmath>
#include <cstdio>

// Define the number of noise symbol if not already done during compilation process
#ifndef N_NOISE
	#define N_NOISE 10
#endif

// Use double type if it is specified else use float 
#ifdef USE_DOUBLE
	typedef real double;
#else
	typedef real float;
#endif

// Define VERBOSE for compiling with std::cout compatibilty

// Define T_MULT for using trivial multiplication of affine forms instead of standard one

#endif // __CONFIG_H__