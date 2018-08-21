#ifndef __CONFIG_OUTER_H__
#define __CONFIG_OUTER_H__

#ifdef USE_MAF1
	#include "aa.h"
	typedef MAF1 AF1;
#elif USE_MAF2
	#include "aa_mod.h"
	typedef MAF2 AF1;
#else
	#include "aa_mod2.h"
	typedef DAF AF1;
#endif

#ifndef N_STATE
	#define N_STATE 2
#endif

#ifndef T_ORDER
	#define T_ORDER 2
#endif

#ifndef TIME_STEP_FIXED
	#define TIME_STEP_FIXED 1
#endif

#define MAX_WIDEN_FACTOR	1.0
#define MIN_WIDEN_FACTOR	0.001

#define WIDEN_MULT			1.2
#define STEP_MULT			0.5 

#define H_MAX 0.1
#define H_MIN 0.0001

#define NOISE_TOL 1e-1

#endif // __CONFIG_OUTER_H__