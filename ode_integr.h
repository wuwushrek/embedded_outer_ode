#ifndef __ODE_INTEGR_H__
#define __ODE_INTEGR_H__

#include "aa.h"
#include "ode_fun.h"

#ifndef N_STATES
	#define N_STATES 2
#endif

#ifndef T_ORDER
	#define T_ORDER 3
#endif

#define MAX_WIDEN_FACTOR	2
#define WIDEN_DIVIDOR		0.6
#define STEP_DIVIDOR		0.6 

#define H_MAX 0.3
#define H_MIN 0.001

class TM_val
{

private:
	real tn;
	real tau_init;
	AF1  *xn; // xn[N_STATES];

	// local variables are stored at class creation
	AF1 der[T_ORDER*N_STATE];
	AF1 rem[N_STATE];
	AF1 apriori_enclosure[N_STATE];

	Interval der_int[T_ORDER*N_STATE];
	Interval xn_int[N_STATE];

public:
	TM_val(AF *x0 , real tau_init = H_MAX , real tn = 0.0f);
	~TM_val();

	void buildAndEval();

	real fixpoint();

};

#endif // __ODE_INTEGR_H__