#ifndef __ODE_INTEGR_H__
#define __ODE_INTEGR_H__

#include "config_outer.h"
#include "ode_dyn.h"

class TM_val
{

public:
	real tn;
	real t_end;
	AF1  *xn; // xn[N_STATES];
	// Interval xn_int[N_STATE];

public:
	// local variables are stored at class creation
	AF1 der[T_ORDER*N_STATE];
	AF1 rem[N_STATE];
	AF1 apriori_enclosure[N_STATE];

	// Interval der_int[T_ORDER*N_STATE];
	// Interval rem_int[N_STATE];
	// Interval apriori_int[N_STATE];

public:
	TM_val(AF1 *x0 , real t_end , real tn = 0.0f);
	~TM_val();

	void buildAndEval();

	real fixpoint();

};

#endif // __ODE_INTEGR_H__