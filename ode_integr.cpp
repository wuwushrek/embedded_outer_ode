#include "ode_integr.h"
#include <cstdio>

TM_val::TM_val(AF1 *xn , real t_end , real tn)
{
	this->tn = tn;
	this->t_end = t_end;
	this->xn = xn;

	// Initialize xn interval values
	// for (uint8_t i=0; i<N_STATE ; i++){
	// 	this->xn_int[i] = xn[i].getInterval();
	// }
}

TM_val::~TM_val(){}

void TM_val::buildAndEval()
{
	// Store the ode T_ORDER derivatives
	ode_derivatives(this->xn , this->der);

	// Use fixpoint iteration to get the step size and the remainder
	real curr_tau = fixpoint();

	// Get remainder
	ode_reminder(this->apriori_enclosure , this->rem);

	// Eval the T_ORDER first derivatives and save it inside xn
	real temp_tau = curr_tau; 
	for (uint8_t i=0 ; i< T_ORDER ; i++){
		for (uint8_t j=0 ; j< N_STATE ; j++){
			this->xn[j] += this->der[i + j*T_ORDER] * temp_tau;
		}
		temp_tau *= curr_tau;
	}

	// Add the remainder to the previous sum
	for(uint8_t i=0 ; i< N_STATE ; i++){
		this->xn[i] += this->rem[i] * temp_tau;
	}

	// Update the new time
	this->tn += curr_tau;

	// Update interval value for xn
	for (uint8_t i=0 ; i<N_STATE ; i++){
        //this->xn[i].print_af();
        //printf("Rayon : %f\n", this->xn[i].getRadius());
		this->xn[i].compress_af(NOISE_TOL);
        //  this->xn[i].print_af();
		// this->xn_int[i] = this->xn[i].getInterval();
	}

}

real TM_val::fixpoint()
{
	// calcul de x satisfaisant x0 + [0,tau][f](x) \subseteq x
    AF1 y1[N_STATE];
    AF1 fx0[N_STATE];
    double step = H_MAX;

    int iter;
    Interval widen, coeff;
    
    ode_function(this->xn , fx0);

    for (int i=0; i<N_STATE ; i++){
        y1[i] = this->xn[i] + fx0[i].getInterval();
    }
    
    iter = 1;
    widen = Interval(-1,1);

    uint8_t found_fixpoint = 0;
    while (iter <=1 || (found_fixpoint == 0))
    {
        if (iter > 25)
            coeff = 1;
        else if (iter > 20)
            coeff = 0.1;
        else if (iter > 15)
            coeff = 0.01;
        else if (iter > 10)
            coeff = 0.001;
        else if (iter > 5)
            coeff = 0.0001;
        else if (iter > 2)
            coeff = 0.00001;
        
        if (iter > 2)
        {
            for (int i=0; i<N_STATE ; i++){
                y1[i] = y1[i] + coeff*widen*y1[i].getInterval();
            }
        }

        for(uint8_t i=0 ; i<N_STATE ; i++){
        	this->apriori_enclosure[i] = y1[i];
        }

        ode_function(this->apriori_enclosure,fx0);

        found_fixpoint = 1;

        for (int i=0; i<N_STATE ; i++){
        	y1[i] = this->xn[i] + Interval(0,step) * fx0[i].getInterval();
        	if ( subseteq(y1[i].getInterval() , this->apriori_enclosure[i].getInterval()) == 0){
        		found_fixpoint = 0;
        		i = N_STATE;
        	}
        }
        iter = iter+1;
    }
    // printf("iter = %d , step = %f \n", (int) iter , step);
    return step;
}