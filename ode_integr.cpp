#include "ode_integr.h"

TM_val::TM_val(AF1 *xn , real tau_init , real tn)
{
	this->tn = tn;
	this->tau_init = tau_init;
	this->xn = xn;

	// Initialize xn interval values
	for (uint8_t i=0; i<N_STATE ; i++){
		this->xn_int[i] = xn[i].getInterval();
	}
}

TM_val::~TM_val(){}

void TM_val::buildAndEval()
{
	// Store the ode T_ORDER derivatives
	ode_function(this->xn , this->der);

	// Use fixpoint iteration to get the step size and the remainder
	real curr_tau = fixpoint();

	// Eval the T_ORDER first derivatives and save it inside xn
	real temp_tau = curr_tau; 
	for (uint8_t i=0 ; i< T_ORDER ; i++){
		for (uint8_t j=0 ; j< N_STATE ; j++){
			this->xn[j] += this->der[j + i*N_STATE] * temp_tau;
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
	for (uint8_t i=0 ; i<N_STATE ; i++)
		this->xn_int[i] = this->xn.getInterval();

}

real TM_val::fixpoint()
{
	uint8_t found_fixpoint = 0;
	uint8_t reduce_enclosure =1;
	real coeff = 1.0;
	real curr_tau;
	Interval enclosurej;
	Interval int_h;

	// Store the interval conversion of every derivatives
	for (uint8_t i =0 ; i< T_ORDER*N_STATE ; i++){
		this->der_int[i] = this->der[i].getInterval();
	}

	do {
		if (reduce_enclosure == 1){
			for(uint8_t index=0 ; index<N_STATE ; index++){
				this->apriori_enclosure[index] = this->xn[index].scaleAF(1.0f + coeff*MAX_WIDEN_FACTOR);
			}
			// Get the T_ORDER derivative of the dynamic function at the a priori enclosure
			ode_reminder(this->apriori_enclosure , this->rem);
			// Reduce initial tau when enclosure is too big
			curr_tau = this->tau_init * coeff;
			reduce_enclosure = 0;
		}
		found_fixpoint = 1;
		for (uint8_t i=0 ; i<N_STATE ; i++){
			int_h = Interval(0 , curr_tau);
			enclosurej = this->xn_int[i];
			for (uint8_t j=0 ; j<T_ORDER ; j++){
				enclosurej += this->der_int[i + j*N_STATE] * int_h;
				int_h.r_elem *= curr_tau;
			}
			enclosurej += this->rem[i] * int_h;
			if (AF1::subseteq(enclosurej , apriori_enclosure[i]) == 0){
				curr_tau *= STEP_DIVIDOR;
				if (curr_tau < H_MIN){
					reduce_enclosure = 1;
					coeff *= WIDEN_DIVIDOR;
				}
				found_fixpoint = 0;
				i = N_STATE; // Break
			}
			// found_fixpoint *= AF1::subseteq(enclosurej , apriori_enclosure[i]);
		}
	} while( found_fixpoint == 0)
	return curr_tau;
}