#include "ode_integr.h"
#include <iostream>

int main()
{
	// Initialize states
	AF1 x0[N_STATE];
	x0[0] = Interval(0.9 , 1);
	x0[1] = Interval(0 , 0.1);
	real t_end = 2;

	// AF1 der[4];
	// ode_derivatives(x0 , der);
	// for (uint8_t i=0 ; i< 4 ; i++){
	// 	print_af(der[i]);
	// 	print_it(der[i].getInterval());
	// }

	// AF1 val[2];
	// ode_function(x0 , val);
	// for(uint8_t i = 0 ; i < 2 ; i++){
	// 	print_af(val[i]);
	// 	print_it(val[i].getInterval());
	// }

	// AF1 rem[2];
	// ode_reminder(x0 , rem);
	// for(uint8_t i=0 ; i<2 ; i++){
	// 	print_af(rem[i]);
	// 	print_it(rem[i].getInterval());
	// }

	// print_af(x0[0]*x0[0]);
	// print_af(x0[0]*x0[1]);
	// print_af(x0[1]*x0[1]);
	// AF1 temp = x0[1]*x0[1];
	// print_af((x0[0]*x0[0] + x0[0]*x0[1] + x0[1]*x0[1])) ;
	// AF1 temp = ((x0[0]*x0[0] + x0[0]*x0[1] )*x0[1]);
	// print_af(temp);

	// print_af(temp * x0[1]) ;
	
	TM_val tm_val(x0 , t_end , 0.0);

	while( tm_val.tn <= t_end){
		tm_val.buildAndEval();
		std::cout << "Time : " << tm_val.tn << std::endl;
		for(uint8_t i=0 ; i< N_STATE ; i++){
			// std::cout << "state = " << (int) i << std::endl;
			// for (uint8_t j=0 ; j< T_ORDER ; j++){
			// 	print_it(tm_val.der_int[j+i*T_ORDER]);
			// }
			// print_af(tm_val.xn[i]);
#ifdef VERBOSE
			tm_val.xn[i].getInterval().print_it();
#endif
			// printf("Nb indexes : %d \n", tm_val.xn[i].nbIndexes);
			//print_af(tm_val.xn[i]);
		}
		printf("--------- \n");
	}
	return 0;
}