#include "ode_integr.h"
#include <iostream>
#include <ctime>

#ifdef VERBOSE
FILE* mfile[N_STATE];
#endif

int main()
{

#ifdef VERBOSE
	for (int i=0 ; i<N_STATE ; i++){
		std::ostringstream stream;
#ifdef USE_MAF1
		stream << "x" << (int) i << "outer_MAF1.out";
#elif USE_MAF2
		stream << "x" << (int) i << "outer_MAF2.out";
#else
		stream << "x" << (int) i << "outer_DAF.out";
#endif
		mfile[i] = fopen(stream.str().c_str(),"w+");
	}
#endif
	// Mesure time
	clock_t begin = clock();

	// Initialize states with the brusselator example
	AF1 x0[N_STATE];
	/*x0[0] = Interval(0.9 , 1);
	x0[1] = Interval(0 , 0.1);
	real t_end = 1.1;*/
	x0[0] = 0.0; //roll
	x0[1] = 0.0; //pitch
	x0[2] = 0.0; // yaw
	x0[3] = Interval(-0.06 , 0.06); // p
	x0[4] = Interval(-0.06 , 0.06); // q
	x0[5] = 0.0; // r
	x0[6] = 0.0; // err p
	x0[7] = 0.0; // err q
	x0[8] = 0.0; // err r
	x0[9] = 0.0; // u
	x0[10] = 0.0;// v
	x0[11] = 0.0;// w
	x0[12] = Interval(0 , 0.1); // z
	x0[13] = 0.0; // integral of Z
	real t_end = 0.3;

	OdeFunc f;
	/*f.params[0] = 1.5;
	f.params[1] = 1;*/
	f.params[0] = 0.0; // p value
	f.params[1] = 0.0; // q value
	f.params[2] = 0.0; // r value
	f.params[3] = 1.0; // z desired setpoint

	// AF1 yp[N_STATE];
	// // AF1 test = x0[0]*x0[0]*x0[1];

	// AF1 test = 1 - (f.params[1]+1)*x0[0] + f.params[0]*x0[0]*x0[0]*x0[1];
	// test.getInterval().print_it();
	// f(yp , x0);
	// for(uint8_t i=0 ; i<N_STATE ; i++){
	// 	yp[i].getInterval().print_it();
	// }
	TM_val tm_val(x0, 0, &f); 
	
//	TM_val tm_val(x0 , t_end , 0.0);

	while( tm_val.tn <= t_end){
#ifdef VERBOSE
		for (uint8_t i = 0 ; i< N_STATE ; i++){
			fprintf(mfile[i] , "%f\t", tm_val.tn);
			tm_val.xn[i].getInterval().print_it(mfile[i]);
		}
#endif
		tm_val.buildAndEval();
	}

#ifdef VERBOSE
	for (uint8_t i = 0 ; i< N_STATE ; i++){
		fprintf(mfile[i] , "%f\t", tm_val.tn);
		tm_val.xn[i].getInterval().print_it(mfile[i]);
	}
#endif
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; //getTotalTime (start_time , end_time );
	std::cout << "[RESULT] elpased time (sec) =" << elapsed_secs << std::endl;
	
	return 0;
}