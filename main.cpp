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

	// Mesure time
	clock_t begin = clock();
#endif

	// Initialize states with the brusselator example
	AF1 x0[N_STATE];
	x0[0] = Interval(0.9 , 1);
	x0[1] = Interval(0 , 0.1);
	real t_end = 1.1;
	
	TM_val tm_val(x0 , t_end , 0.0);

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
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; //getTotalTime (start_time , end_time );
	std::cout << "[RESULT] elpased time (sec) =" << elapsed_secs << std::endl;
#endif
	
	return 0;
}