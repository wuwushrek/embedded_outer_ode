#include "aa.h"
#include <iostream>

void print_af(const AF1 &af)
{
	std::cout << "---------" <<std::endl;
	std::cout << "center = " << af.getCenter() << std::endl;
	for (uint8_t i =0 ; i<N_NOISE ; i++){
		std::cout << "---> eps[" << (int) i << "] = " << af[i] <<std::endl;
	}
	std::cout << "err_term = " << af.getErrorTerm() << std::endl;
	std::cout << "---------" <<std::endl;
}

void print_it(const Interval &it)
{
	std::cout << "it = [" << it.getMin() << " , " << it.getMax() <<"]" << std::endl;
}

int main()
{
	AF1 af1(Interval(-1 , 1));
	AF1 af2(Interval(-1 , 1));
	AF1 af4(Interval(2 , 3));
	AF1 af3(4);

	print_af(af1);
	print_it(af2.getInterval());

	// Addition
	AF1 res = af1 + af2;
	print_af(res);

	// Substraction
	res = af1 - af2;
	print_af(res);

	// multiplication
	res = af1 * af2;
	print_af(res);

	//Division
	res = af1 / af4;
	print_af(res);
	print_it(res.getInterval());

	res = af4 / af4;
	print_af(res);

	res = af1 / af3;
	print_af(res);

	// pow
	res = af1*af1;
	print_af(res);
	print_it(res.getInterval());

	//sin
	res = sin((af1+2) * (M_PI / 20));
	print_af(res);
	print_it(res.getInterval());

	//cos
	res = cos((af1+2) * (M_PI / 20));
	print_af(res);
	print_it(res.getInterval());

	//tan
	res = tan((af1+2) * (M_PI / 20));
	print_af(res);
	print_it(res.getInterval());

	res = sin((af1+2) * (M_PI / 20)) / cos((af1+2) * (M_PI / 20));
	print_af(res);
	print_it(res.getInterval());

	return 0;
}