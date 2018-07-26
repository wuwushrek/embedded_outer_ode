#ifndef __INTERVAL_H__
#define __INTERVAL_H__

#include "config.h"

class Interval
{
private:
	real l_elem;
	real r_elem;

public:
	Interval(real val = 0.0f);
	Interval(const Interval&);
	Interval(real val1 , real val2);

	Interval & operator = (const Interval &);

	real getCenter() const;
	real getMin() const;
	real getMax() const;
	real getRadius() const;
	real getRange() const;
};

#endif //__INTERVAL_H__