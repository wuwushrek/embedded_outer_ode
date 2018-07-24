#ifndef __AA_H__
#define __AA_H__

#include "config.h"
#include "interval.h"

class AF1 {

private:
	real center;
	real deviations[N_NOISE];
	real err_term;

public:
	// Constructor
	AF1(real center = 0.0f);
	AF1(const AF1 &);
	AF1(const Interval &);

	// Destructor
	~AF1();

	// Overloading operators
	AF1 & operator = (const AF1 &);
	AF1 & operator = (const real);

	AF1 operator + (const AF1 &) const;
	AF1 operator - (const AF1 &) const;
	AF1 operator * (const AF1 &) const;
	AF1 operator / (const AF1 &) const;
	AF1 operator ^ (const uint8_t) const;

	// unary operator minus
	AF1 operator - () const;

	AF1 operator * (const real) const;
	AF1 operator / (const real) const;

#ifdef VERBOSE
	void af1print(const char *) const;
#endif

	real getCenter() const;
	interval getRange() const;
	real getErrorTerm() const;

	friend AF1 cos(const AF1 &);
	friend AF1 sin(const AF1 &);
	friend AF1 tan(const AF1 &);
	friend AF1 power(const uint8_t);
	friend AF1 inv(const AF1 &);
};

AF1 operator * (real , const AF1 &);
AF1 operator / (real , const AF1 &);
AF1 operator + (real , const AF1 &);
AF1 operator - (real , const AF1 &);

#endif //__AA_H__