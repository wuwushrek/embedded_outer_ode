#ifndef __AA_H__
#define __AA_H__

#include "config.h"
#include "interval.h"

class MAF1 {

private:
	real center;
	real deviations[N_NOISE];
	real err_term;
	static uint8_t index_last;

public:
	// Constructor
	MAF1(real center = 0.0f);
	MAF1(const MAF1 &);
	MAF1(const Interval&);

	// Destructor
	~MAF1();

	// Overloading operators
	MAF1 & operator = (const MAF1 &);
	MAF1 & operator = (const real);

	real & operator[](uint8_t );
	real operator[](uint8_t ) const;

	MAF1 operator + (const MAF1 &) const;
	MAF1 operator - (const MAF1 &) const;
	MAF1 operator * (const MAF1 &) const;
	MAF1 operator / (const MAF1 &) const;

	MAF1 & operator += (const MAF1 &);
	MAF1 & operator -= (const MAF1 &);
	// MAF1 operator ^ (const uint8_t) const;

	// unary operator minus
	MAF1 operator - () const;

	MAF1 operator * (const real) const;
	MAF1 operator / (const real) const;

	real getCenter() const;
	Interval getInterval() const;
	real getErrorTerm() const;

	real getMax() const;
	real getMin() const;
	real getRadius() const;

	void compress_af(real tol);

#ifdef VERBOSE
	void print_af(FILE *f = NULL);
#endif
	
	friend MAF1 cos(const MAF1 &);
	friend MAF1 sin(const MAF1 &);
	friend MAF1 tan(const MAF1 &);
	friend MAF1 inv(const MAF1 &);
};

MAF1 operator * (real , const MAF1 &);
MAF1 operator / (real , const MAF1 &);
MAF1 operator + (real , const MAF1 &);
MAF1 operator - (real , const MAF1 &);

#endif //__AA_H__