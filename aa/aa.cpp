#include "aa.h"
#include <cstring>

using namespace std;

uint8_t MAF1::index_last = 0;

/************************************************************/
/* Common definition of constructor + getter of affine form */
/************************************************************/
MAF1::MAF1(real center)
{
	this->center = center;;
	for (uint8_t index = 0; index < N_NOISE ; index++){
		this->deviations[index] = 0.0f;
	}
	this->err_term = 0.0f;
}

MAF1::MAF1(const MAF1 &af)
{
	this->center = af.center;
	memcpy(this->deviations , af.deviations , sizeof(af.deviations));
	this->err_term = af.err_term;
}

MAF1::MAF1(const Interval &it)
{
	this->center = it.getCenter();
	for (uint8_t index = 0; index < N_NOISE ; index++){
		this->deviations[index] = 0.0f;
	}
	if (index_last >= N_NOISE)
		this->err_term = it.getRadius();
	else {
		this->deviations[index_last] = it.getRadius();
		this->err_term = 0.0f;
		index_last++;
	}
	// assert_af(index_last <= N_NOISE);
}

MAF1::~MAF1(){}

MAF1 & MAF1::operator = (const real center)
{
	this->center = center;
	for (uint8_t index = 0; index < N_NOISE ; index++){
		this->deviations[index] = 0.0f;
	}
	this->err_term = 0.0f;
	return *this;
}

MAF1 & MAF1::operator = (const MAF1 &af)
{
	this->center = af.center;
	memcpy(this->deviations , af.deviations , sizeof(af.deviations));
	this->err_term = af.err_term;
	return *this;
}

real & MAF1::operator[](uint8_t index) 
{
	assert_af(index>=0 && index < N_NOISE);
	return this->deviations[index];
}

real MAF1::operator[](uint8_t index) const
{
	assert_af(index>=0 && index < N_NOISE);
	return this->deviations[index];
}

real MAF1::getCenter() const
{
	return this->center;
}

Interval MAF1::getInterval() const
{
	real radius = this->getRadius();
	return Interval(this->center - radius , this->center + radius);
}

real MAF1::getErrorTerm() const
{
	return this->err_term;
}

real MAF1::getRadius() const
{
	real radius = this->err_term;
	for (uint8_t index = 0 ; index < N_NOISE ; index ++){
		radius += abs(this->deviations[index]);
	}
	return radius;
}

real MAF1::getMax() const
{
	return this->getCenter() + this->getRadius();
}

real MAF1::getMin() const
{
	return this->getCenter() - this->getRadius();
}

void MAF1::compress_af(real tol)
{
	return ;
}

#ifdef VERBOSE
void MAF1::print_af(FILE *file)
{
	if (file != NULL){
		fprintf(file,"Center = %f \n", this->center);
		for (uint8_t i=0 ; i < N_NOISE ; i++){
			fprintf(file,"---> eps[%d] = %f \n", (int) i , this->deviations[i]);
		}
		fprintf(file,"-------------------------------\n");
		return ;
	}
	printf("Center = %f \n", this->center);
	for (uint8_t i=0 ; i < N_NOISE ; i++){
		printf("---> eps[%d] = %f \n", (int) i , this->deviations[i]);
	}
	printf("-------------------------------\n");
}
#endif

/************************************************************/
/* 	Affine form basic arithmetic operations definition   	*/
/************************************************************/

MAF1 MAF1::operator + (const MAF1 &other) const
{
	MAF1 temp(*this);
	temp.center += other.center;
	for (uint8_t index=0 ; index< N_NOISE ; index++){
		temp[index] += other.deviations[index];
	}
	temp.err_term += other.err_term;
	return temp;
}


MAF1 MAF1::operator - (const MAF1 &other) const
{
	MAF1 temp(*this);
	temp.center -= other.center;
	for (uint8_t index=0 ; index< N_NOISE ; index++){
		temp[index] -= other.deviations[index];
	}
	temp.err_term += other.err_term;
	return temp;
}

MAF1 MAF1::operator - () const
{
	MAF1 temp(*this);
	temp.center *= -1;
	for (uint8_t index=0 ; index< N_NOISE ; index++){
		temp[index] *= -1;
	}
	// err_term doesn't change
	return temp;
}

MAF1 MAF1::operator * (const real val) const
{
	MAF1 temp(*this);
	temp.center *= val;
	for (uint8_t index=0 ; index< N_NOISE ; index++){
		temp[index] *= val;
	}
	temp.err_term *= abs(val);
	return temp;
}

MAF1 MAF1::operator / (const real val) const
{
	assert_af(val != 0);

	MAF1 temp(*this);
	temp.center /= val;
	for (uint8_t index=0 ; index< N_NOISE ; index++){
		temp[index] /= val;
	}
	temp.err_term /= abs(val);
	return temp;
}

MAF1 & MAF1::operator += (const MAF1 &other)
{
	this->center += other.center;
	for (uint8_t i=0 ; i< N_NOISE ; i++){
		this->deviations[i] += other[i];
	}
	this->err_term += other.err_term;
	return *this;
}

MAF1 & MAF1::operator -= (const MAF1 &other)
{
	this->center -= other.center;
	for (uint8_t i=0 ; i< N_NOISE ; i++){
		this->deviations[i] -= other[i];
	}
	this->err_term += other.err_term;
	return *this;
}

/************************************************************/
/* 	Affine form basic arithmetic operations approximation  	*/
/************************************************************/

MAF1 MAF1::operator * (const MAF1 &other) const
{
	uint8_t index;
	MAF1 temp;

	// Center of the new affine form update
#ifndef FAST_MULT
	for (index = 0; index < N_NOISE ; index++){
		temp.center += other[index] * this->deviations[index];
	}
	temp.center /= 2;
#endif
	temp.center += other.center * this->center;

	// Noise terms update
	for(index = 0 ; index< N_NOISE ; index++){
		temp[index] = other.center* this->deviations[index] + this->center * other[index];
	}

	// accumulation Error term update
#ifndef FAST_MULT
	for(index = 0; index < N_NOISE ; index++){
		temp.err_term += abs(other[index] * this->deviations[index]);
	}
	temp.err_term /= -2;
#endif
	temp.err_term += other.err_term * abs(this->center) + this->err_term * abs(other.center) + other.getRadius()*this->getRadius();

	return temp;
}

MAF1 MAF1::operator / (const MAF1 &other) const
{
#ifndef FAST_DIV
	real r;
	real a , b;
	r = other.getRadius();
	assert_af(other.center != 0);

	if (r == 0)
		return *this / other.center;

	a = other.center - r;
	b = other.center + r;
	assert_af(a*b > 0);

	MAF1 temp;
	real curr_center = this->center / other.center;
	for (uint8_t index = 0 ; index <N_NOISE ; index++){
		temp[index] = this->deviations[index] - curr_center*other[index];
	}
	temp.err_term = abs(this->err_term - this->center*other.err_term);
	return temp * inv(other) + curr_center;
#else
	return (*this) * inv(other);
#endif
}

/*MAF1 MAF1::operator ^ (const uint8_t n) const
{
	if (n == 0)
		return MAF1(1.0);
	else if (n == 1)
		return *this;

	real a , b;
	real r;

	r = this->getRadius();

	if ( r == 0){
		a = 1;
		for(uint8_t i =0 ; i< n ; i++){
			a *= this->center;
		}
		return MAF1(a);
	}

	a = this->center - r;
	b = this->center + r;

	real fa , fb;
	real alpha , dzeta , delta;

	fa = 1;
	fb = 1;
	for (uint8_t i = 0 ; i< n ; i++){
		fa *= a;
		fb *= b;
	}

	if (r > MIN_RAD){
		alpha = (fb - fa)/(b-a);
	} else {
		alpha = n*fa/(a+EPS_ZERO);
	}

	real x_1;
	if ( n == 2)
		x_1 = - abs(alpha/2);
	else
		x_1 = -pow((real) (abs(alpha/n)), (real) (1.0/(n-1)));
	real x_2 = -x_1;
	real fx_1, fx_2;
	// check valid points
	if (x_1 > a){
		// x_1 is valid
		fx_1 = 1;
		for (uint8_t i = 0 ; i< n ; i++){
			fx_1 *= x_1;
		}
	}
	else {
		x_1 = a;
		fx_1 = fa;
	}

	if (x_2 < b) {
		fx_2 = 1;
		for (uint8_t i = 0 ; i< n ; i++){
			fx_2 *= x_2;
		}
	}
	else {
		x_2 = b;
		fx_2 = fb;
	}

	double y_2 = fx_2 - alpha*x_2;
	double y_1 = fx_1 - alpha*x_1;

	delta = 0.5*(y_1 - y_2);
	dzeta = 0.5*(y_1 + y_2);

	// inverse of other
	MAF1 temp;
	temp.center = alpha * this->center + dzeta;
	for (uint8_t index = 0 ; index < N_NOISE ; index++){
		temp[index] = this->deviations[index] * alpha;
	}
	temp.err_term = this->err_term + delta;

	return temp;
}*/

MAF1 inv(const MAF1 &other)
{
	real a , b;
	real r;

	r = other.getRadius();

	if (r == 0)
		return MAF1(1.0/other.center);

	a = other.center - r;
	b = other.center + r;

	real fa , fb;
	real alpha ,dzeta , delta;

	assert_af(a*b > 0.0f);

	fa = 1/a;
	fb = 1/b;

	alpha = - fa * fb;
	if (a > 0){
		delta = (fa + fb)/2  - (1/sqrt(a*b));
		dzeta = fa + fb - delta;
	}else {
		delta = -((fa + fb)/2  + (1/sqrt(a*b)));
		dzeta = fa + fb + delta;
	}

	// inverse of other
	MAF1 temp;
	temp.center = alpha * other.center + dzeta;
	for (uint8_t index = 0 ; index < N_NOISE ; index++){
		temp[index] = other[index] * alpha;
	}
	temp.err_term = other.err_term + delta;

	return temp;
}

MAF1 operator * (real val, const MAF1 &af)
{
	return af * val;
}

MAF1 operator / (real val, const MAF1 &af)
{
	return af / val;
}

MAF1 operator + (real val , const MAF1 &af)
{
	return af + val;
}

MAF1 operator - (real val, const MAF1 &af)
{
	return af - val;
}


/************************************************************/
/* 	Affine form approximation of trigonometry functions  	*/
/************************************************************/

/* TODO: NEED to optimize PI/2 PI/4 2PI by pre calculating them */
MAF1 sin(const MAF1 &other)
{
	real a , b;
	real r;

	r = other.getRadius();

	if (r == 0)
		return MAF1(sin(other.center));

	assert_af( r < M_PI/4); // For the moment

	/*if ( r >= 2*M_PI ){
		MAF1 res;
		res.err_term = 1.0f;
		return res;
	}*/

	a = other.center - r;
	b = other.center + r;

	if (a> M_PI || a < -M_PI){
		real temp = floor(a / (2*M_PI));
		a -= temp * 2 * M_PI;
		a -= (a > M_PI) ? 2*M_PI : 0;
		b = a + 2*r;
	}

	real fa , fb;
	real alpha ,dzeta , delta;

	fa = sin(a);
	fb = sin(b);

	if( b < 0 || (b < M_PI  && a > 0)){ // chebyshev approx in  this case
		if (r > MIN_RAD){
			alpha = (fb - fa) / (b - a);
			real sol = ((a > 0) - (a < 0)) * acos(alpha);
			real fsol = sqrt(1 - alpha*alpha); // fast computation of sin(acos)
			dzeta = (fa + fsol - alpha * (a + sol)) / 2;
			delta = abs(fsol - fa - alpha * (sol - a)) / 2;
		}
		else{
			alpha = cos(a);
			dzeta = fa - alpha*a;
			delta = 0.0f;
		}

	} else { // min range optimization since the derivative cant be zero because of the radius constraints
		if ( a <= 0 ){
			alpha = 1.0f;
			delta = (-fb + fa - alpha * (a - b)) / 2;
		} else {
			alpha = -1.0f;
			delta = (fb - fa + alpha * (a - b)) / 2;
		}
		dzeta = (fa + fb - alpha * (a + b)) / 2;
	}

	MAF1 temp;
	temp.center = alpha * other.center + dzeta;
	for (uint8_t index = 0 ; index < N_NOISE ; index++){
		temp[index] = other[index] * alpha;
	}
	temp.err_term = other.err_term + delta;

	return temp;
}

MAF1 cos(const MAF1 &other)
{
	return sin(other + M_PI/2);
}

MAF1 tan(const MAF1 &other)
{
	real a , b;
	real r;

	r = other.getRadius();

	if (r == 0)
		return MAF1(tan(other.center));

	a = other.center - r;
	b = other.center + r;

	// In case the angles are not between -PI/2 and PI/2
	if (a>= M_PI/2 || a <= -M_PI/2){
		real temp = floor(a / (M_PI));
		a -= temp *  M_PI;
		a -= (a > (M_PI/2)) ? M_PI : 0;
		b = a + 2*r;
	}

	assert_af(! ( (a <= M_PI/2 && b >= M_PI/2 ) || (a <= -M_PI/2 && b >= -M_PI/2 ) ));

	real fa , fb;
	real alpha ,dzeta , delta;

	fa = tan(a);
	fb = tan(b);

	if (r > MIN_RAD)
		alpha = (fb - fa) / (b-a);
	else
		alpha =  1 + fa * fb;

	real aux = sqrt(alpha - 1.0);
	real eps = atan(aux);
	// real eps2 = -eps1;

	if ( eps >= a && eps <= b && -eps >= a && -eps <= b){
		dzeta = 0;
		delta = abs(aux - alpha * eps);
	} else {
		if (-eps >= a && -eps <= b){
			eps = -eps;
			aux = -aux;
		}
		dzeta = (fa + aux - alpha*(a + eps)) / 2;
		delta = abs(aux - fa - alpha * (eps - a)) / 2;
	}

	MAF1 temp;
	temp.center = alpha * other.center + dzeta;
	for (uint8_t index = 0 ; index < N_NOISE ; index++){
		temp[index] = other[index] * alpha;
	}
	temp.err_term = other.err_term + delta;

	return temp;
}
