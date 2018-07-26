#include "interval.h"

Interval::Interval(real val)
{
	this->l_elem = val;
	this->r_elem = val;
}

Interval::Interval(real lo , real hi)
{
	assert_af(lo <= hi);
	this->l_elem = lo;
	this->r_elem = hi;
}

Interval::Interval(const Interval &it)
{
	this->l_elem = it.l_elem;
	this->r_elem = it.r_elem;
}

real Interval::getCenter() const
{
	return (this->l_elem + this->r_elem) / 2;
}

real Interval::getMin() const
{
	return this->l_elem;
}

real Interval::getMax() const
{
	return this->r_elem;
}

real Interval::getRadius() const
{
	return (this->r_elem - this->l_elem)/2;
}

real Interval::getRange() const
{
	return (this->r_elem - this->l_elem);
}

Interval & Interval::operator = (const Interval &it)
{

	this->l_elem = it.l_elem;
	this->r_elem = it.r_elem;
	return *this;
}