/*
 * aa_aaf.h -- Affine Arithmetic class
 * Copyright (C) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
 * Copyright (c) 2004 LIRIS (University Claude Bernard Lyon 1)
 * Copyright (C) 2009 LUH (Leibniz Universitaet Hannover)
 *
 * This file is part of aaflib.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with libaa; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef __AA_MOD2_H__
#define __AA_MOD2_H__

#include "interval.h"
#include "config.h"


//#include "vnode.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <list>

using namespace std;

// Affine Arithmetic Form
class DAF
{

 private:

  // central value
  real cvalue;
  // length of indexes
  unsigned length;
  // array size of indexes and deviations
  unsigned size;

#ifdef FAST_RAD
  real radius;
#endif

  // At creation we don't store null deviations

  // values of parial deviations
  real* deviations;
  // indexes of partial deviations
  unsigned * indexes;

  // highest deviation symbol in use
  static unsigned last;

 public:
  
  // constructors
  DAF(real v0 = 0.0);
  DAF(real, const real*, const unsigned *, unsigned);
  DAF(const DAF &);
  // ajout SP
  DAF(const Interval);

  // destructor
  ~DAF();

  real operator[](unsigned) const;

  DAF & operator = (const DAF &);
  DAF & operator = (const real);
  DAF operator + (const DAF &) const;
  DAF operator - (const DAF &) const;
  DAF operator * (const DAF &) const;
  DAF operator / (const DAF &) const;

  DAF operator - () const;
  DAF & operator += (real);
  DAF & operator -= (real);
  DAF & operator *= (real);
  DAF & operator /= (real);
  DAF & operator += (const DAF &);
  DAF & operator -= (const DAF &);
  DAF operator * (real);
  DAF operator / (real) const;

  friend std::ostream & operator << (std::ostream &, const DAF &);

  friend DAF cos(const DAF &);
  friend DAF sin(const DAF &);
  friend DAF tan(const DAF &);
  friend DAF inv(const DAF &);
 
  unsigned getlength() const;
  real getCenter() const;
  real getRadius() const;

#ifdef VERBOSE
  void print_af(FILE *f = NULL);
#endif

  Interval getInterval() const;
  void compress_af(real tol);

private:
  static unsigned inclast();
};

// binary operators
DAF operator * (real, const DAF);
DAF operator / (real, const DAF);
DAF operator + (real, const DAF);
DAF operator - (real, const DAF);


// DAF inline methods

/************************************************************
 * Method:        getCenter
 * Author & Date: ??? - ???
 * Description:   
 *   Returns the central value of the DAF 
 *
 *   Input  : -
 *   Output : real: center value
 ************************************************************/
inline real DAF::getCenter() const
{
  return cvalue;
}

/************************************************************
 * Method:        inclast
 * Author & Date: ??? - ???
 * Description:   
 *   increases the highest symbol
 *
 *   Input  : -
 *   Output : unsigned : highest symbol to use
 ************************************************************/
inline unsigned DAF::inclast()
{
  return ++last;
}

/************************************************************
 * Method:        getlength
 * Author & Date: ??? - ???
 * Description:   
 *   returns the number of partial deviations
 *
 *   Input  : -
 *   Output : unsigned : number of partial deviations
 ************************************************************/
inline unsigned DAF::getlength() const
{
  return length;
}

#endif  // __AA_MOD2_H__
