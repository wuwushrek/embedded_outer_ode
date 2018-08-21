/*
 * aa_aafarithm.c -- Affine arithmetical operations
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
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


#include "aa_mod2.h"
#include <cmath>
#include <algorithm>

#include <iostream>
#include <ostream>
#include <fstream>
#include <cstdio>

using namespace std;

// static variables
// defualt starting index for deviations
unsigned DAF::last = 0;

/************************************************************
 * Method:        DAF
 * Author & Date: ??? - ???
 * Description:   
 *   Default constructor creating an DAF without deviations
 *
 *   Input  : real : corresponding center value
 *   Output : -
 ************************************************************/
DAF::DAF(real v0):
  cvalue(v0), 
  length(0),
  size(0),
#ifdef FAST_RAD
  radius(0.0),
#endif 
  deviations(NULL),
  indexes(NULL)
{
}


/************************************************************
 * Method:        DAF
 * Author & Date: ??? - ?DAF operator + (real, const DAF);??
 * Description:   
 *   Default constructor creating an DAF with deviations
 *
 *   Input  : real     : corresponding center value
 *            real *   : array of index values
 *            unsigned * : array of indices
 *            unsigned   : # of indices
 *   Output : -
 ************************************************************/
DAF::DAF(real v0, const real * t1, const unsigned * t2, unsigned T):
  cvalue(v0), 
  length(T),
  size(T),
#ifdef FAST_RAD
  radius(0.0),
#endif
  deviations(new real[T]),
  indexes(new unsigned[T])
{
  for (unsigned i = 0; i < length; i++)
  {
    deviations[i] = t1[i];
#ifdef FAST_RAD
    radius += fabs(deviations[i]);
#endif
    indexes[i] = t2[i];
  }

  if (indexes[length-1] > last) 
    last = indexes[length-1];
}


/************************************************************
 * Method:        DAF
 * Author & Date: ??? - ???
 * Description:   
 *   Copy constructor
 *
 *   Input  : DAF : DAF to be copied
 *   Output : -
 ************************************************************/
DAF::DAF(const DAF &P):
  cvalue(P.cvalue), 
  length(P.length),
  size(P.size),
#ifdef FAST_RAD
  radius(P.radius),
#endif 
  deviations(NULL),
  indexes(NULL)
{
  if (size)
  {
    deviations = new real [size];
    indexes = new unsigned [size];
  }

  for (unsigned i = 0; i < length; i++)
  {
    deviations[i] = P.deviations[i];
    indexes[i] = P.indexes[i];
  }

}


/************************************************************
 * Method:        DAF
 * Author & Date: ??? - ???
 * Description:   
 *   Constructor creating DAF from interval data
 *
 *   Input  : AAInterval           : Interval
 *   Output : -
 ************************************************************/
DAF::DAF(const Interval iv):
  cvalue(iv.getCenter()), 
  length(1),
  size(1),
#ifdef FAST_RAD
  radius(0.0),
#endif 
  deviations(new real[1]),
  indexes(new unsigned[1])
{  
  deviations[0] = iv.getRadius();

  indexes[0] = inclast();

#ifdef FAST_RAD
  radius = fabs(deviations[0]);
#endif

}


/************************************************************
 * Method:        ~DAF
 * Author & Date: ??? - ???
 * Description:   
 *   Destructor
 *
 *   Input  : -
 *   Output : -
 ************************************************************/
DAF::~DAF()
{
  if (size)
  {
    delete [] deviations;
    delete [] indexes;
  }
}


/************************************************************
 * Method:        =
 * Author & Date: ??? - ???
 * Description:   
 *   Affectation operator
 *
 *   Input  : DAF
 *   Output : -
 ************************************************************/
DAF & DAF::operator = (const real d)
{
  length = 0;
  cvalue = d;

#ifdef FAST_RAD
  radius = 0.0;
#endif

  return *this;
}

/************************************************************
 * Method:        =
 * Author & Date: ??? - ???
 * Description:   
 *   Affectation operator
 *
 *   Input  : DAF
 *   Output : -
 ************************************************************/
DAF & DAF::operator = (const DAF & P)
{
  unsigned plength = P.getlength();

  if (&P != this)
  {
    if (size < plength)
    {
      if (size)
      {
  delete [] deviations;
  delete [] indexes;
      }
      size = plength;
      if (size)
      {
  deviations = new real [size];
  indexes = new unsigned [size];
      }
    }

    cvalue = P.cvalue;
    length = plength;

    for (unsigned i = 0; i < length; i++)
    {
      deviations[i] = P.deviations[i];
      indexes[i] = P.indexes[i];
    }

#ifdef FAST_RAD
    radius = P.radius;
#endif
  }  

  return *this;
}

std::ostream & operator << (std::ostream & s, const DAF & P)
{

  // s.setf(0, ios_base::floatfield);
  s << "-------------\n";
  s << "Length = " << P.length << "\n";
  s << "v0     = " << P.cvalue << "\n";
#ifdef FAST_RAD
  s << "Radius = " << P.getRadius() << "\n";
#endif

  for (unsigned i = 0; i < P.length ; i++)
  {
    s << "e" << P.indexes[i] << " -> " << P.deviations[i] << "\n";
  }
  s << "-------------\n";

  return s;
}

#ifdef VERBOSE
void DAF::print_af(FILE *file)
{
  if (file != NULL){
    std::ostringstream stream;
    stream << *this << std::endl;
    std::string str =  stream.str();
    fprintf(file, "%s" , str.c_str());
    return;
  }
  std::cout << *this << std::endl;
}
#endif

void DAF::compress_af(real tol)
{
  real r = tol * this->getRadius();
  real rest = 0.0f;
  uint16_t index = 0;
  for(uint16_t i=0 ; i < length ; i++){
    real temp = abs(this->deviations[i]);
    if(temp >= r){
      this->deviations[index] = this->deviations[i];
      this->indexes[index] = this->indexes[i];
      index++;
    }else {
      rest += temp;
    }
  }
  if (index < length){
    deviations[index] = rest;
    indexes[index] = inclast();
    length = index+1;
  }
}

real DAF::getRadius() const
{
#ifdef FAST_RAD
  return radius;
#else
  real sum = 0.0;

  for (unsigned i = 0; i < length; i++)
  {
    if (deviations[i] >= 0.0)
      sum+=deviations[i];
    else
      sum+=-deviations[i];
  }

  return sum;
#endif
}

Interval DAF::getInterval() const
{
  real r = getRadius();
  return Interval(getCenter()-r , getCenter()+r);
}
/********************************************************************/


/************************************************************
 * Operator:      +=
 * Author & Date: ??? - ???
 * Description:   
 *   Addition is an affine operation
 *
 *   Input  : real : real value to be added
 *   Output : DAF
 ************************************************************/
DAF & DAF::operator += (real cst)
{
  cvalue += cst;
  return (*this);
}


/************************************************************
 * Operator:      -=
 * Author & Date: ??? - ???
 * Description:   
 *   Subtraction is an affine operation
 *
 *   Input  : real : real value to be subtracted
 *   Output : DAF
 ************************************************************/
DAF & DAF::operator -= (real cst)
{
  cvalue -= cst;
  return (*this);
}


/************************************************************
 * Operator:      *=
 * Author & Date: ??? - ???
 * Description:   
 *   Scalar multiplication is an affine operation
 *
 *   Input  : real : real value to be multiplied by
 *   Output : DAF
 ************************************************************/
DAF & DAF::operator *= (real cst)
{
  cvalue *= cst;
  for (unsigned int i = 0; i < length; i++)
    deviations[i] *= cst;

#ifdef FAST_RAD
  radius *= fabs(cst);
#endif

  return (*this);
}

/************************************************************
 * Operator:      /=
 * Author & Date: ??? - ???
 * Description:   
 *   Scalar division is an affine operation
 *
 *   Input  : real : real value to be divided by
 *   Output : DAF
 ************************************************************/
DAF & DAF::operator /= (real cst)
{
  cst = 1.0/cst;
  cvalue *= cst;
  for (unsigned int i = 0; i < length; i++)
    deviations[i] *= cst;

#ifdef FAST_RAD
  radius *= fabs(cst);
#endif

  return (*this);
}


/************************************************************
 * Operator:      +
 * Author & Date: ??? - ???
 * Description:   
 *   Affine addition
 *
 *   Input  : const DAF : DAF to be added
 *   Output : DAF
 ************************************************************/
DAF DAF::operator + (const DAF & P) const
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    DAF Temp(cvalue+P.cvalue);
    return (Temp);
  }

  if (l1 == 0)
  {
    DAF Temp(P);
    Temp += cvalue;
    return (Temp);
  }

  if (l2 == 0)
  {
    DAF Temp(*this);
    Temp += P.cvalue;
    return (Temp);
  }

  // Create our resulting DAF
  DAF Temp(cvalue + P.cvalue);

  unsigned * id1=indexes;
  unsigned * id2=P.indexes;

  real * va1=deviations;
  real * va2=P.deviations;

  unsigned * pu1=id1;
  unsigned * pu2=id2;

  if (l1+l2)
    Temp.indexes = new unsigned [l1+l2]; // the indexes of the result
  unsigned * idtemp = Temp.indexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin-idtemp;

  if (ltemp)
    Temp.deviations = new real [ltemp];
  real * vatempg = Temp.deviations;

  Temp.length = ltemp;
  Temp.size = ltemp;

  // Fill the deviations array
  // of the resulting DAF

  for (unsigned i = 0; i < ltemp; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i]=va2[b];  // va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i]=va1[a];  // va1[a]+0
      pu1++;
      continue;
    }
    
    vatempg[i]=va1[a] + va2[b];
    pu1++;
    pu2++;
  }

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    Temp.radius += fabs(vatempg[i]);
#endif
  
  return Temp;
}


/************************************************************
 * Operator:      +=
 * Author & Date: Darius Grabowski - 07/2007
 * Description:   
 *   Affine addition
 *
 *   Input  : const DAF : DAF to be added
 *   Output : DAF &     : *this
 ************************************************************/
DAF & DAF::operator += (const DAF & P)
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    cvalue += P.cvalue;
    return (*this);
  }

  if (l1 == 0)
  {
    real c = cvalue;
    *this = P;
    cvalue += c;
    return (*this);
  }

  if (l2 == 0)
  {
    cvalue += P.cvalue;
    return (*this);
  }

  // Create our resulting DAF

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  real * va1 = deviations;
  real * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  unsigned * tempIndexes = NULL;
  real * tempDeviations = NULL;

  if (l1+l2)
    tempIndexes = new unsigned [l1+l2]; // the indexes of the result
  unsigned * idtemp = tempIndexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin-idtemp;

  if (ltemp)
    tempDeviations = new real [ltemp];
  real * vatempg = tempDeviations;

  // Fill the deviations array
  // of the resulting DAF

  for (unsigned i = 0; i < ltemp; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i]=va2[b];  // va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i]=va1[a];  // va1[a]+0
      pu1++;
      continue;
    }
    
    vatempg[i]=va1[a] + va2[b];
    pu1++;
    pu2++;
  }

  // set new properties
  length = ltemp;
  size = ltemp;
  delete [] deviations;
  delete [] indexes;
  deviations = tempDeviations;
  indexes = tempIndexes;
  cvalue += P.cvalue;

#ifdef FAST_RAD
  radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    radius += fabs(vatempg[i]);
#endif
  
  return *this;
}


/************************************************************
 * Operator:      -
 * Author & Date: ??? - ???
 * Description:   
 *   Affine subtraction
 *
 *   Input  : const DAF : DAF to be subtracted
 *   Output : DAF
 ************************************************************/
DAF DAF::operator - (const DAF & P) const
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    DAF Temp(cvalue-P.cvalue);
    return (Temp);
  }

  if (l1 == 0)
  {
    DAF Temp(P);
    Temp *= -1.0;
    Temp += cvalue;
    return (Temp);
  }

  if (l2 == 0)
  {
    DAF Temp(*this);
    Temp -= P.cvalue;
    return (Temp);
  }

  // Create our resulting DAF
  DAF Temp(cvalue-P.cvalue);

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  real * va1 = deviations;
  real * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  Temp.indexes = new unsigned [l1+l2];
  unsigned * idtemp = Temp.indexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin-idtemp;

  Temp.deviations = new real [ltemp];
  real * vatempg = Temp.deviations;

  Temp.length = ltemp;
  Temp.size = ltemp;

  // Fill the deviations array
  // of the resulting DAF
  
  for (unsigned i=0;i<ltemp;i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i] = -va2[b];  // 0-va2[b]
      pu2++;
      continue;
    }
    
    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i] = va1[a];  // va1[a]-0
      pu1++;
      continue;
    }

    vatempg[i]=va1[a]-va2[b];
    pu1++;
    pu2++;
  }

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    Temp.radius += fabs(vatempg[i]);
#endif
  
  return Temp;
}


/************************************************************
 * Operator:      -=
 * Author & Date: Darius Grabowski - 07/2007
 * Description:   
 *   Affine addition
 *
 *   Input  : const DAF : DAF to be subtracted
 *   Output : DAF &     : *this
 ************************************************************/
DAF & DAF::operator -= (const DAF & P)
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    cvalue -= P.cvalue;
    return (*this);
  }

  if (l1 == 0)
  {
    real c = cvalue;
    *this = -P;
    cvalue += c;
    return (*this);
  }

  if (l2 == 0)
  {
    cvalue -= P.cvalue;
    return (*this);
  }

  // Create our resulting DAF

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  real * va1 = deviations;
  real * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  unsigned * tempIndexes = NULL;
  real * tempDeviations = NULL;

  if (l1+l2)
    tempIndexes = new unsigned [l1+l2]; // the indexes of the result
  unsigned * idtemp = tempIndexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin-idtemp;

  if (ltemp)
    tempDeviations = new real [ltemp];
  real * vatempg = tempDeviations;

  // Fill the deviations array
  // of the resulting DAF

  for (unsigned i = 0; i < ltemp; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i]=-va2[b];  // -va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i]=va1[a];  // va1[a]-0
      pu1++;
      continue;
    }
    
    vatempg[i]=va1[a] - va2[b];
    pu1++;
    pu2++;
  }

  // set new properties
  length = ltemp;
  size = ltemp;
  delete [] deviations;
  delete [] indexes;
  deviations = tempDeviations;
  indexes = tempIndexes;
  cvalue -= P.cvalue;

#ifdef FAST_RAD
  radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    radius += fabs(vatempg[i]);
#endif
  
  return *this;
}


/************************************************************
 * Operator:      -
 * Author & Date: ??? - ???
 * Description:   
 *   Unary operator
 *
 *   Input  : const DAF : DAF to be multiplied by -1
 *   Output : DAF
 ************************************************************/
DAF DAF::operator - () const
{
  DAF Temp(*this);

  Temp.cvalue = -(Temp.cvalue);
  for (unsigned i = 0; i < length; i++)
  {
    Temp.deviations[i] = -(Temp.deviations[i]);
  }

  return Temp;
}


/************************************************************
 * Operator:      *
 * Author & Date: ??? - ???
 * Description:   
 *   Mul by a constant (on right), Affine operation
 *
 *   Input  : real : real to be multiplyed by
 *   Output : DAF
 ************************************************************/
DAF DAF::operator * (real cst)
{
  DAF Temp(*this);
  Temp.cvalue = cst*cvalue;

  for (unsigned i = 0; i < length; i++)
  {
    Temp.deviations[i] = cst*(Temp.deviations[i]);
  }

#ifdef FAST_RAD
  Temp.radius = fabs(cst)*radius;
#endif

  return Temp;
}


/************************************************************
 * Operator:      /
 * Author & Date: ??? - ???
 * Description:   
 *   Division by a constant (on right), Affine operation
 *
 *   Input  : real : real to be divided by
 *   Output : DAF
 ************************************************************/
DAF DAF::operator / (real cst) const
{
  DAF Temp(*this);
  Temp.cvalue = cvalue/cst;

  for (unsigned i = 0; i < length; i++)
  {
    Temp.deviations[i] = (Temp.deviations[i])/cst;
  }

#ifdef FAST_RAD
  Temp.radius = radius/fabs(cst);
#endif
  return Temp;
}

// -- Non member DAF functions --


/************************************************************
 * Operator:      *
 * Author & Date: ??? - ???
 * Description:   
 *   Mul by a constant (the left case)
 *
 *   Input  : real : real factor
 *            DAF    : DAF factor
 *   Output : DAF
 ************************************************************/
DAF operator * (real cst, const DAF P)
{
  DAF Temp(P);
  return Temp*cst;
}


/************************************************************
 * Operator:      +
 * Author & Date: ??? - ???
 * Description:   
 *   Add a constant (the left case)
 *
 *   Input  : real : real summand
 *            DAF    : DAF summand
 *   Output : DAF
 ************************************************************/
DAF operator + (real cst, const DAF P)
{
  DAF Temp(P);
  Temp += cst;
  return (Temp);
}


/************************************************************
 * Operator:      +
 * Author & Date: ??? - ???
 * Description:   
 *   Sub a constant (the left case)
 *
 *   Input  : real : real 
 *            DAF    : DAF 
 *   Output : DAF
 ************************************************************/
DAF operator - (real cst, DAF P)
{
  DAF Temp = -P;
  Temp += cst;
  return (Temp);
}

DAF DAF::operator * (const DAF & P) const
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    DAF Temp(cvalue*P.cvalue);
    return Temp;
  }
  if (l1 == 0)
  {
    // if *this is real
    DAF Temp(P);
    Temp *= cvalue;
    return Temp;
  }
  if (l2 == 0)
  {
    // if P is real
    DAF Temp(*this);
    Temp *= P.cvalue;
    return Temp;
  }

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  real * va1 = deviations;
  real * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  DAF Temp(cvalue*P.cvalue);  // Create our resulting DAF

  Temp.indexes = new unsigned [l1+l2+1];
  
  unsigned * idtemp=Temp.indexes;

  // Fill the indexes array

  unsigned * fin = std::set_union(id1, id1 + l1, id2, id2 + l2, idtemp);
  unsigned ltemp = fin - idtemp;

  Temp.deviations = new real [ltemp + 1];
  real * vatempg = Temp.deviations;

  Temp.length = ltemp + 1;
  Temp.size = Temp.length;
  
  real commonTermCenter = 0.0;
  real commonTermDeviation = 0.0;

  // Fill the deviations array

  for (unsigned i = 0; i < ltemp; i++)
  {
    unsigned a = pu1 - id1;
    unsigned b = pu2 - id2;
    
    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i] = cvalue*va2[b];  // cvalue*va2[b]+(P.cvalue)*0
      pu2++;
      continue;
    }
    
    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i] = (P.cvalue)*va1[a];  // cvalue*0+(P.cvalue)*va1[a]
      pu1++;
      continue;
    }
    
    vatempg[i] = cvalue*va2[b] + (P.cvalue)*va1[a];
    commonTermCenter += va2[b]*va1[a];
    commonTermDeviation += fabs(va2[b]*va1[a]);
    pu1++;
    pu2++;
  }

  // Compute the error
  // in a new deviation symbol
  real delta = getRadius()*P.getRadius();

  Temp.indexes[ltemp] = inclast();

    Temp.deviations[ltemp] = delta;
  
  // consider deviations occuring in both expressions  
  commonTermCenter *= 0.5;
  commonTermDeviation *= 0.5;
  Temp.cvalue += commonTermCenter;
  Temp.deviations[ltemp] -= commonTermDeviation;

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < Temp.length; i++)
    Temp.radius += fabs(Temp.deviations[i]);
#endif
  
  return Temp;
}

DAF DAF::operator / (const DAF & P) const
{
  if (this == &P)
    return DAF(1.0);
  else
    return (*this)*inv(P);
}

DAF inv(const DAF & P)
{
  real a, b;
  real fa, fb;
  real r;
  real alpha, dzeta, delta;

  if (P.length == 0)
  {
    DAF Temp(1.0/(P.cvalue));
    return Temp;
  }

  r = P.getRadius();

  a = P.cvalue - r;
  b = P.cvalue + r;

  assert_af(a*b < 0);
  

  fa = 1/a;
  fb = 1/b;

    alpha = -fa*fb;    
    real u = sqrt(a*b);

    if (a > 0)
    {
      delta = +0.5*(fa+fb-2.0/u);
      dzeta = fa+fb-delta;
    }
    else
    {
      delta = -0.5*(fa+fb+2.0/u);
      dzeta = fa+fb+delta;
    }
  // z0 = alpha*x0 + dzeta

  DAF Temp(alpha*(P.cvalue) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new real [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi

  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }  
  
  // Compute the error in a new deviation symbol
  // zk = delta
  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif

  return Temp;
}