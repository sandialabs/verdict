/*=========================================================================

  Module:    VerdictVector.cpp

  Copyright (c) 2006 Sandia Corporation.
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/


/*
 *
 * VerdictVector.cpp contains implementation of Vector operations
 *
 * This file is part of VERDICT
 *
 */


#include "verdict.h"
#include <math.h>
#include "VerdictVector.hpp"
#include <float.h>

#if defined(__BORLANDC__)
#pragma warn -8004 /* "assigned a value that is never used" */
#endif

namespace VERDICT_NAMESPACE
{

const double TWO_VERDICT_PI = 2.0 * VERDICT_PI;

VerdictVector &VerdictVector::length(const double new_length)
{
  double len = this->length();
  xVal *= new_length / len;
  yVal *= new_length / len;
  zVal *= new_length / len;
  return *this;
}

double VerdictVector::interior_angle(const VerdictVector &otherVector)
{
  double cosAngle=0., angleRad=0., len1, len2=0.;
  
  if (((len1 = this->length()) > 0) && ((len2 = otherVector.length()) > 0))
    cosAngle = (*this % otherVector)/(len1 * len2);
  else
  {
    assert(len1 > 0);
    assert(len2 > 0);
  }
  
  if ((cosAngle > 1.0) && (cosAngle < 1.0001))
  {
    cosAngle = 1.0;
    angleRad = acos(cosAngle);
  }
  else if (cosAngle < -1.0 && cosAngle > -1.0001)
  {
    cosAngle = -1.0;
    angleRad = acos(cosAngle);
  }
  else if (cosAngle >= -1.0 && cosAngle <= 1.0)
    angleRad = acos(cosAngle);
  else
  {
    assert(cosAngle < 1.0001 && cosAngle > -1.0001);
  }
  
  return( (angleRad * 180.) / VERDICT_PI );
}

VerdictVector::VerdictVector(const double xyz[3]) 
  : xVal(xyz[0]), yVal(xyz[1]), zVal(xyz[2])
{}

} // namespace verdict
