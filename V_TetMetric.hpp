/*=========================================================================
  Copyright (c) 2006 Sandia Corporation.
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/


/*
 *
 * V_TetMetric.hpp contains declarations of tet related shape quantities
 *
 * This file is part of VERDICT
 *
 */

#ifndef VERDICT_TET_METRIC_HPP
#define VERDICT_TET_METRIC_HPP

#include <assert.h>
#include <math.h>
#include "VerdictVector.hpp"


namespace verdict {

//
//  want to avoid direct dependence of verdict on STK for now, this is the one function it really needs so just
//  copy it here for now.
//
    inline double stk_sqrt(double x) {
      return std::sqrt(x);
    }
    template <typename T>
    inline T stk_sqrt(T x) {
      return stk::math::sqrt(x);
    }










  //
  //  Compute the characteristic length scale of an element, this specifically influnces the explicit time step on the element
  //    edge0 = c1-c0;
  //    edge1 = c2-c0;
  //    normal = Cross(edge0, edge1)  # Note, not a unit normal
  //    p = centroid-c0
  //    normal^ = UnitVector(normal)
  //    dist = Dot(p, normal)
  //    
  //



template <typename T>
T tet_dimension( int num_nodes, const T* coordinates ) {

  const T* c0 = coordinates;
  const T* c1 = coordinates+3;
  const T* c2 = coordinates+6;
  const T* c3 = coordinates+9;

  T centroid[3];
  centroid[0] = 0.25*(c0[0] + c1[0] + c2[0] + c3[0]);
  centroid[1] = 0.25*(c0[1] + c1[1] + c2[1] + c3[1]);
  centroid[2] = 0.25*(c0[2] + c1[2] + c2[2] + c3[2]);
  //
  //  Most simplified calculation, distance from centroid to each face should be identical, so just calculate one
  //
  T edge0[3];
  edge0[0] = c1[0]-c0[0];
  edge0[1] = c1[1]-c0[1];
  edge0[2] = c1[2]-c0[2];

  T edge1[3];
  edge1[0] = c2[0]-c0[0];
  edge1[1] = c2[1]-c0[1];
  edge1[2] = c2[2]-c0[2];

  T normal[3];
  normal[0] = edge0[1]*edge1[2] - edge0[2]*edge1[1];
  normal[1] = edge0[2]*edge1[0] - edge0[0]*edge1[2];
  normal[2] = edge0[0]*edge1[1] - edge0[1]*edge1[0];

  T normalLen = stk_sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  normal[0] /= normalLen;
  normal[1] /= normalLen;
  normal[2] /= normalLen;

  T p[3];
  p[0] = centroid[0]-c0[0];
  p[1] = centroid[1]-c0[1]; 
  p[2] = centroid[2]-c0[2];

  return 2.0*(normal[0]*p[0] + normal[1]*p[1] + normal[2]*p[2]);
}




}


#endif
