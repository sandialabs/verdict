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
  /*
    inline double stk_sqrt(double x) {
      return std::sqrt(x);
    }
    template <typename T>
    inline T stk_sqrt(T x) {
      return stk::math::sqrt(x);
    }
    */









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



double tet_dimension( int num_nodes, double coordinates[][3] ) {
 
    VerdictVector c0(coordinates[0][0], coordinates[0][1], coordinates[0][2]);
    VerdictVector c1(coordinates[1][0], coordinates[1][1], coordinates[1][2]);
    VerdictVector c2(coordinates[2][0], coordinates[2][1], coordinates[2][2]);
    VerdictVector c3(coordinates[3][0], coordinates[3][1], coordinates[3][2]);

    VerdictVector centroid = 0.25 * (c0 + c1 + c2 + c3);    
    
    VerdictVector edge0 = c1 - c0;
    VerdictVector edge1 = c2 - c0;
    VerdictVector normal = edge0*edge1;
    normal.normalize();
    VerdictVector p = centroid - c0;
    double radius0 = normal%p; 

    edge0 = c3 - c0;
    edge1 = c1 - c0;
    normal = edge0*edge1;
    normal.normalize();
    p = centroid - c0;    
    double radius1 = normal%p; 

    edge0 = c2 - c0;
    edge1 = c3 - c0;
    normal = edge0*edge1;
    normal.normalize();
    p = centroid - c0;    
    double radius2 = normal%p; 

    edge0 = c3 - c1;
    edge1 = c2 - c1;
    normal = edge0*edge1;
    normal.normalize();
    p = centroid - c1;
    double radius3 = normal%p; 

    return 0.5*(radius0 + radius1 + radius2 + radius3 );
}




}


#endif
