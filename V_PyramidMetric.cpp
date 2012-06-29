/*=========================================================================

  Module:    V_PyramidMetric.cpp

  Copyright (c) 2006 Sandia Corporation.
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/


/*
 *
 * PyramidMetrics.cpp contains quality calculations for Pyramids
 *
 * This file is part of VERDICT
 *
 */


#include "verdict.h"
#include "VerdictVector.hpp"
#include "verdict.h"
#include <memory.h> 

// local method
void v_make_pyramid_tets(double coordinates[][3], double tet1_coords[][3], double tet2_coords[][3],
                                                  double tet3_coords[][3], double tet4_coords[][3]);
void v_make_pyramid_faces(double coordinates[][3], double base[][3], double tri1[][3],
                          double tri2[][3], double tri3[][3], double tri4[][3]);
/*
  the pyramid element

       5
       ^
       |\ 
      /| \\_
     |  \   \
     |  | \_ \_
     /   \/4\  \
    |   /|    \ \_
    | /  \      \ \
    /     |       \
  1 \_    |      _/3
      \_   \   _/
        \_ | _/
          \_/
          2

    a quadrilateral base and a pointy peak like a pyramid
          
*/


/*!
  the volume of a pyramid

  the volume is calculated by dividing the pyramid into
  2 tets and summing the volumes of the 2 tets.
*/

C_FUNC_DEF double v_pyramid_volume( int num_nodes, double coordinates[][3] )
{
    
  double volume = 0;
  VerdictVector side1, side2, side3;
  
  if (num_nodes == 5)
  {
    // divide the pyramid into 2 tets and calculate each

    side1.set( coordinates[1][0] - coordinates[0][0],
        coordinates[1][1] - coordinates[0][1],
        coordinates[1][2] - coordinates[0][2] );
    
    side2.set( coordinates[3][0] - coordinates[0][0],
        coordinates[3][1] - coordinates[0][1],
        coordinates[3][2] - coordinates[0][2] );
    
    side3.set( coordinates[4][0] - coordinates[0][0],
        coordinates[4][1] - coordinates[0][1], 
        coordinates[4][2] - coordinates[0][2] );
    
    // volume of the first tet
    volume = (side3 % (side1 * side2 ))/6.0;
    
    
    side1.set( coordinates[3][0] - coordinates[2][0],
        coordinates[3][1] - coordinates[2][1],
        coordinates[3][2] - coordinates[2][2] );
    
    side2.set( coordinates[1][0] - coordinates[2][0],
        coordinates[1][1] - coordinates[2][1],
        coordinates[1][2] - coordinates[2][2] );
    
    side3.set( coordinates[4][0] - coordinates[2][0],
        coordinates[4][1] - coordinates[2][1],
        coordinates[4][2] - coordinates[2][2] );
    
    // volume of the second tet
    volume += (side3 % (side1 * side2 ))/6.0;
 
  }   
  return (double)volume;
    
}

C_FUNC_DEF double v_pyramid_jacobian( int num_nodes, double coordinates[][3] )
{
  // break the pyramid into four tets return the minimum scaled jacobian of the two tets
  double tet1_coords[4][3];
  double tet2_coords[4][3];
  double tet3_coords[4][3];
  double tet4_coords[4][3];

  v_make_pyramid_tets(coordinates, tet1_coords, tet2_coords, tet3_coords, tet4_coords);

  double j1 = v_tet_jacobian(4, tet1_coords);
  double j2 = v_tet_jacobian(4, tet2_coords);
  double j3 = v_tet_jacobian(4, tet3_coords);
  double j4 = v_tet_jacobian(4, tet4_coords);

  double p1 = j1 < j2 ? j1 : j2;
  double p2 = j3 < j4 ? j3 : j4;

  return p1 < p2 ? p1 : p2;
}

C_FUNC_DEF double v_pyramid_scaled_jacobian( int num_nodes, double coordinates[][3] )
{
  // ideally there will be four equilateral triangles and one square.  
  // Test each face
  double base[4][3];
  double tri1[3][3];
  double tri2[3][3];
  double tri3[3][3];
  double tri4[3][3];

  v_make_pyramid_faces(coordinates, base, tri1, tri2, tri3, tri4);

  double s1 = v_quad_scaled_jacobian(4, base);
  double s2 = v_tri_scaled_jacobian(3, tri1);
  double s3 = v_tri_scaled_jacobian(3, tri2);
  double s4 = v_tri_scaled_jacobian(3, tri3);
  double s5 = v_tri_scaled_jacobian(3, tri4);

  return .2*(s1 + s2 + s3 + s4 + s5);
}

C_FUNC_DEF double v_pyramid_shape( int num_nodes, double coordinates[][3] )
{
  // ideally there will be four equilateral triangles and one square.  
  // Test each face
  double base[4][3];
  double tri1[3][3];
  double tri2[3][3];
  double tri3[3][3];
  double tri4[3][3];

  v_make_pyramid_faces(coordinates, base, tri1, tri2, tri3, tri4);

  double s1 = v_quad_shape(4, base);
  double s2 = v_tri_shape(3, tri1);
  double s3 = v_tri_shape(3, tri2);
  double s4 = v_tri_shape(3, tri3);
  double s5 = v_tri_shape(3, tri4);

  return .2*(s1 + s2 + s3 + s4 + s5);
}

void v_make_pyramid_tets(double coordinates[][3], double tet1_coords[][3], double tet2_coords[][3],
                                                  double tet3_coords[][3], double tet4_coords[][3])
{
   // tet1
  tet1_coords[0][0] = coordinates[0][0];
  tet1_coords[0][1] = coordinates[0][1];
  tet1_coords[0][2] = coordinates[0][2];

  tet1_coords[1][0] = coordinates[1][0];
  tet1_coords[1][1] = coordinates[1][1];
  tet1_coords[1][2] = coordinates[1][2];

  tet1_coords[2][0] = coordinates[2][0];
  tet1_coords[2][1] = coordinates[2][1];
  tet1_coords[2][2] = coordinates[2][2];

  tet1_coords[3][0] = coordinates[4][0];
  tet1_coords[3][1] = coordinates[4][1];
  tet1_coords[3][2] = coordinates[4][2];

  // tet2
  tet2_coords[0][0] = coordinates[0][0];
  tet2_coords[0][1] = coordinates[0][1];
  tet2_coords[0][2] = coordinates[0][2];

  tet2_coords[1][0] = coordinates[2][0];
  tet2_coords[1][1] = coordinates[2][1];
  tet2_coords[1][2] = coordinates[2][2];

  tet2_coords[2][0] = coordinates[3][0];
  tet2_coords[2][1] = coordinates[3][1];
  tet2_coords[2][2] = coordinates[3][2];

  tet2_coords[3][0] = coordinates[4][0];
  tet2_coords[3][1] = coordinates[4][1];
  tet2_coords[3][2] = coordinates[4][2];

  // tet3
  tet3_coords[0][0] = coordinates[0][0];
  tet3_coords[0][1] = coordinates[0][1];
  tet3_coords[0][2] = coordinates[0][2];

  tet3_coords[1][0] = coordinates[1][0];
  tet3_coords[1][1] = coordinates[1][1];
  tet3_coords[1][2] = coordinates[1][2];

  tet3_coords[2][0] = coordinates[3][0];
  tet3_coords[2][1] = coordinates[3][1];
  tet3_coords[2][2] = coordinates[3][2];

  tet3_coords[3][0] = coordinates[4][0];
  tet3_coords[3][1] = coordinates[4][1];
  tet3_coords[3][2] = coordinates[4][2];

  // tet4
  tet4_coords[0][0] = coordinates[1][0];
  tet4_coords[0][1] = coordinates[1][1];
  tet4_coords[0][2] = coordinates[1][2];

  tet4_coords[1][0] = coordinates[2][0];
  tet4_coords[1][1] = coordinates[2][1];
  tet4_coords[1][2] = coordinates[2][2];

  tet4_coords[2][0] = coordinates[3][0];
  tet4_coords[2][1] = coordinates[3][1];
  tet4_coords[2][2] = coordinates[3][2];

  tet4_coords[3][0] = coordinates[4][0];
  tet4_coords[3][1] = coordinates[4][1];
  tet4_coords[3][2] = coordinates[4][2];
}

void v_make_pyramid_faces(double coordinates[][3], double base[][3], double tri1[][3],
                          double tri2[][3], double tri3[][3], double tri4[][3])
{
   // base
  base[0][0] = coordinates[0][0];
  base[0][1] = coordinates[0][1];
  base[0][2] = coordinates[0][2];

  base[1][0] = coordinates[1][0];
  base[1][1] = coordinates[1][1];
  base[1][2] = coordinates[1][2];

  base[2][0] = coordinates[2][0];
  base[2][1] = coordinates[2][1];
  base[2][2] = coordinates[2][2];

  base[3][0] = coordinates[3][0];
  base[3][1] = coordinates[3][1];
  base[3][2] = coordinates[3][2];

  // tri1
  tri1[0][0] = coordinates[0][0];
  tri1[0][1] = coordinates[0][1];
  tri1[0][2] = coordinates[0][2];

  tri1[1][0] = coordinates[1][0];
  tri1[1][1] = coordinates[1][1];
  tri1[1][2] = coordinates[1][2];

  tri1[2][0] = coordinates[4][0];
  tri1[2][1] = coordinates[4][1];
  tri1[2][2] = coordinates[4][2];

  // tri2
  tri2[0][0] = coordinates[1][0];
  tri2[0][1] = coordinates[1][1];
  tri2[0][2] = coordinates[1][2];

  tri2[1][0] = coordinates[2][0];
  tri2[1][1] = coordinates[2][1];
  tri2[1][2] = coordinates[2][2];

  tri2[2][0] = coordinates[4][0];
  tri2[2][1] = coordinates[4][1];
  tri2[2][2] = coordinates[4][2];

  // tri3
  tri3[0][0] = coordinates[2][0];
  tri3[0][1] = coordinates[2][1];
  tri3[0][2] = coordinates[2][2];

  tri3[1][0] = coordinates[3][0];
  tri3[1][1] = coordinates[3][1];
  tri3[1][2] = coordinates[3][2];

  tri3[2][0] = coordinates[4][0];
  tri3[2][1] = coordinates[4][1];
  tri3[2][2] = coordinates[4][2];

   // tri4
  tri4[0][0] = coordinates[3][0];
  tri4[0][1] = coordinates[3][1];
  tri4[0][2] = coordinates[3][2];

  tri4[1][0] = coordinates[0][0];
  tri4[1][1] = coordinates[0][1];
  tri4[1][2] = coordinates[0][2];

  tri4[2][0] = coordinates[4][0];
  tri4[2][1] = coordinates[4][1];
  tri4[2][2] = coordinates[4][2];
}


C_FUNC_DEF void v_pyramid_quality( int num_nodes, double coordinates[][3], 
    unsigned int metrics_request_flag, PyramidMetricVals *metric_vals )
{
  memset( metric_vals, 0, sizeof( PyramidMetricVals ) );

  if(metrics_request_flag & V_PYRAMID_VOLUME)
    metric_vals->volume = v_pyramid_volume(num_nodes, coordinates);
  else if(metrics_request_flag & V_PYRAMID_JACOBIAN)
    metric_vals->jacobian = v_pyramid_jacobian(num_nodes, coordinates);
  else if(metrics_request_flag & V_PYRAMID_SCALED_JACOBIAN)
    metric_vals->scaled_jacobian = v_pyramid_scaled_jacobian(num_nodes, coordinates);
  else if(metrics_request_flag & V_PYRAMID_SHAPE)
    metric_vals->shape = v_pyramid_shape(num_nodes, coordinates);
}

