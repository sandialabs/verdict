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
#include "verdict.h"
#include "VerdictVector.hpp"


namespace verdict {


  double tet_dimension( int num_nodes, double coordinates[][3] ) { 
    
    //area1 (0,1,2)
    double a1 = verdict::tri_area(3, coordinates);
    
    //area2 (0,3,1)
    double tmp_coords[3][3];
    tmp_coords[0][0] = coordinates[0][0];
    tmp_coords[0][1] = coordinates[0][1];
    tmp_coords[0][2] = coordinates[0][2];

    tmp_coords[1][0] = coordinates[3][0];
    tmp_coords[1][1] = coordinates[3][1];
    tmp_coords[1][2] = coordinates[3][2];

    tmp_coords[2][0] = coordinates[1][0];
    tmp_coords[2][1] = coordinates[1][1];
    tmp_coords[2][2] = coordinates[1][2];

    double a2 = verdict::tri_area(3, tmp_coords);

    //area3 (0,2,3)
    tmp_coords[1][0] = coordinates[2][0];
    tmp_coords[1][1] = coordinates[2][1];
    tmp_coords[1][2] = coordinates[2][2];

    tmp_coords[2][0] = coordinates[3][0];
    tmp_coords[2][1] = coordinates[3][1];
    tmp_coords[2][2] = coordinates[3][2];
    
    double a3 = verdict::tri_area(3, tmp_coords);

    //area4 (1,3,2)
    tmp_coords[0][0] = coordinates[1][0];
    tmp_coords[0][1] = coordinates[1][1];
    tmp_coords[0][2] = coordinates[1][2];

    tmp_coords[1][0] = coordinates[3][0];
    tmp_coords[1][1] = coordinates[3][1];
    tmp_coords[1][2] = coordinates[3][2];

    tmp_coords[2][0] = coordinates[2][0];
    tmp_coords[2][1] = coordinates[2][1];
    tmp_coords[2][2] = coordinates[2][2];
    
    double a4 = verdict::tri_area(3, tmp_coords);

    double tet_volume = verdict::tet_volume(4, coordinates);

    double inradius = 3*tet_volume/(a1 + a2 + a3 + a4);

    return inradius;
}
  


}


#endif
