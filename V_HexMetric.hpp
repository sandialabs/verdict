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
 * V_HexMetric.hpp contains declarations of hex related shape quantities
 *
 * This file is part of VERDICT
 *
 */

#ifndef VERDICT_HEX_METRIC_HPP
#define VERDICT_HEX_METRIC_HPP

#include <assert.h>

extern "C" double hex_scaled_jacobian( int num_nodes, double coordinates[][3]);
extern "C" double hex_shape          ( int num_nodes, double coordinates[][3]);
extern "C" double hex_nodal_jacobian_ratio( int num_nodes, double coordinates[][3]);


#endif
