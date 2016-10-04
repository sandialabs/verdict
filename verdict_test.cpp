/*=========================================================================

Module:    verdict_test.cpp

Copyright (c) 2006 Sandia Corporation.
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*
 *
 * verdict_test.cpp provides routines for testing the quality metrics code
 *
 * This file is part of VERDICT
 *
 */
#define VERDICT_EXPORTS

#define VERDICT_THREAD_SAFE 0

#include "verdict.h"
#include "v_vector.h"
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>


#define MAX_NODES_PER_ELEMENT 27
#define MAX_TESTS_PER_ELEMENT 20


#define VERDICT_SIGNIFICANT_FIG 5    // 7 significant figures for doubles

struct test_case
{
  const char* testname;
  VerdictFunction function[MAX_TESTS_PER_ELEMENT];
  int num_nodes;
  // note: the 1st dim. of coords must bigger than the maximum num_nodes 
  // for any one element being tested
  double coords[MAX_NODES_PER_ELEMENT][3];
  double answer[MAX_TESTS_PER_ELEMENT];
};


int main( )
{
  // all test cases go here  
  test_case testcases[] = {
  {
  "edge calc 1",
  {v_edge_length, 0},
  2,
  { { 0,0,0 }, {1,1,1} },
  { 1.732050807568877, 0 }
  },

  {
  "edge calc 2",
  {v_edge_length, 0 },
  2,
  { { 0,0,0 }, { 1,0,0 } },
  { 1.0, 0 }
  },

  {
  "edge calc 3",
  {v_edge_length, 0 },
  2,
  { { 0,0,0 }, { 0,0,0 } },
  { 0, 0 }
  },

  { 
  "simple wedge" ,
  {v_wedge_volume, 0},
  6,
  { { 0,0,0}, {-1,1,0}, {-1,0,0}, {0,0,1}, {-1,1,1}, {-1,0,1} },
  { 0.5, 0 }
  },

  {
  "singularity wedge",
  {v_wedge_volume, 0},
  6,
  { {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0} },
  { 0 }
  },

  {
  "simple tri",
  { v_tri_area, v_tri_minimum_angle, v_tri_maximum_angle, v_tri_condition, 
  v_tri_shape, v_tri_shape_and_size, v_tri_distortion, 0},
  3,
  { 
  {0,0,0}, 
  {5,0,0}, 
  {2.5, 4.330127, 0}, 
  },

  { 10.825317,60,60,1.0,1.0,.00853333,1.0}
  },
    {
      "singular tri",
      {v_tri_area, v_tri_aspect_ratio, v_tri_condition, 
       v_tri_distortion, v_tri_minimum_angle, v_tri_maximum_angle,
       v_tri_relative_size_squared, v_tri_shape, v_tri_shape_and_size,
       0}, 
      3,
      { 
        {0,0,0}, 
        {0.5,0.8660254037,0}, 
        {1,0,0} },
      { .433013, 1, 1, 1, 60, 60, .1875, 1, .1875,0}
    },
  {
  "simple quad",
  { v_quad_skew, 0},
  4,
  {
  {0,0,0},
  {1,0,0},
  {1,7,0}, 
  {0,7,0 } 
  }, 
  { 0, 0 }
  },

  {
  "simple quad",
  { v_quad_aspect_ratio, v_quad_skew, v_quad_taper, v_quad_warpage, v_quad_area, 
  v_quad_stretch, v_quad_minimum_angle, v_quad_maximum_angle,
  v_quad_condition, v_quad_jacobian, v_quad_shear,
  v_quad_shape, v_quad_shape_and_size, v_quad_shear_and_size,
  v_quad_distortion, 0},
  4,
  {
  {2,0,0},    //1
  {1,1,2},   //2
  {0,1,0 },      //3
  {0,0,0},   //0



  }, 
  { 1.42996, .09245, .745356, .008, 2.69258, .57735, 56.7891, 90, 2.30793, 1.11417 , .557086, .433289, .059764, .0768395, .56268, 0}
  },
  {
  "tet test",
  { v_tet_volume, v_tet_condition, v_tet_jacobian, 
  v_tet_shape, v_tet_shape_and_size, v_tet_distortion, 0 },
  4,
  {
  {-5, -5, -5 },
  {-5, 5, -5 },
  {-5, -5, 5 },
  {5, -5, -5 },
 
  },

  {166.66666,1.22474,1000,.839947,.0000302381,1,0}
  },
  {
  "hex test",
  { v_hex_skew, v_hex_taper, v_hex_volume, v_hex_stretch, v_hex_diagonal,
  v_hex_dimension, v_hex_condition, v_hex_jacobian, v_hex_shear,
  v_hex_shape, v_hex_shear_and_size, v_hex_shape_and_size,
  v_hex_distortion, 0 },
  8,

  { 
  { -.2, -.7, -.3 },  //1
  { -.7, .4, -.6 },   //2
  { -.5, .5, .3  },  //3
  { -.3, -.5, .5  }, //0
  
  { .5, -.8, -.2 },   //5
  { .4, .4, -.6  },    //6
  { .2, .5, .2   },   //7 
  { .5, -.3, .8  }  //4
  },


  {0.24589,0.178458,0.813062,0.62097,0.689622,0.524594,1.27306,.477,0.77785,0.789785,0.524143,0.532186,.584798}
  },

  {
  "hex27 test",
  { v_hex_jacobian, 0 },
  27,
  {
  {-0.5,-0.5, 0.5},
  {-0.5,-0.5,-0.5},
  {-0.5,0.5, -0.5},
  {-0.5, 0.5,0.5},
  {0.5,-0.5,0.5},
  {0.5,-0.5,-0.5},
  {0.5,0.5,-0.5},
  {0.5,0.5,0.5},
  {-0.59093,-0.554599,0.031478},
  {-0.551216,-0.052901, -0.59258},
  {-0.589451,0.534972,-0.00525},
  {-0.575503,-0.008395, 0.604191},
  {0.035399,-0.527937,0.582587},
  {-0.015396,-0.537798,-0.606826},
  {-0.027222,0.51188,-0.576597},
  {-0.006207, 0.565371,0.585294},
  {0.562112,-0.548542, -0.000139},
  {0.570746,0.008609,-0.595057},
  {0.563111,0.550215,-0.024469},
  {0.537894,0.050161, 0.606614},
  {-0.02166575,-0.0022401,0.0023115},
  {-0.610022, -0.053727,-0.007623},
  {0.593828,-0.007149,-0.019913},
  {0.019742,0.031299,0.64416},
  {0.033398,-0.032358, -0.596668},
  {0.021178,-0.581132,0.013062},
  {0.000338, 0.60863,-0.017854}
  },
  {0.560215}
  },


  {
  "wedge21 test",
  { v_wedge_jacobian, 0 },
  21,
  {
    {0, 0, -1},
    {1, 0, -1},
    {0, 1, -1},
    {0, 0, 1},
    {1, 0, 1},
    {0, 1, 1},
    {0.5, 0, -1},
    {0.5, 0.5, -1},
    {0, 0.5, -1},
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {0.5, 0, 1},
    {0.5, 0.5, 1},
    {0, 0.5, 1},
    {.333333, .333333, 0},
    {.333333, .333333, -1},
    {.333333, .333333, 1},
    {.5, .5, 0},
    {0, .5, 0},
    {0.5, 0, 0}
  },
  {0.999998}
  },

  {
  "tet15 test",
  { v_tet_jacobian, 0 },
  15,
  {
    {0, 0,0},
    {0, 2,0},
    {1, 1,0},
    {6.12323399573677e-17, 1,-1},
    {0, 1,0},
    {0.5,1.5 ,0},
    {0.5, 0.5,0},
    {3.06161699786838e-17, 0.5,-0.5},
    {3.06161699786838e-17, 1.5,-0.5},
    {0.707106781186548, 1,-0.707106781186547},
    {0.301776695296637, 1,-0.301776695296637},
    {0.333333333333333, 1,0},
    {0.425380791638466, 1.33333333333333,-0.425380791638466},
    {0.425380791638466, 0.666666666666667,-0.425380791638466},
    {2.04107799857892e-17, 1,-0.333333333333333}
  },
  {2.0}
  },


    // keep this one last
    { 0, {0} , 0, {{0}} , {0} } };



  v_set_tri_size(1);
  v_set_quad_size(1);
  v_set_tet_size(1);
  v_set_hex_size(1);
      
    
  bool passed = true; // have all the tests performed so far passed?

  std::cout.setf( std::ios::scientific, std::ios::floatfield );
  std::cout.precision(VERDICT_SIGNIFICANT_FIG);

  std::cout << std::endl;

   
  // loop through each test
  for (int i = 0; testcases[i].testname != 0; i++ )
    {
     
    for (int j = 0;  testcases[i].function[j] != 0; j++ )
      {       
      double answer_from_lib =
        (testcases[i].function[j])
        (testcases[i].num_nodes, testcases[i].coords);
       
      std::stringstream expected;
      expected.setf(std::ios::scientific, std::ios::floatfield);
      expected.precision(VERDICT_SIGNIFICANT_FIG);
      expected << testcases[i].answer[j];

      std::stringstream answer;
      answer.setf(std::ios::scientific, std::ios::floatfield);
      answer.precision(VERDICT_SIGNIFICANT_FIG);
      answer << answer_from_lib;

      if ( expected.str() != answer.str() )
        {
        std::cout << std::endl;
        std::cout << "Test case \"" << testcases[i].testname
             << "\" #" << j+1 << " FAILED" << std::endl;

        std::cout    << "answer calculated was    "
                << answer.str() << std::endl
                << "answer expected was      " 
                << expected.str()
                << std::endl << std::endl;
        passed = false;
        }
      else
        {
        std::cout << "Test case \"" << testcases[i].testname
             << "\" #" << j+1 << " passed" << std::endl;
//           << "answer calculated was    " 
//                << answer_from_lib << endl
//           << "answer expected was      " 
//           << testcases[i].answer[j] << endl;
//    //     cout << endl;
        }
      }
    }

  std::cout << std::endl;

  return passed ? 0 : 1;
}
