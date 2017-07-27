/*!
 * \brief Unittests for verdict element quality metrics
 *
 * \date 7/26/2017
 *
 * \author Matt Staten
 */

#include "gtest/gtest.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>

#include <verdict.h>

#define MAX_NODES_PER_ELEMENT 27
#define VERDICT_SIGNIFICANT_FIG 5    // 7 significant figures for doubles

struct metric_and_answer
{
  VerdictFunction mFunction;
  double mAnswer;
};

struct test_case
{
  const char* mTestName;
  std::vector<metric_and_answer> mTestData;
  int mNumNodes;
  double mCoords[MAX_NODES_PER_ELEMENT][3];
};

void runtest(test_case &this_testcase)
{
  const char *testname = this_testcase.mTestName;
  int num_nodes = this_testcase.mNumNodes;

  size_t itest = 0;
  for (metric_and_answer & data : this_testcase.mTestData)
  {
    itest++;
    // Compute the answer from the verdict library
    double answer_from_lib = data.mFunction(num_nodes, this_testcase.mCoords);

    std::stringstream expected;
    expected.setf(std::ios::scientific, std::ios::floatfield);
    expected.precision(VERDICT_SIGNIFICANT_FIG);
    expected << data.mAnswer;

    std::stringstream answer;
    answer.setf(std::ios::scientific, std::ios::floatfield);
    answer.precision(VERDICT_SIGNIFICANT_FIG);
    answer << answer_from_lib;

    EXPECT_EQ(expected.str(), answer.str());
    if (expected.str() != answer.str())
    {
      std::cout << std::endl;
      std::cout << "Test case \"" << testname
        << "\" #" << itest << " FAILED" << std::endl;

      std::cout << "answer calculated was    "
        << answer.str() << std::endl
        << "answer expected was      "
        << expected.str()
        << std::endl << std::endl;
      EXPECT_TRUE(false);
    }
    else
    {
      std::cout << "Test case \"" << testname
        << "\" #" << itest << " passed" << std::endl;
    }
  }
}

TEST(verdict, edge_calc_1)
{
  test_case testcase = {
    "edge_calc_1",
    { { v_edge_length, 1.732050807568877 } },
    2,
    {
      { 0, 0, 0 },
      { 1, 1, 1 }
    }
  };

  runtest(testcase);
}

TEST(verdict, edge_calc_2)
{
  test_case testcase = {
    "edge_calc_2",
    { { v_edge_length, 1.0 } },
    2,
    {
      { 0, 0, 0 },
      { 1, 0, 0 }
    }
  };

  runtest(testcase);
}


TEST(verdict, edge_calc_3)
{
  test_case testcase = {
    "edge_calc_3",
    { { v_edge_length, 0.0 } },
    2,
    {
      { 0, 0, 0 },
      { 0, 0, 0 }
    }
  };

  runtest(testcase);
}

TEST(verdict, simple_wedge_volume)
{
  test_case testcase = {
    "simple_wedge_volume",
    { { v_wedge_volume, 0.5 } },
    6,
    {
      { 0, 0, 0 },
      { -1, 1, 0 },
      { -1, 0, 0 },
      { 0, 0, 1 },
      { -1, 1, 1 },
      { -1, 0, 1 }
    }
  };

  runtest(testcase);
}

TEST(verdict, singularity_wedge)
{
  test_case testcase = {
    "singularity wedge",
    { { v_wedge_volume, 0.0 } },
    6,
    {
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 } }
  };

  runtest(testcase);
}

TEST(verdict, simple_tri)
{
  test_case testcase = {
    "simple_tri",
    {
      { v_tri_area, 10.825317 },
      { v_tri_minimum_angle, 60 },
      { v_tri_maximum_angle, 60 },
      { v_tri_condition, 1.0 },
      { v_tri_shape, 1.0 },
      //    { v_tri_shape_and_size, 0.00853333 },
      { v_tri_distortion, 1.0 }
    },
    3,
    {
      { 0, 0, 0 },
      { 5, 0, 0 },
      { 2.5, 4.330127, 0 }
      }
  };

  runtest(testcase);
}

TEST(verdict, singular_tri)
{
  test_case testcase = {
    "singular_tri",
    {
      { v_tri_area, 0.433013 },
      { v_tri_aspect_ratio, 1 },
      { v_tri_condition, 1 },
      { v_tri_distortion, 1 },
      { v_tri_minimum_angle, 60 },
      { v_tri_maximum_angle, 60 },
      //       {v_tri_relative_size_squared, 0.1875},
      //      { v_tri_shape_and_size, 0.1875 },
      { v_tri_shape, 1 }
    },
    3,
    {
      { 0, 0, 0 },
      { 0.5, 0.8660254037, 0 },
      { 1, 0, 0 } }
  };

  runtest(testcase);
}

TEST(verdict, simple_quad1)
{
  test_case testcase = {
    "simple_quad1",
    { { v_quad_skew, 0 } },
    4,
    {
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 7, 0 },
      { 0, 7, 0 }
    }
  };

  runtest(testcase);
}

TEST(verdict, simple_quad2)
{
  test_case testcase = {
    "simple_quad2",
    {
      { v_quad_aspect_ratio, 1.42996 },
      { v_quad_skew, 0.09245 },
      { v_quad_taper, 0.745356 },
      { v_quad_warpage, 0.008 },
      { v_quad_area, 2.69258 },
      { v_quad_stretch, 0.57735 },
      { v_quad_minimum_angle, 56.7891 },
      { v_quad_maximum_angle, 90 },
      { v_quad_condition, 2.30793 },
      { v_quad_jacobian, 1.11417 },
      { v_quad_shear, 0.557086 },
      { v_quad_shape, 0.433289 },
      //{ v_quad_shape_and_size,  0.059764   },
      //{ v_quad_shear_and_size,  0.0768395    },
      { v_quad_distortion, 0.56268 }
    },
    4,
    {
      { 2, 0, 0 },    //1
      { 1, 1, 2 },   //2
      { 0, 1, 0 },      //3
      { 0, 0, 0 },   //0
      }
  };

  runtest(testcase);
}

TEST(verdict, tet_test1)
{
  test_case testcase = {
    "tet_test1",
    {
      { v_tet_volume, 166.66666 },
      { v_tet_condition, 1.22474 },
      { v_tet_jacobian, 1000 },
      { v_tet_shape, 0.839947 },
      //{v_tet_shape_and_size,  0.0000302381  },
      { v_tet_distortion, 1 }
    },
    4,
    {
      { -5, -5, -5 },
      { -5, 5, -5 },
      { -5, -5, 5 },
      { 5, -5, -5 },
      }
  };

  runtest(testcase);
}

TEST(verdict, hex_test1)
{
  test_case testcase = {
    "hex_test1",
    {
      { v_hex_skew, 0.24589 },
      { v_hex_taper, 0.178458 },
      { v_hex_volume, 0.815667 },
      { v_hex_stretch, 0.62097 },
      { v_hex_diagonal, 0.689622 },
      { v_hex_dimension, 0.524594 },
      { v_hex_condition, 1.27306 },
      { v_hex_jacobian, 0.477 },
      { v_hex_shear, 0.77785 },
      { v_hex_shape, 0.789785 },
      //{ v_hex_shear_and_size, 0.524143 },
      //{ v_hex_shape_and_size, 0.532186 },
      { v_hex_distortion, 0.584798 }
    },
    8,
    {
      { -.2, -.7, -.3 },  //1
      { -.7, .4, -.6 },   //2
      { -.5, .5, .3 },  //3
      { -.3, -.5, .5 }, //0
      { .5, -.8, -.2 },   //5
      { .4, .4, -.6 },    //6
      { .2, .5, .2 },   //7 
      { .5, -.3, .8 }  //4
      }

  };

  runtest(testcase);
}

TEST(verdict, hex27_test1)
{
  test_case testcase = {
    "hex27_test1",
    { { v_hex_jacobian, 0.560215 } },
    27,
    {
      { -0.5, -0.5, 0.5 },
      { -0.5, -0.5, -0.5 },
      { -0.5, 0.5, -0.5 },
      { -0.5, 0.5, 0.5 },
      { 0.5, -0.5, 0.5 },
      { 0.5, -0.5, -0.5 },
      { 0.5, 0.5, -0.5 },
      { 0.5, 0.5, 0.5 },
      { -0.59093, -0.554599, 0.031478 },
      { -0.551216, -0.052901, -0.59258 },
      { -0.589451, 0.534972, -0.00525 },
      { -0.575503, -0.008395, 0.604191 },
      { 0.035399, -0.527937, 0.582587 },
      { -0.015396, -0.537798, -0.606826 },
      { -0.027222, 0.51188, -0.576597 },
      { -0.006207, 0.565371, 0.585294 },
      { 0.562112, -0.548542, -0.000139 },
      { 0.570746, 0.008609, -0.595057 },
      { 0.563111, 0.550215, -0.024469 },
      { 0.537894, 0.050161, 0.606614 },
      { -0.02166575, -0.0022401, 0.0023115 },
      { -0.610022, -0.053727, -0.007623 },
      { 0.593828, -0.007149, -0.019913 },
      { 0.019742, 0.031299, 0.64416 },
      { 0.033398, -0.032358, -0.596668 },
      { 0.021178, -0.581132, 0.013062 },
      { 0.000338, 0.60863, -0.017854 }
    }
  };

  runtest(testcase);
}

TEST(verdict, wedge_21_test1)
{
  test_case testcase = {
    "wedge21_test1",
    { { v_wedge_jacobian, 0.999998 } },
    21,
    {
      { 0, 0, -1 },
      { 1, 0, -1 },
      { 0, 1, -1 },
      { 0, 0, 1 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { 0.5, 0, -1 },
      { 0.5, 0.5, -1 },
      { 0, 0.5, -1 },
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0.5, 0, 1 },
      { 0.5, 0.5, 1 },
      { 0, 0.5, 1 },
      { .333333, .333333, 0 },
      { .333333, .333333, -1 },
      { .333333, .333333, 1 },
      { .5, .5, 0 },
      { 0, .5, 0 },
      { 0.5, 0, 0 }
    }
  };

  runtest(testcase);
}

TEST(verdict, tet15_test1)
{
  test_case testcase = {
    "tet15_test1",
    { { v_tet_jacobian, 2.0 } },
    15,
    {
      { 0, 0, 0 },
      { 0, 2, 0 },
      { 1, 1, 0 },
      { 6.12323399573677e-17, 1, -1 },
      { 0, 1, 0 },
      { 0.5, 1.5, 0 },
      { 0.5, 0.5, 0 },
      { 3.06161699786838e-17, 0.5, -0.5 },
      { 3.06161699786838e-17, 1.5, -0.5 },
      { 0.707106781186548, 1, -0.707106781186547 },
      { 0.301776695296637, 1, -0.301776695296637 },
      { 0.333333333333333, 1, 0 },
      { 0.425380791638466, 1.33333333333333, -0.425380791638466 },
      { 0.425380791638466, 0.666666666666667, -0.425380791638466 },
      { 2.04107799857892e-17, 1, -0.333333333333333 }
    }
  };

  runtest(testcase);
}

