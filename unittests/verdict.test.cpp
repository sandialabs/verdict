/*=========================================================================

  Module:    V_EdgeMetric.cpp

  Copyright 2003,2006,2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
  Under the terms of Contract DE-NA0003525 with NTESS,
  the U.S. Government retains certain rights in this software.

  See LICENSE for details.

=========================================================================*/

/*!
 * \brief Unittests for verdict element quality metrics
 *
 * \date 7/26/2017
 *
 * \author Matt Staten
 */

#include "verdict.h"

#include "gtest/gtest.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#define MAX_NODES_PER_ELEMENT 27
// 7 significant figures for doubles, use 5 to make tests less sensitive
// #define VERDICT_SIGNIFICANT_FIG 5
// absolute and relative tolerances for comparing floating point answers
#define VERDICT_ABSOLUTE_TOL 1e-05
#define VERDICT_RELATIVE_TOL 1e-05

struct metric_and_answer
{
  std::function<double(int, double[][3])> mFunction;
  double mAnswer;
};

struct test_case
{
  const char* mTestName;
  std::vector<metric_and_answer> mTestData;
  int mNumNodes;
  double mCoords[MAX_NODES_PER_ELEMENT][3];
};

void runtest(test_case& this_testcase)
{
  const char* testname = this_testcase.mTestName;
  int num_nodes = this_testcase.mNumNodes;

  size_t itest = 0;
  for (metric_and_answer& data : this_testcase.mTestData)
  {
    itest++;

    // Compute the answer from the verdict library
    const double calculated_answer = data.mFunction(num_nodes, this_testcase.mCoords);
    const double expected_answer = data.mAnswer;

    // list the current subcase, otherwise there is no way to tell which metric failed
    std::cout << "Test case \"" << testname << "\" #" << itest
              << " Calculated = " << calculated_answer << " Expected = " << expected_answer
              << std::endl;

    EXPECT_NEAR(calculated_answer, expected_answer,
      std::abs(expected_answer) * VERDICT_RELATIVE_TOL + VERDICT_ABSOLUTE_TOL);

    // old test using strings
    /*
    std::stringstream expected;
    expected.setf(std::ios::scientific, std::ios::floatfield);
    expected.precision(VERDICT_SIGNIFICANT_FIG);
    expected << data.mAnswer;

    std::stringstream answer;
    answer.setf(std::ios::scientific, std::ios::floatfield);
    answer.precision(VERDICT_SIGNIFICANT_FIG);
    answer << answer_from_lib;

    // for expected values of zero, we only care about absolute tolerance, not relative tolerance,
    std::string expected_v_str = expected.str();
    std::string answer_v_str = answer.str();
    if ( std::abs(answer_from_lib) < VERDICT_ABSOLUTE_TOL )
    {
      std::stringstream expected_abs;
      expected_abs.setf(std::ios::scientific, std::ios::floatfield);
      expected_abs.precision(VERDICT_SIGNIFICANT_FIG);
      expected_abs << 1. + data.mAnswer;
      expected_v_str = expected_abs.str();

      std::stringstream answer_abs;
      answer_abs.setf(std::ios::scientific, std::ios::floatfield);
      answer_abs.precision(VERDICT_SIGNIFICANT_FIG);
      answer_abs << 1. + answer_from_lib;
      answer_v_str = answer_abs.str();
    }

    EXPECT_EQ(expected_v_str, answer_v_str);
    if (expected_v_str!= answer_v_str)
    {
      std::cout << std::endl;
      std::cout << "Test case \"" << testname
        << "\" #" << itest << " FAILED" << std::endl;

      std::cout << "answer calculated was    "
        << answer.str() << " ( " << answer_v_str << " ) " << std::endl
        << "answer expected was      "
        << expected.str() << " ( " << expected_v_str << " ) "
        << std::endl << std::endl;
      EXPECT_TRUE(false);
    }
    else
    {
      std::cout << "Test case \"" << testname
        << "\" #" << itest << " passed" << std::endl;
    }
    */
  }
}

TEST(verdict, tet_incircle_right)
{
  test_case testcase = { "tet_incircle_right",
    { { verdict::tet_inradius, 0.5 - 1. / std::sqrt(12) } }, 4,
    { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };

  runtest(testcase);
}

TEST(verdict, tet_incircle_right2)
{
  test_case testcase = { "tet_incircle_right2",
    { { verdict::tet_inradius, 2.0 * (0.5 - 1. / std::sqrt(12)) } }, 4,
    { { 0, 0, 0 }, { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } } };

  runtest(testcase);
}

TEST(verdict, tet_incircle_equilateral)
{
  const double pi = verdict::VERDICT_PI;

  // equilateral tet, with side length 1
  double l2 = std::sin(30 * pi / 180) / std::sin(120 * pi / 180);
  double h = std::sqrt(1 - l2 * l2);

  test_case testcase = { "tet_incircle_equilateral",
    { { verdict::tet_inradius, std::sqrt(6.) / 12. } }, 4,
    { { 0, 0, 0 }, { 1, 0, 0 }, { std::cos(60 * pi / 180), std::sin(60 * pi / 180), 0 },
      { l2 * std::cos(30 * pi / 180), l2 * std::sin(30 * pi / 180), h } } };

  runtest(testcase);
}

TEST(verdict, tet_incircle_equilateral2)
{
  const double pi = 2. * std::asin(1.);

  // equilateral tet, with side length 2
  double l2 = 2. * std::sin(30 * pi / 180) / std::sin(120 * pi / 180);
  double h = std::sqrt(4. - l2 * l2);

  test_case testcase = { "tet_incircle_equilateral2",
    { { verdict::tet_inradius, 2 * std::sqrt(6.) / 12. } }, 4,
    { { 0, 0, 0 }, { 2, 0, 0 }, { 2 * std::cos(60 * pi / 180), 2 * std::sin(60 * pi / 180), 0 },
      { l2 * std::cos(30 * pi / 180), l2 * std::sin(30 * pi / 180), h } } };

  runtest(testcase);
}

TEST(verdict, tet_incircle_flat)
{
  test_case testcase = { "tet_incircle_flat", { { verdict::tet_inradius, 0.0 } }, 4,
    { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0.5, 0.5, 0 } } };

  runtest(testcase);
}

TEST(verdict, tet_incircle_inverted)
{
  test_case testcase = { "tet_incircle_inverted",
    { { verdict::tet_inradius, -0.5 + 1. / std::sqrt(12) } }, 4,
    { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 } } };

  runtest(testcase);
}

TEST(verdict, tet_incircle_too_few_nodes)
{
  test_case testcase = { "tet_incircle_too_few_nodes", { { verdict::tet_inradius, 0. } },
    3, // not enough nodes
    {
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_incircle_regression)
{
  test_case testcase = { "tet_incircle_regression", { { verdict::tet_inradius, 3.2723844167 } }, 4,
    {
      // some arbitrary vectors, ensure verdict returns the value it used to
      { 12, 44, 13.3 },
      { 17, 0.003, 66 },
      { 3, 33, 234 },
      { 14, -123, 21 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_incircle_regression2)
{
  test_case testcase = { "tet_incircle_regression2",
    { { verdict::tet_inradius, 2 * 3.27238441675 } }, 4,
    {
      // some arbitrary vectors, ensure verdict returns the value it used to
      { 2 * 12, 2 * 44, 2 * 13.3 },
      { 2 * 17, 2 * 0.003, 2 * 66 },
      { 2 * 3, 2 * 33, 2 * 234 },
      { 2 * 14, 2 * -123, 2 * 21 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_equiangle_skew_1)
{
  test_case testcase = { "tet_equiangle_skew_1",
    { { verdict::tet_equiangle_skew, 0.25 }, { verdict::tet_equivolume_skew, 0.5 } }, 4,
    { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };

  runtest(testcase);
}

TEST(verdict, edge_calc_1)
{
  test_case testcase = { "edge_calc_1", { { verdict::edge_length, 1.732050807568877 } }, 2,
    { { 0, 0, 0 }, { 1, 1, 1 } } };

  runtest(testcase);
}

TEST(verdict, edge_calc_2)
{
  test_case testcase = { "edge_calc_2", { { verdict::edge_length, 1.0 } }, 2,
    { { 0, 0, 0 }, { 1, 0, 0 } } };

  runtest(testcase);
}

TEST(verdict, edge_calc_3)
{
  test_case testcase = { "edge_calc_3", { { verdict::edge_length, 0.0 } }, 2,
    { { 0, 0, 0 }, { 0, 0, 0 } } };

  runtest(testcase);
}

TEST(verdict, wedge_simple)
{
  test_case testcase = { "wedge_simple",
    { { verdict::wedge_volume, 0.5 }, { verdict::wedge_equiangle_skew, 0.25 },
      { verdict::wedge_edge_ratio, std::sqrt(2.) },
      { verdict::wedge_max_aspect_frobenius, 1.1357042248 },
      { verdict::wedge_mean_aspect_frobenius, 1.0978474174 }, { verdict::wedge_distortion, 1 },
      { verdict::wedge_max_stretch, 1 } },
    6, { { 0, 0, 0 }, { -1, 1, 0 }, { -1, 0, 0 }, { 0, 0, 1 }, { -1, 1, 1 }, { -1, 0, 1 } } };

  runtest(testcase);
}

TEST(verdict, wedge_singularity)
{
  test_case testcase = { "wedge_singularity",
    { { verdict::wedge_volume, 0.0 }, { verdict::wedge_equiangle_skew, 1 },
      /*  3 */ { verdict::wedge_edge_ratio, verdict::VERDICT_DBL_MAX },
      /*  4 */ { verdict::wedge_max_aspect_frobenius, verdict::VERDICT_DBL_MAX },
      /*  5 */ { verdict::wedge_mean_aspect_frobenius, verdict::VERDICT_DBL_MAX },
      /*  6 */ { verdict::wedge_distortion, verdict::VERDICT_DBL_MAX },
      /*  7 */ { verdict::wedge_max_stretch, verdict::VERDICT_DBL_MAX } },
    6, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } };

  runtest(testcase);
}

TEST(verdict, tri_simple)
{
  test_case testcase = { "tri_simple",
    { { verdict::tri_area, 10.8253175 }, { verdict::tri_aspect_ratio, 1 },
      { verdict::tri_condition, 1.0 }, { verdict::tri_distortion, 1.0 },
      { verdict::tri_minimum_angle, 60 }, { verdict::tri_maximum_angle, 60 },
      { verdict::tri_shape, 1.0 },
      /*  8 */ { verdict::tri_edge_ratio, 1 },
      /*  9 */ { verdict::tri_aspect_frobenius, 1 },
      /* 10 */ { verdict::tri_equiangle_skew, 1.8069363487e-09 },
      /* 11 */ { verdict::tri_normalized_inradius, 1.0 } },
    3, { { 0, 0, 0 }, { 5, 0, 0 }, { 2.5, 4.330127, 0 } } };

  runtest(testcase);
}

TEST(verdict, tri_singular)
{
  test_case testcase = { "tri_singular",
    { { verdict::tri_area, 0.43301270185 }, { verdict::tri_aspect_ratio, 1 },
      { verdict::tri_condition, 1 }, { verdict::tri_distortion, 1 },
      { verdict::tri_minimum_angle, 60 }, { verdict::tri_maximum_angle, 60 },
      { verdict::tri_shape, 1 },
      /*  8 */ { verdict::tri_edge_ratio, 1 },
      /*  9 */ { verdict::tri_aspect_frobenius, 1 },
      /* 10 */ { verdict::tri_equiangle_skew, 4.0316550098e-11 },
      /* 11 */ { verdict::tri_normalized_inradius, 1.0 } },
    3, { { 0, 0, 0 }, { 0.5, 0.8660254037, 0 }, { 1, 0, 0 } } };

  runtest(testcase);
}

// tri_distortion is not well covered, we need to test it with a six noded triangle */
TEST(verdict, tri_six_nodes)
{
  test_case testcase = { "tri_six_nodes",
    { { verdict::tri_area, 0.54772255751 }, { verdict::tri_aspect_ratio, 1.3268079265 },
      { verdict::tri_condition, 1.1173381066 }, { verdict::tri_distortion, -3.3198288345 },
      { verdict::tri_minimum_angle, 45.406778838 }, { verdict::tri_maximum_angle, 85.823122909 },
      { verdict::tri_shape, 0.89498424344 },
      /*  8 */ { verdict::tri_edge_ratio, 1.4005493428 },
      /*  9 */ { verdict::tri_aspect_frobenius, 1.1173381066 },
      /* 10 */ { verdict::tri_equiangle_skew, 0.24322035270 },
      /* 11 */ { verdict::tri_normalized_inradius, 0.603527 } },
    6,
    { { 0, 0, 0 }, { 1, 0, 0.2 }, { 0, 1, 0.4 }, { 0, 0.5, 0.1 }, { 0.5, 0.5, 0.3 },
      { 0.5, 0, 0.2 } } };

  runtest(testcase);
}

TEST(verdict, quad_simple1)
{
  test_case testcase = { "quad_simple1",
    {
      /* 1 */ {verdict::quad_skew, 0},
      /* 2 */ {verdict::quad_aspect_ratio, 4.0},
      /* 3 */ {verdict::quad_max_edge_ratio, 7.0}
    },
    4,
    {
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 7, 0 },
      { 0, 7, 0 }
  } };

  runtest(testcase);
}

TEST(verdict, quad_simple2)
{
  test_case testcase = { "quad_simple2",
    { /*  1 */ { verdict::quad_aspect_ratio, 1.74792 },
      /*  2 */ { verdict::quad_skew, 0.092450032704 },
      /*  3 */ { verdict::quad_taper, 0.74535599250 },
      /*  4 */ { verdict::quad_warpage, 0.008 },
      /*  5 */ { verdict::quad_area, 2.6925824036 },
      /*  6 */ { verdict::quad_stretch, 0.57735026919 },
      /*  7 */ { verdict::quad_minimum_angle, 56.789089239 },
      /*  8 */ { verdict::quad_maximum_angle, 90 },
      /*  9 */ { verdict::quad_condition, 2.3079277745 },
      /* 10 */ { verdict::quad_jacobian, 1.1141720291 },
      /* 11 */ { verdict::quad_shear, 0.55708601453 },
      /* 12 */ { verdict::quad_shape, 0.43328912241 },
      /* 13 */ { verdict::quad_distortion, 0.56267957295 },
      /* 14 */ { verdict::quad_equiangle_skew, 0.36901011957 },
      /* 15 */ { verdict::quad_radius_ratio, 1.7320508076 },
      /* 16 */ { verdict::quad_med_aspect_frobenius, 1.2274682929 },
      /* 17 */ { verdict::quad_max_aspect_frobenius, 1.3416407865 },
      /* 18 */ { verdict::quad_oddy, 1.6000000000 } },
    4,
    {
      { 2, 0, 0 }, // 1
      { 1, 1, 2 }, // 2
      { 0, 1, 0 }, // 3
      { 0, 0, 0 }, // 0
    } };

  runtest(testcase);
}

TEST(verdict, quad_chevron)
{
  test_case testcase = { "quad_chevron",
    { /*  1 */ { verdict::quad_aspect_ratio, 3.00675 },
      /*  2 */ { verdict::quad_skew, 6.5193481489e-01 },
      /*  3 */ { verdict::quad_taper, 9.8606743061e-01 },
      /*  4 */ { verdict::quad_warpage, -9.9459686559e-01 },
      /*  5 */ { verdict::quad_area, 3.5562234182 },
      /*  6 */ { verdict::quad_stretch, 7.4945733623e-01 },
      /*  7 */ { verdict::quad_minimum_angle, 19.319495622 },
      /*  8 */ { verdict::quad_maximum_angle, 226.10065328 },
      /*  9 */ { verdict::quad_condition, verdict::VERDICT_DBL_MAX },
      /* 10 */ { verdict::quad_jacobian, -3.4157021569 },
      /* 11 */ { verdict::quad_shear, 0 },
      /* 12 */ { verdict::quad_shape, 0 },
      /* 13 */ { verdict::quad_distortion, -9.0249845339e-01 },
      /* 14 */ { verdict::quad_equiangle_skew, 1.5122294809 },
      /* 15 */ { verdict::quad_radius_ratio, 2.9914983186 },
      /* 16 */ { verdict::quad_med_aspect_frobenius, 1.8922421445 },
      /* 17 */ { verdict::quad_max_aspect_frobenius, 3.5710779160 },
      /* 18 */ { verdict::quad_oddy, 23.505194964 } },
    4,
    {
      // reflex angle at node 2
      { 1, 1, 1 },     // 1
      { 3, 1, 2 },     // 2
      { 5, 0.8, 1.3 }, // 3
      { 2, 1.1, 3.7 }, // 0
    } };

  runtest(testcase);
}

TEST(verdict, quad_collapsed)
{
  test_case testcase = { "quad_collapsed",
    { /*  1 */ { verdict::quad_aspect_ratio, verdict::VERDICT_DBL_MAX },
      /*  2 */ { verdict::quad_skew, 1 },
      /*  3 */ { verdict::quad_taper, 0 },
      /*  4 */ { verdict::quad_warpage, 0 },
      /*  5 */ { verdict::quad_area, 0 },
      /*  6 */ { verdict::quad_stretch, 4.7140452079e-01 },
      /*  7 */ { verdict::quad_minimum_angle, 0 },
      /*  8 */ { verdict::quad_maximum_angle, 180 },
      /*  9 */ { verdict::quad_condition, verdict::VERDICT_DBL_MAX },
      /* 10 */ { verdict::quad_jacobian, 0 },
      /* 11 */ { verdict::quad_shear, 0 },
      /* 12 */ { verdict::quad_shape, 0 },
      /* 13 */ { verdict::quad_distortion, 0 },
      /* 14 */ { verdict::quad_equiangle_skew, 1 },
      /* 15 */ { verdict::quad_radius_ratio, verdict::VERDICT_DBL_MAX },
      /* 16 */ { verdict::quad_med_aspect_frobenius, verdict::VERDICT_DBL_MAX },
      /* 17 */ { verdict::quad_max_aspect_frobenius, verdict::VERDICT_DBL_MAX },
      /* 18 */ { verdict::quad_oddy, verdict::VERDICT_DBL_MAX } },
    4,
    {
      // colinear nodes
      { 3, 1, 1 },  // 1
      { 6, 1, 2 },  // 2
      { 12, 1, 4 }, // 3
      { 9, 1, 3 },  // 0
    } };

  runtest(testcase);
}

TEST(verdict, quad_bowtie)
{
  test_case testcase = { "quad_bowtie",
    { /*  1 */ { verdict::quad_aspect_ratio, 15.132908 },
      /*  2 */ { verdict::quad_skew, 6.0840291338e-01 },
      /*  3 */ { verdict::quad_taper, 4.9416251426e+00 },
      /*  4 */ { verdict::quad_warpage, -9.9972444493e-01 },
      /*  5 */ { verdict::quad_area, 3.4869452964e+00 },
      /*  6 */ { verdict::quad_stretch, 1.7893591011e+00 },
      /*  7 */ { verdict::quad_minimum_angle, 2.0794246576e+01 },
      /*  8 */ { verdict::quad_maximum_angle, 3.2636439841e+02 },
      /*  9 */ { verdict::quad_condition, verdict::VERDICT_DBL_MAX },
      /* 10 */ { verdict::quad_jacobian, -2.2480385649e+01 },
      /* 11 */ { verdict::quad_shear, 0 },
      /* 12 */ { verdict::quad_shape, 0 },
      /* 13 */ { verdict::quad_distortion, -2.3654968475e+00 },
      /* 14 */ { verdict::quad_equiangle_skew, 2.6262710934e+00 },
      /* 15 */ { verdict::quad_radius_ratio, 2.7490501584e+00 },
      /* 16 */ { verdict::quad_med_aspect_frobenius, 2.3523244623e+00 },
      /* 17 */ { verdict::quad_max_aspect_frobenius, 2.8418376590e+00 },
      /* 18 */ { verdict::quad_oddy, 1.4152082560e+01 } },
    4,
    {
      { -1, -1, -1 },     // 1
      { 6, 2, -1.1 },     // 2
      { 3, -2.5, -1.15 }, // 3
      { 4.5, 4.3, -0.9 }, // 0
    } };

  runtest(testcase);
}

TEST(verdict, tet_test1)
{
  test_case testcase = { "tet_test1",
    { { verdict::tet_volume, 166.66666667 }, { verdict::tet_condition, 1.2247448714 },
      { verdict::tet_jacobian, 1000 }, { verdict::tet_shape, 0.83994736660 },
      { verdict::tet_distortion, 1 },
      /*  6 */ { verdict::tet_edge_ratio, std::sqrt(2.) },
      /*  7 */ { verdict::tet_radius_ratio, 1.3660254038 },
      /*  8 */ { verdict::tet_aspect_ratio, 1.3660254038 },
      /*  9 */ { verdict::tet_aspect_frobenius, 1.1905507890 },
      /* 10 */ { verdict::tet_minimum_angle, 54.735610317 },
      /* 11 */ { verdict::tet_collapse_ratio, 0.40824829046 },
      /* 12 */ { verdict::tet_equivolume_skew, 0.5 },
      /* 13 */ { verdict::tet_normalized_inradius, 0.61401440738235424 } },
    4,
    {
      { -5, -5, -5 },
      { -5, 5, -5 },
      { -5, -5, 5 },
      { 5, -5, -5 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_test2_singular)
{
  test_case testcase = { "tet_test2_singular",
    { { verdict::tet_volume, 0 }, { verdict::tet_condition, verdict::VERDICT_DBL_MAX },
      { verdict::tet_jacobian, 0 }, { verdict::tet_shape, 0 }, { verdict::tet_distortion, 1 },
      /*  6 */ { verdict::tet_edge_ratio, verdict::VERDICT_DBL_MAX },
      /*  7 */ { verdict::tet_radius_ratio, verdict::VERDICT_DBL_MAX },
      /*  8 */ { verdict::tet_aspect_ratio, verdict::VERDICT_DBL_MAX },
      /*  9 */ { verdict::tet_aspect_frobenius, verdict::VERDICT_DBL_MAX },
      /* 10 */ { verdict::tet_minimum_angle, verdict::VERDICT_DBL_MAX },
      /* 11 */ { verdict::tet_collapse_ratio, verdict::VERDICT_DBL_MAX },
      /* 12 */ { verdict::tet_equivolume_skew, verdict::VERDICT_DBL_MAX },
      /* 13 */ { verdict::tet_normalized_inradius, verdict::VERDICT_DBL_MAX } },
    4,
    { // all points are 0, so the tet has zero edge length, volume, etc.
      { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } };

  runtest(testcase);
}

TEST(verdict, tet10_test1)
{
  test_case testcase = { "tet10_test1",
    {
      /* 0 */ { verdict::tet_distortion, 1 },
      /* 1 */ { verdict::tet_normalized_inradius, 0.61401440738235424 },
    },
    10,
    {
      { -5, 5, 5 },
      { -5, -5, -5 },
      { 5, -5, 5 },
      { -5, -5, 5 },
      { -5, 0, 0 },
      { 0, -5, 0 },
      { 0, 0, 5 },
      { -5, 0, 5 },
      { -5, -5, 0 },
      { 0, -5, 5 },
    } };

  runtest(testcase);
}

TEST(verdict, tet10_test2_singular)
{
  test_case testcase = { "tet10_test2_singular",
    {
      /* 0 */ { verdict::tet_distortion, verdict::VERDICT_DBL_MAX },
      /* 1 */ { verdict::tet_normalized_inradius, verdict::VERDICT_DBL_MAX },
    },
    10,
    {
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
    } };

  runtest(testcase);
}

TEST(verdict, hex_test1)
{
  test_case testcase = { "hex_test1_regression",
    { { verdict::hex_skew, 0.24588970374 }, { verdict::hex_taper, 0.17845765256 },
      { verdict::hex_volume, 0.81566666667 }, { verdict::hex_stretch, 0.62097029970 },
      { verdict::hex_diagonal, 0.68962192988 }, { verdict::hex_dimension, 0.52459420058 },
      { verdict::hex_condition, 1.2730598257 }, // same as verdict::hex_max_aspect_frobenius
      { verdict::hex_med_aspect_frobenius, 1.1435711972 }, { verdict::hex_jacobian, 0.477 },
      { verdict::hex_shear, 0.77784951012 }, { verdict::hex_shape, 0.78978528102 },
      { verdict::hex_distortion, 0.58479771148 },
      { verdict::hex_nodal_jacobian_ratio, 0.42513368984 }, { verdict::hex_oddy, 2.0988338297 },
      { verdict::hex_edge_ratio, 1.7944358445 }, { verdict::hex_equiangle_skew, 0.37159771083 } },
    8,
    {
      { -.2, -.7, -.3 }, // 1
      { -.7, .4, -.6 },  // 2
      { -.5, .5, .3 },   // 3
      { -.3, -.5, .5 },  // 0
      { .5, -.8, -.2 },  // 5
      { .4, .4, -.6 },   // 6
      { .2, .5, .2 },    // 7
      { .5, -.3, .8 }    // 4
    } };

  runtest(testcase);
}

TEST(verdict, hex_test2_perfect_cube)
{
  test_case testcase = { "hex_test2_perfect_cube",
    { { verdict::hex_skew, 0. }, { verdict::hex_taper, 0. }, { verdict::hex_volume, 1.0 },
      { verdict::hex_stretch, 1.0 }, { verdict::hex_diagonal, 1.0 },
      { verdict::hex_dimension, 1. / std::sqrt(3.) }, { verdict::hex_condition, 1.0 },
      { verdict::hex_med_aspect_frobenius, 1 }, { verdict::hex_jacobian, 1.0 },
      { verdict::hex_shear, 1.0 }, { verdict::hex_shape, 1.0 }, { verdict::hex_distortion, 1.0 },
      { verdict::hex_nodal_jacobian_ratio, 1.0 }, { verdict::hex_oddy, 0. },
      { verdict::hex_edge_ratio, 1 }, { verdict::hex_equiangle_skew, 0 } },
    8,
    { // perfect unit hex
      { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 },
      { 0, 1, 1 } } };

  runtest(testcase);
}

TEST(verdict, hex_test3_flat)
{
  test_case testcase = { "hex_test3_flat",
    { { verdict::hex_skew, verdict::VERDICT_DBL_MAX },
      { verdict::hex_taper, verdict::VERDICT_DBL_MAX }, { verdict::hex_volume, 0. },
      { verdict::hex_stretch, 0.19245008973 }, { verdict::hex_diagonal, 1.0 },
      { verdict::hex_dimension, 0. }, { verdict::hex_condition, verdict::VERDICT_DBL_MAX },
      { verdict::hex_med_aspect_frobenius, verdict::VERDICT_DBL_MAX }, { verdict::hex_jacobian, 0 },
      { verdict::hex_shear, 0 }, { verdict::hex_shape, 0 },
      { verdict::hex_distortion, verdict::VERDICT_DBL_MAX },
      { verdict::hex_nodal_jacobian_ratio, -verdict::VERDICT_DBL_MAX },
      { verdict::hex_oddy, verdict::VERDICT_DBL_MAX },
      { verdict::hex_edge_ratio, 1.0 / std::sqrt(0.02) }, { verdict::hex_equiangle_skew, 0.5 } },
    8,
    { // squashed flat
      { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0.1, 0.1, 0 }, { 0.9, 0.1, 0 },
      { 0.9, 0.9, 0 }, { 0.1, 0.9, 0 } } };

  runtest(testcase);
}

TEST(verdict, hex_test4_inside_out)
{
  test_case testcase = { "hex_test4_inside_out",
    { // != indicates a value that is different than the right-side-out perfect cube
      /*  1 */ { verdict::hex_skew, 0 },
      /*  2 */ { verdict::hex_taper, 0 },
      /*  3 */ { verdict::hex_volume, -1. }, // !=
      /*  4 */ { verdict::hex_stretch, 1 },
      /*  5 */ { verdict::hex_diagonal, 1.0 },
      /*  6 */ { verdict::hex_dimension, 1. / std::sqrt(3.) },
      /*  7 */ { verdict::hex_condition, verdict::VERDICT_DBL_MAX },            // !=
      /*  8 */ { verdict::hex_med_aspect_frobenius, verdict::VERDICT_DBL_MAX }, // !=
      /*  9 */ { verdict::hex_jacobian, -1. },                                  // !=
      /* 10 */ { verdict::hex_shear, 0 },                                       // !=
      /* 11 */ { verdict::hex_shape, 0 },                                       // !=
      /* 12 */ { verdict::hex_distortion, 1.0 },
      /* 13 */ { verdict::hex_nodal_jacobian_ratio, -verdict::VERDICT_DBL_MAX }, // !=
      /* 14 */ { verdict::hex_oddy, verdict::VERDICT_DBL_MAX },                  // !=
      /* 15 */ { verdict::hex_edge_ratio, 1.0 },
      /* 16 */ { verdict::hex_equiangle_skew, 0 } },
    8,
    { // perfect unit hex, but ordered inside out!
      { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 }, { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 },
      { 0, 1, 0 } } };

  runtest(testcase);
}

TEST(verdict, hex27_test1)
{
  test_case testcase = { "hex27_test1", { { verdict::hex_jacobian, 0.56021468927 } }, 27,
    { { -0.5, -0.5, 0.5 }, { -0.5, -0.5, -0.5 }, { -0.5, 0.5, -0.5 }, { -0.5, 0.5, 0.5 },
      { 0.5, -0.5, 0.5 }, { 0.5, -0.5, -0.5 }, { 0.5, 0.5, -0.5 }, { 0.5, 0.5, 0.5 },
      { -0.59093, -0.554599, 0.031478 }, { -0.551216, -0.052901, -0.59258 },
      { -0.589451, 0.534972, -0.00525 }, { -0.575503, -0.008395, 0.604191 },
      { 0.035399, -0.527937, 0.582587 }, { -0.015396, -0.537798, -0.606826 },
      { -0.027222, 0.51188, -0.576597 }, { -0.006207, 0.565371, 0.585294 },
      { 0.562112, -0.548542, -0.000139 }, { 0.570746, 0.008609, -0.595057 },
      { 0.563111, 0.550215, -0.024469 }, { 0.537894, 0.050161, 0.606614 },
      { -0.02166575, -0.0022401, 0.0023115 }, { -0.610022, -0.053727, -0.007623 },
      { 0.593828, -0.007149, -0.019913 }, { 0.019742, 0.031299, 0.64416 },
      { 0.033398, -0.032358, -0.596668 }, { 0.021178, -0.581132, 0.013062 },
      { 0.000338, 0.60863, -0.017854 } } };

  runtest(testcase);
}

TEST(verdict, hex20_test1)
{
  test_case testcase = { "hex20_test1",
    {
      { verdict::hex_jacobian, 1 },
      { verdict::hex_distortion, 3.4668482987e-01 },
    },
    20,
    { { -0.5, -0.5, 0.5 }, { -0.5, -0.5, -0.5 }, { -0.5, 0.5, -0.5 }, { -0.5, 0.5, 0.5 },
      { 0.5, -0.5, 0.5 }, { 0.5, -0.5, -0.5 }, { 0.5, 0.5, -0.5 }, { 0.5, 0.5, 0.5 },
      { -0.59093, -0.554599, 0.031478 }, { -0.551216, -0.052901, -0.59258 },
      { -0.589451, 0.534972, -0.00525 }, { -0.575503, -0.008395, 0.604191 },
      { 0.035399, -0.527937, 0.582587 }, { -0.015396, -0.537798, -0.606826 },
      { -0.027222, 0.51188, -0.576597 }, { -0.006207, 0.565371, 0.585294 },
      { 0.562112, -0.548542, -0.000139 }, { 0.570746, 0.008609, -0.595057 },
      { 0.563111, 0.550215, -0.024469 }, { 0.537894, 0.050161, 0.606614 } } };

  runtest(testcase);
}

TEST(verdict, wedge_21_test1)
{
  test_case testcase = { "wedge21_test1", { { verdict::wedge_jacobian, 0.99999775 } }, 21,
    { { 0, 0, -1 }, { 1, 0, -1 }, { 0, 1, -1 }, { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 },
      { 0.5, 0, -1 }, { 0.5, 0.5, -1 }, { 0, 0.5, -1 }, { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 },
      { 0.5, 0, 1 }, { 0.5, 0.5, 1 }, { 0, 0.5, 1 }, { .333333, .333333, 0 },
      { .333333, .333333, -1 }, { .333333, .333333, 1 }, { .5, .5, 0 }, { 0, .5, 0 },
      { 0.5, 0, 0 } } };

  runtest(testcase);
}

TEST(verdict, tet15_test1)
{
  test_case testcase = { "tet15_test1", { { verdict::tet_jacobian, 2.0 } }, 15,
    { { 0, 0, 0 }, { 0, 2, 0 }, { 1, 1, 0 }, { 6.12323399573677e-17, 1, -1 }, { 0, 1, 0 },
      { 0.5, 1.5, 0 }, { 0.5, 0.5, 0 }, { 3.06161699786838e-17, 0.5, -0.5 },
      { 3.06161699786838e-17, 1.5, -0.5 }, { 0.707106781186548, 1, -0.707106781186547 },
      { 0.301776695296637, 1, -0.301776695296637 }, { 0.333333333333333, 1, 0 },
      { 0.425380791638466, 1.33333333333333, -0.425380791638466 },
      { 0.425380791638466, 0.666666666666667, -0.425380791638466 },
      { 2.04107799857892e-17, 1, -0.333333333333333 } } };

  runtest(testcase);
}

TEST(verdict, knife_test1)
{
  test_case testcase = { "knife_test1", { { verdict::knife_volume, 2. / 3. } }, 7,
    { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0.5, 0.5, 1 },
      { 1, 1, 1 } } };

  runtest(testcase);
}

TEST(verdict, tet_meanRatio_perfect_tet)
{
  test_case testcase = { "tet_meanRatio_perfect_tet",
    {
      { verdict::tet_mean_ratio, 1.0 },
    },
    4,
    {
      { 0, 0, 0 },
      { 0, 1, 1 },
      { 1, 0, 1 },
      { 1, 1, 0 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_meanRatio_degenerate_tet)
{
  test_case testcase = { "tet_meanRatio_degenerate_tet",
    {
      { verdict::tet_mean_ratio, 0.0 },
    },
    4,
    {
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_meanRatio_unit_right_angle_tet)
{
  test_case testcase = { "tet_meanRatio_unit_right_angle_tet",
    {
      { verdict::tet_mean_ratio, 0.769800 },
    },
    4,
    {
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_meanRatio_nearly_flat_right_angle_tet)
{
  test_case testcase = { "tet_meanRatio_nearly_flat_right_angle_tet",
    {
      { verdict::tet_mean_ratio, 0.01414107 },
    },
    4,
    {
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, .01 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_meanRatio_inverted_nearly_flat_right_angle_tet)
{
  test_case testcase = { "tet_meanRatio_nearly_flat_right_angle_tet",
    {
      { verdict::tet_mean_ratio, -0.014141 },
    },
    4,
    {
      { 0, 0, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, -.01 },
    } };

  runtest(testcase);
}

TEST(verdict, tet10_meanRatio_perfect_tet)
{
  test_case testcase = {
    "tet10_meanRatio_perfect_tet",
    {
      { verdict::tet_mean_ratio, 1.0},
    },
    10,
    {
      { 0, 0, 0 },
      { 0, 1, 1 },
      { 1, 0, 1 },
      { 1, 1, 0 },
      { 0, 0.5, 0.5 },
      { 0.5, 0.5, 1 },
      { 0.5, 0, 0.5 },
      { 0.5, 0.5, 0 },
      { 0.5, 1, 0.5 },
      { 1, 0.5, 0.5 },      
    }
  };

  runtest(testcase);
}

TEST(verdict, tet10_meanRatio_imperfect_tet)
{
  test_case testcase = {
    "tet10_meanRatio_imperfect_tet",
    {
      { verdict::tet_mean_ratio, 0.39802657202178215 },
    },
    10,
    {
      {-0.05, -0.1, 0.125},
      {0.9, 0.025, -0.075},
      {-0.1, 1.025, -0.075},
      {0.125, 0.1, 1.025},
      {0.525, 0.025, -0.05},
      {0.5, 0.525, 0.05},
      {0.125, 0.45, -0.075},
      {0.1, 0.075, 0.575},
      {0.475, -0.075, 0.55},
      {-0.075, 0.6, 0.45},
    }
  };

  runtest(testcase);
}


TEST(verdict, tet_normalized_inradius_perfect_tet10)
{
  test_case testcase = { "tet_normalized_subtet_perfect_tet",
    {
      { verdict::tet_normalized_inradius, 1.0 },
    },
    10,
    {
      { 0.000000, 0.000000, 0.000000 },
      { 0.000000, 1.000000, 1.000000 },
      { 1.000000, 0.000000, 1.000000 },
      { 1.000000, 1.000000, 0.000000 },
      { 0.000000, 0.500000, 0.500000 },
      { 0.500000, 0.500000, 1.000000 },
      { 0.500000, 0.000000, 0.500000 },
      { 0.500000, 0.500000, 0.000000 },
      { 0.500000, 1.000000, 0.500000 },
      { 1.000000, 0.500000, 0.500000 },
    } };

  runtest(testcase);
}

TEST(verdict, tet_normalized_inradius_deformed_tet10)
{
  test_case testcase = { "tet_normalized_subtet_perfect_tet",
    {
      { verdict::tet_normalized_inradius, 0.49373899958502704 },
    },
    10,
    {
      { -0.05, -0.1, 0.125 },
      { 0.9, 0.025, -0.075 },
      { -0.1, 1.025, -0.075 },
      { 0.125, 0.1, 1.025 },
      { 0.525, 0.025, -0.05 },
      { 0.5, 0.525, 0.05 },
      { 0.125, 0.45, -0.075 },
      { 0.1, 0.075, 0.575 },
      { 0.475, -0.075, 0.55 },
      { -0.075, 0.6, 0.45 },
    } };

  runtest(testcase);
}

TEST(verdict, tet8_volume)
{
  test_case testcase = { "Volume of tet8 element",
    {
      { verdict::tet_volume, 0.2103399378 },
    },
    8,
    {
      { -0.05, -0.1, 0.125 },
      { 0.9, 0.025, -0.075 },
      { -0.1, 1.025, -0.075 },
      { 0.125, 0.1, 1.025 },
      { 0.224518, 0.314238, 0.088176 },
      { 0.450716, 0.495582, 0.292855 },
      { 0.026687, 0.383860, 0.337325 },
      { 0.362298, -0.064834, 0.403524 },
    } };

  runtest(testcase);
}

TEST(verdict, tet10_volume)
{
  test_case testcase = { "Volume of tet10 element",
    {
      { verdict::tet_volume, 0.153515625 },
    },
    10,
    {
      { -0.05, -0.1, 0.125 },
      { 0.9, 0.025, -0.075 },
      { -0.1, 1.025, -0.075 },
      { 0.125, 0.1, 1.025 },
      { 0.525, 0.025, -0.05 },
      { 0.5, 0.525, 0.05 },
      { 0.125, 0.45, -0.075 },
      { 0.1, 0.075, 0.575 },
      { 0.475, -0.075, 0.55 },
      { -0.075, 0.6, 0.45 },
    } };

  runtest(testcase);
}

TEST(verdict, tet14_volume)
{
  test_case testcase = { "Volume of tet14 element",
    {
      { verdict::tet_volume, 0.1680788275 },
    },
    14,
    {
      { -0.05, -0.1, 0.125 },
      { 0.9, 0.025, -0.075 },
      { -0.1, 1.025, -0.075 },
      { 0.125, 0.1, 1.025 },
      { 0.525, 0.025, -0.05 },
      { 0.5, 0.525, 0.05 },
      { 0.125, 0.45, -0.075 },
      { 0.1, 0.075, 0.575 },
      { 0.475, -0.075, 0.55 },
      { -0.075, 0.6, 0.45 },
      { 0.392383, 0.4289156667, -0.007145333333 },
      { 0.3433533333, 0.4255263333, 0.2706586667 },
      { 0.02896466667, 0.2684996667, 0.4035243333 },
      { 0.299518, 0.005904333333, 0.4548423333 },
    } };

  runtest(testcase);
}

TEST(verdict, tet15_volume)
{
  test_case testcase = { "Volume of tet15 element",
    {
      { verdict::tet_volume, 0.1680788275 },
    },
    15,
    {
      { -0.05, -0.1, 0.125 },
      { 0.9, 0.025, -0.075 },
      { -0.1, 1.025, -0.075 },
      { 0.125, 0.1, 1.025 },
      { 0.525, 0.025, -0.05 },
      { 0.5, 0.525, 0.05 },
      { 0.125, 0.45, -0.075 },
      { 0.1, 0.075, 0.575 },
      { 0.475, -0.075, 0.55 },
      { -0.075, 0.6, 0.45 },
      { 0.303125, 0.26875, 0.25 },
      { 0.392383, 0.4289156667, -0.007145333333 },
      { 0.3433533333, 0.4255263333, 0.2706586667 },
      { 0.02896466667, 0.2684996667, 0.4035243333 },
      { 0.299518, 0.005904333333, 0.4548423333 },
    } };

  runtest(testcase);
}

TEST(verdict, hex20_volume)
{
  test_case testcase = { "Volume of hex20 element",
    {
      { verdict::hex_volume, 0.9425 },
    },
    20,
    {
      { -0.5, -0.5, 0.5 },
      { -0.5, -0.5, -0.5 },
      { -0.5, 0.5, -0.5 },
      { -0.5, 0.5, 0.4 },
      { 0.5, -0.5, 0.5 },
      { 0.5, -0.5, -0.5 },
      { 0.6, 0.7, -0.5 },
      { 0.5, 0.5, 0.5 },
      { -0.5, -0.5, 0.0 },
      { -0.5, 0.0, -0.5 },
      { -0.5, 0.5, 0.0 },
      { -0.5, 0.0, 0.5 },
      { 0.3, -0.5, 0.5 },
      { 0.0, -0.5, -0.5 },
      { 0.0, 0.5, -0.5 },
      { 0.0, 0.4, 0.5 },
      { 0.3, -0.5, 0.0 },
      { 0.6, 0.0, -0.5 },
      { 0.5, 0.5, 0.0 },
      { 0.5, 0.0, 0.4 },
    } };

  runtest(testcase);
}

TEST(verdict, hex27_volume)
{
  test_case testcase = { "Volume of hex27 element",
    {
      { verdict::hex_volume, 1.242167 },
    },
    27,
    { { -0.5, -0.5, 0.5 }, { -0.5, -0.5, -0.5 }, { -0.5, 0.5, -0.6 }, { -0.5, 0.5, 0.4 },
      { 0.5, -0.5, 0.5 }, { 0.5, -0.5, -0.5 }, { 0.6, 0.7, -0.5 }, { 0.5, 0.5, 0.5 },
      { -0.5, -0.5, 0.0 }, { -0.5, 0.0, -0.5 }, { -0.5, 0.5, 0.0 }, { -0.5, 0.0, 0.5 },
      { 0.3, -0.5, 0.5 }, { 0.0, -0.5, -0.5 }, { 0.0, 0.5, -0.5 }, { 0.0, 0.5, 0.8 },
      { 0.8, -0.5, 0.0 }, { 0.5, 0.0, -0.5 }, { 0.7, 0.5, 0.0 }, { 0.5, 0.0, 0.6 },
      { 0.0, 0.0, 0.0 }, { -0.5, 0.0, 0.0 }, { 0.5, 0.0, 0.0 }, { 0.0, 0.0, 0.5 },
      { 0.0, 0.0, -0.5 }, { 0.0, -0.5, 0.0 }, { 0.3, 0.8, 0.0 } } };

  runtest(testcase);
}
