/*=========================================================================

  Module:    V_EdgeMetric.cpp

  Copyright 2003,2006,2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
  Under the terms of Contract DE-NA0003525 with NTESS,
  the U.S. Government retains certain rights in this software.

  See LICENSE for details.

=========================================================================*/

#include "gtest/gtest.h"

void InitializeTest(int argc, char** argv);
void CloseTest(void);

int main(int argc, char** argv)
{
  try
  {
    ::testing::InitGoogleTest(&argc, argv);
    int return_status = RUN_ALL_TESTS();
    return return_status;
  }
  catch (std::exception& e)
  {
    printf("Exception Caught \"%s\": Test Failed\n", e.what());
    return 1;
  }
  catch (...)
  {
    printf("Unexpected Exception Caught: Test Failed\n");
    return 1;
  }
}
