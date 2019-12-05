#include "gtest/gtest.h"

void InitializeTest( int argc, char **argv );
void CloseTest( void );

int main(int argc, char **argv)
{
  try
  {
    ::testing::InitGoogleTest(&argc, argv);
    int return_status = RUN_ALL_TESTS();
    return return_status;
  }
  catch(std::exception& e)
  {
    printf("Exception Caught \"%s\": Test Failed\n", e.what());
    return 1;
  }
  catch(...)
  {
    printf("Unexpected Exception Caught: Test Failed\n");
    return 1;
  }
}

