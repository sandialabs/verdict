
SET(TEST_SRCS
    unittest_main.cpp
    verdict.test.cpp
   )

ADD_EXECUTABLE(unittests_verdict ${TEST_SRCS})
TARGET_LINK_LIBRARIES(unittests_verdict verdict GTest::GTest GTest::Main)
