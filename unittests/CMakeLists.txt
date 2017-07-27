PROJECT(unittests_verdict)

INCLUDE_DIRECTORIES( SYSTEM
                     ${gtest_SOURCE_DIR}
                     ${gtest_SOURCE_DIR}/include)

SET(TEST_SRCS
    unittest_main.cpp
    verdict.test.cpp
   )

ADD_EXECUTABLE(unittests_verdict ${TEST_SRCS})
TARGET_LINK_LIBRARIES(unittests_verdict verdict gtest gtest_main)
