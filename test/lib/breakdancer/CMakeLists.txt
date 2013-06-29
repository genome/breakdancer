cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS breakdancer ${Boost_LIBRARIES})
include_directories(${GTEST_INCLUDE_DIRS})

add_unit_tests(TestBdLib
    TestBamReader.cpp
    TestBreakDancer.cpp
    TestConfig.cpp
    TestPoisson.cpp
    TestRead.cpp
    TestReadCountsByLib.cpp
    TestUtility.cpp
)
