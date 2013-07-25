cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS breakdancer ${Boost_LIBRARIES})
include_directories(${GTEST_INCLUDE_DIRS})

add_unit_tests(TestBdLib
    TestBreakDancer.cpp
    TestReadCountsByLib.cpp
)
