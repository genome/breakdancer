cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS config ${Boost_LIBRARIES})
include_directories(${GTEST_INCLUDE_DIRS})

add_unit_tests(TestConfigLib
    TestBamConfig.cpp
    TestBamConfigEntry.cpp
)
