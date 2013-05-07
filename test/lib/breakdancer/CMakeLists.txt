cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS breakdancer ${Boost_LIBRARIES})
include_directories(${GTEST_INCLUDE_DIRS})

#def_test(Region)
def_test(BamReader)
def_test(BreakDancer)
def_test(Config)
def_test(Poisson)
def_test(Read)
def_test(Utility)
