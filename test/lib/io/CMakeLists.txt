cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS io ${Boost_LIBRARIES})
include_directories(${GTEST_INCLUDE_DIRS})

add_unit_tests(TestIoLib
    TestBam.cpp
    TestBamConfig.cpp
    TestBamConfigEntry.cpp
    TestBamIo.cpp
    TestBamMerger.cpp
    TestBamReader.cpp
    TestIlluminaPEReadClassifier.cpp
    TestLibraryFlagDistribution.cpp
    TestRead.cpp
    TestRegionLimitedBamReader.cpp
)
