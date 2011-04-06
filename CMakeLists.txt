cmake_minimum_required(VERSION 2.8)

project(breakdancer)

set(CMAKE_INSTALL_PREFIX "/usr")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/build-common/cmake")
include(TestHelper)
include(VersionHelper)
include(ProjectHelper)

configure_file(version/version.h.in version/version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/version)

# TODO: replace with FindSamtools.cmake
set(SAMTOOLS_INC /gsc/pkg/bio/samtools/samtools-0.1.6)
set(SAMTOOLS_LIB /gsc/pkg/bio/samtools/samtools-0.1.6)
include_directories(${SAMTOOLS_INC})
link_directories(${SAMTOOLS_LIB})

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-std=c99 -Wall -ansi -pedantic")
endif ()

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wno-long-long -Wall -ansi -pedantic")
endif ()

find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# make sure to pick up headers from library dirs
include_directories("src/lib")

# unit tests
find_package(GTest)
if(GTEST_FOUND)
    message("Google Test framework found, building unit tests")
    enable_testing(true)
    add_projects(test/lib)
else()
    message("Google Test framework not found, no tests will be built")
    message("GTest is available at http://code.google.com/p/googletest/")
    message("Ubuntu users can likely sudo apt-get install gtest")
endif()

# main project

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
add_projects(src/lib)
add_projects(src/exe)
add_subdirectory(perl)

include(PackageInfo.cmake)
