cmake_minimum_required(VERSION 2.8)

project(breakdancer)

if(NOT BD_SCORE_FLOAT_TYPE)
    set(BD_SCORE_FLOAT_TYPE double)
endif()
message("Score floating point type: ${BD_SCORE_FLOAT_TYPE}")
add_definitions("-DSCORE_FLOAT_TYPE=${BD_SCORE_FLOAT_TYPE}")

if(NOT CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "/usr")
endif()
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/build-common/cmake")
include(TestHelper)
include(VersionHelper)
include(ProjectHelper)

set(EXE_VERSION_SUFFIX ${FULL_VERSION})

configure_file(version/version.h.in version/version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/version)

set(SAMTOOLS_ADDITIONAL_HEADERS sam.h)
find_package(Samtools 0.1.6 REQUIRED)
include_directories(${Samtools_INCLUDE_DIRS})

find_package(Boost 1.48 REQUIRED filesystem program_options system chrono regex)
link_libraries(rt)
include_directories(${Boost_INCLUDE_DIRS})

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-std=c99 -Wall -ansi -pedantic")
endif ()

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wno-long-long -Wall -ansi -pedantic -Wno-variadic-macros")
endif ()

# make sure to pick up headers from library dirs
include_directories("src/lib")

# libgtest requires threads, so link thread libraries everywhere
find_package(Threads)
link_libraries(${CMAKE_THREAD_LIBS_INIT})

enable_testing(true)
add_projects(test/lib)

# integration testing
add_subdirectory(integration-test)


# main project

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
add_projects(src/lib)
add_projects(src/exe)
add_subdirectory(perl)

include(PackageInfo.cmake)
