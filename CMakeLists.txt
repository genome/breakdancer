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

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/build-common/cmake
)

include(TestHelper)
include(VersionHelper)
include(ProjectHelper)
include(BuildDeps)
include(PlatformChecks)

set(EXE_VERSION_SUFFIX ${FULL_VERSION})

configure_file(version/version.h.in version/version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/version)


###########################################################################
# Setup compiler flags
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU|Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic -fPIC")
endif()

if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -fPIC ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor -Wno-deprecated-register")
else()
    message(
        WARNING
        "Unknown compiler: ${CMAKE_C_COMPILER_ID}, we'll see how this goes..."
    )
endif()

# Find C++11 flags
find_cxx11_flags(CXX11_FLAGS CXX11_FOUND)
if(NOT CXX11_FOUND)
    message(FATAL_ERROR
        "Failed to find C++11 support with compiler:
         ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}).")
endif(NOT CXX11_FOUND)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_FLAGS}")

###########################################################################
# Build dependencies (samtools and boost)
add_custom_target(deps ALL)

set(SAMTOOLS_URL ${CMAKE_SOURCE_DIR}/vendor/samtools-0.1.19.tar.gz)
set(BOOST_URL ${CMAKE_SOURCE_DIR}/vendor/boost-1.54-breakdancer.tar.gz)
set(REQUIRED_BOOST_LIBS filesystem program_options serialization system chrono regex)

build_boost(${BOOST_URL} ${CMAKE_BINARY_DIR}/vendor/boost ${REQUIRED_BOOST_LIBS})
include_directories(${Boost_INCLUDE_DIRS})

build_samtools(${SAMTOOLS_URL} ${CMAKE_BINARY_DIR}/vendor/samtools)
include_directories(${Samtools_INCLUDE_DIRS})

# make sure to pick up headers from library dirs
include_directories("src/lib")

###########################################################################
## Main project
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
add_projects(src/lib)
add_projects(src/exe)
add_subdirectory(perl)

###########################################################################
## Testing
set(TEST_DATA_DIRECTORY ${CMAKE_SOURCE_DIR}/test-data)

configure_file(test-data/TestData.hpp.in test-data/TestData.hpp @ONLY)
include_directories(${PROJECT_BINARY_DIR}/test-data)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/test-bin)
enable_testing(true)
add_projects(test/lib)

# integration testing
add_subdirectory(integration-test)


###########################################################################
## Packaging
include(PackageInfo.cmake)
