cmake_minimum_required(VERSION 2.8)

include(ExternalProject)

set(SAMTOOLS_ROOT ${CMAKE_BINARY_DIR}/vendor/samtools)

ExternalProject_Add(
    samtools
    URL ${CMAKE_SOURCE_DIR}/vendor/samtools-0.1.19.tar.gz
    SOURCE_DIR ${SAMTOOLS_ROOT}
    BINARY_DIR ${SAMTOOLS_ROOT}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)

set(Samtools_INCLUDE_DIRS ${SAMTOOLS_ROOT})
set(Samtools_LIBRARIES ${SAMTOOLS_ROOT}/${CMAKE_FIND_LIBRARY_PREFIXES}bam${CMAKE_STATIC_LIBRARY_SUFFIX} m z)
