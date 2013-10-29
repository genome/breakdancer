cmake_minimum_required(VERSION 2.8)

include(ExternalProject)

set(SAMTOOLS_ROOT ${CMAKE_BINARY_DIR}/vendor/samtools)
set(SAMTOOLS_LIB ${SAMTOOLS_ROOT}/${CMAKE_FIND_LIBRARY_PREFIXES}bam${CMAKE_STATIC_LIBRARY_SUFFIX})
set(SAMTOOLS_BIN ${SAMTOOLS_ROOT}/samtools)

ExternalProject_Add(
    samtools-0.1.19
    URL ${CMAKE_SOURCE_DIR}/vendor/samtools-0.1.19.tar.gz
    SOURCE_DIR ${SAMTOOLS_ROOT}
    BINARY_DIR ${SAMTOOLS_ROOT}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND cp ${SAMTOOLS_ROOT}/samtools ${CMAKE_BINARY_DIR}/bin
)

add_library(bam STATIC IMPORTED)
set_property(TARGET bam PROPERTY IMPORTED_LOCATION ${SAMTOOLS_LIB})

set(Samtools_INCLUDE_DIRS ${SAMTOOLS_ROOT})
set(Samtools_LIBRARIES bam m z)
