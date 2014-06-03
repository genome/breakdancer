cmake_minimum_required(VERSION 2.8)

include(ExternalProject)

function(build_boost BOOST_URL BUILD_DIR)
    set(REQUIRED_BOOST_LIBS ${ARGN})

    set(BOOST_ROOT ${BUILD_DIR}/boost)
    set(BOOST_SRC ${BUILD_DIR}/boost-src)
    set(BOOST_BUILD_LOG ${BOOST_SRC}/build.log)

    foreach(libname ${REQUIRED_BOOST_LIBS})
        set(BOOST_BUILD_LIBS ${BOOST_BUILD_LIBS} --with-${libname})
        set(BOOST_LIBS ${BOOST_LIBS}
            ${BOOST_ROOT}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}boost_${libname}${CMAKE_STATIC_LIBRARY_SUFFIX}
            )
    endforeach(libname ${REQUIRED_BOOST_LIBS})
    set(Boost_LIBRARIES ${BOOST_LIBS} PARENT_SCOPE)

    ExternalProject_Add(
        boost-libs
        URL ${BOOST_URL}
        SOURCE_DIR ${BOOST_SRC}
        BINARY_DIR ${BOOST_SRC}
        CONFIGURE_COMMAND "./bootstrap.sh"
        BUILD_COMMAND
            ./b2 --prefix=${BOOST_ROOT} --layout=system link=static threading=multi install
                ${BOOST_BUILD_LIBS} ${BOOST_BUILD_OPTS} > ${BOOST_BUILD_LOG}
        INSTALL_COMMAND ""
    )

    add_dependencies(deps boost-libs)
    set(Boost_INCLUDE_DIRS ${BOOST_ROOT}/include PARENT_SCOPE)

endfunction(build_boost BOOST_VERSION BUILD_DIR)

function(build_samtools SAMTOOLS_URL BUILD_DIR)
    set(SAMTOOLS_ROOT ${BUILD_DIR})
    set(SAMTOOLS_LOG ${SAMTOOLS_ROOT}/build.log)
    set(SAMTOOLS_LIB ${SAMTOOLS_ROOT}/${CMAKE_FIND_LIBRARY_PREFIXES}bam${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(SAMTOOLS_BIN ${SAMTOOLS_ROOT}/samtools)

    ExternalProject_Add(
        samtools-lib
        URL ${SAMTOOLS_URL}
        SOURCE_DIR ${SAMTOOLS_ROOT}
        BINARY_DIR ${SAMTOOLS_ROOT}
        CONFIGURE_COMMAND echo "Building samtools, build log at ${SAMTOOLS_LOG}"
        BUILD_COMMAND make > ${SAMTOOLS_LOG} 2>&1
        INSTALL_COMMAND ""
    )

    find_package(ZLIB)
    if (NOT ZLIB_FOUND)
        message(FATAL_ERROR
            "Unable to find zlib development headers/library."
            " On Debian-based distributions, you can install these with"
            " sudo apt-get install zlib1g-dev")
    endif (NOT ZLIB_FOUND)

    add_library(bam STATIC IMPORTED)
    set_property(TARGET bam PROPERTY IMPORTED_LOCATION ${SAMTOOLS_LIB})

    set(Samtools_INCLUDE_DIRS ${ZLIB_INCLUDE_DIRS};${SAMTOOLS_ROOT} PARENT_SCOPE)
    set(Samtools_LIBRARIES bam m ${ZLIB_LIBRARIES} PARENT_SCOPE)

    add_dependencies(deps samtools-lib)
endfunction(build_samtools SAMTOOLS_URL BUILD_DIR)
