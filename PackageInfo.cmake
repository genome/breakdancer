cmake_minimum_required(VERSION 2.8)

# .deb packaging
set(ARCH "i686")
if(${CMAKE_SIZEOF_VOID_P} MATCHES 8)
    set(ARCH "x86_64")
endif ()

# The format of the description field is a short summary line followed by a
# longer paragraph indented by a single space on each line
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SV detection from paired end read mapping in Illumina Data
 BreakDancer performs genome-wide detection of structural variants from next generation paired-end sequencing reads. breakdancer-max predicts five types of structural variants: insertions, deletions, inversions, inter- and intra-chromosomal translocations from Illumina next-generation short paired-end sequencing reads by using read pairs that are mapped with unexpected separation distances or orientation. Please read our paper for detailed algorithmic description. http://www.nature.com/nmeth/journal/v6/n9/abs/nmeth.1363.html
 ")

set(CPACK_PACKAGE_NAME "breakdancer${EXE_VERSION_SUFFIX}")
set(CPACK_PACKAGE_VENDOR "The Genome Institute")
set(CPACK_PACKAGE_VERSION ${FULL_VERSION}${PACKAGE_VERSION_SUFFIX})
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "David Larson <dlarson@genome.wustl.edu> and Travis Abbott <tabbott@genome.wustl.edu>")
set(CPACK_SYSTEM_NAME "Linux-${ARCH}")
set(CPACK_TOPLEVEL_TAG "Linux-${ARCH}")
set(CPACK_DEBIAN_PACKAGE_SECTION science)
set(CPACK_DEBIAN_PACKAGE_PRIORITY optional)
set(CPACK_DEBIAN_PACKAGE_REPLACES "breakdancer1.0")
set(CPACK_DEBIAN_PACKAGE_CONFLICTS "breakdancer1.0")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.4), libgcc1 (>= 1:4.1.1-21), libstdc++6 (>= 4.2.1-4), libstatistics-descriptive-perl (>=2.9-1), libgd-graph-histogram-perl (>=1.1-1), zlib1g")
if (CMAKE_BUILD_TYPE MATCHES package)
    set(CPACK_GENERATOR "DEB")
else(CMAKE_BUILD_TYPE MATCHES package)
    set(CPACK_GENERATOR "TGZ")
endif(CMAKE_BUILD_TYPE MATCHES package)

configure_file(debian/postinst.in debian/postinst @ONLY)
configure_file(debian/prerm.in debian/prerm @ONLY)
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "debian/postinst;debian/prerm")

include(CPack)

