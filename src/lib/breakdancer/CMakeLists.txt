project(breakdancer)

set(SOURCES
    AlignmentFilter.hpp
    BDConfig.cpp
    BDConfig.hpp
    BamConfig.cpp
    BamConfig.hpp
    BamIO.hpp
    BamMerger.cpp
    BamMerger.hpp
    BamReader.cpp
    BamReader.hpp
    BasicRegion.hpp
    BreakDancer.cpp
    BreakDancer.hpp
    IBamReader.hpp
    LibraryInfo.hpp
    Options.cpp
    Options.hpp
    Read.cpp
    Read.hpp
    ReadCountsByLib.hpp
    ReadFlags.hpp
    RegionLimitedBamReader.hpp
    namespace.hpp
    utility.hpp
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer ${Samtools_LIBRARIES} z m)
