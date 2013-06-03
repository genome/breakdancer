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
    FlagCodec.hpp
    IBamReader.hpp
    IlluminaLIFlagCodec.hpp
    IlluminaStandardFlagCodec.hpp
    LibraryInfo.hpp
    Options.cpp
    Options.hpp
    Read.cpp
    Read.hpp
    ReadCountsByLib.hpp
    ReadFlags.cpp
    ReadFlags.hpp
    RegionLimitedBamReader.hpp
    namespace.hpp
    utility.hpp
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer ${Samtools_LIBRARIES} z m)
