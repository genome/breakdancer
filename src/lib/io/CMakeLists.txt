project(breakdancer)

set(SOURCES
    Alignment.cpp
    Alignment.hpp
    AlignmentFilter.hpp
    BamConfig.cpp
    BamConfig.hpp
    BamConfigEntry.cpp
    BamConfigEntry.hpp
    BamIo.cpp
    BamIo.hpp
    BamMerger.cpp
    BamMerger.hpp
    BamReader.hpp
    BamReaderBase.hpp
    BamSummary.cpp
    BamSummary.hpp
    BamWriter.cpp
    BamWriter.hpp
    ConfigLoader.cpp
    ConfigLoader.hpp
    FastqWriter.cpp
    FastqWriter.hpp
    IAlignmentClassifier.hpp
    IlluminaPEReadClassifier.cpp
    IlluminaPEReadClassifier.hpp
    LibraryConfig.hpp
    LibraryFlagDistribution.cpp
    LibraryFlagDistribution.hpp
    LibraryInfo.hpp
    RawBamEntry.hpp
    RegionLimitedBamReader.hpp
)

add_library(io ${SOURCES})
target_link_libraries(io common ${Samtools_LIBRARIES} z m)
