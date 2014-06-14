project(breakdancer)

set(SOURCES
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
    IReadClassifier.hpp
    IlluminaPEReadClassifier.cpp
    IlluminaPEReadClassifier.hpp
    LibraryConfig.hpp
    LibraryFlagDistribution.cpp
    LibraryFlagDistribution.hpp
    LibraryInfo.hpp
    RawBamEntry.hpp
    Read.cpp
    Read.hpp
    RegionLimitedBamReader.hpp
)

add_library(io ${SOURCES})
target_link_libraries(io common ${Samtools_LIBRARIES} z m)
