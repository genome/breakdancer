project(breakdancer)

set(SOURCES
    AlignmentFilter.hpp
    BamIo.cpp
    BamIo.hpp
    BamMerger.cpp
    BamMerger.hpp
    BamReader.hpp
    BamReaderBase.hpp
    BamWriter.cpp
    BamWriter.hpp
    FastqWriter.cpp
    FastqWriter.hpp
    IReadClassifier.hpp
    IlluminaPEReadClassifier.cpp
    IlluminaPEReadClassifier.hpp
    RawBamEntry.hpp
    Read.cpp
    Read.hpp
    RegionLimitedBamReader.hpp
)

add_library(io ${SOURCES})
target_link_libraries(io common ${Samtools_LIBRARIES} z m)
