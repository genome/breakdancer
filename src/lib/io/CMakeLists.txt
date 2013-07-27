project(breakdancer)

set(SOURCES
    AlignmentFilter.hpp
    BamIo.cpp
    BamIo.hpp
    BamMerger.cpp
    BamMerger.hpp
    BamReader.hpp
    BamReaderBase.hpp
    FastqWriter.cpp
    FastqWriter.hpp
    RawBamEntry.hpp
    Read.cpp
    Read.hpp
    RegionLimitedBamReader.hpp
)

add_library(io ${SOURCES})
target_link_libraries(io common ${Samtools_LIBRARIES} z m)
