project(breakdancer)

set(SOURCES
    AlignmentFilter.hpp
    BamIO.hpp
    BamMerger.cpp
    BamMerger.hpp
    BamReader.cpp
    BamReader.hpp
    FastqWriter.cpp
    FastqWriter.hpp
    IBamReader.hpp
    Read.cpp
    Read.hpp
    RegionLimitedBamReader.hpp
)

add_library(io ${SOURCES})
target_link_libraries(io common ${Samtools_LIBRARIES} z m)
