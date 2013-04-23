project(breakdancer)

set(SOURCES
    AlnParser.cpp
    AlnParser.h
    BDConfig.cpp
    BDConfig.hpp
    BDTypedefs.hpp
    BamMerger.cpp
    BamMerger.hpp
    BamReader.cpp
    BamReader.hpp
    IBamReader.hpp
    LegacyBamReader.cpp
    LegacyBamReader.hpp
    Read.cpp
    Read.hpp
    Region.cpp
    Region.hpp
    samtools.h
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer ${Samtools_LIBRARIES} z m)
