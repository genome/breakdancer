project(breakdancer)

set(SOURCES
    AlnParser.cpp
    AlnParser.h
    BDTypedefs.hpp
    LegacyBamReader.cpp
    LegacyBamReader.hpp
    Region.cpp
    Region.hpp
    samtools.h
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer ${Samtools_LIBRARIES} z m)
