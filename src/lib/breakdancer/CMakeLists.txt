project(breakdancer)

set(SOURCES
    AlnParser.cpp
    AlnParser.h
    samtools.h
    Region.hpp
    Region.cpp
    BDTypedefs.hpp
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer ${Samtools_LIBRARIES} z m)
