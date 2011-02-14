project(breakdancer)

set(SOURCES
    AlnParser.cpp
    AlnParser.h
    Poisson.cpp
    Poisson.h
    samtools.h
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer bam z m)
