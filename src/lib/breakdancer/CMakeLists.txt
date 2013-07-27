project(breakdancer)

set(SOURCES
    BasicRegion.hpp
    BedWriter.cpp
    BedWriter.hpp
    BreakDancer.cpp
    BreakDancer.hpp
    ReadCountsByLib.hpp
    ReadRegionData.cpp
    ReadRegionData.hpp
    SvBuilder.cpp
    SvBuilder.hpp
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer common io config ${Samtools_LIBRARIES} z m)
