project(breakdancer)

set(SOURCES
    BDConfig.cpp
    BDConfig.hpp
    BamConfig.cpp
    BamConfig.hpp
    BamMerger.cpp
    BamMerger.hpp
    BamReader.cpp
    BamReader.hpp
    BasicRegion.hpp
    BreakDancer.hpp
    IBamReader.hpp
    Options.hpp
    Read.cpp
    Read.hpp
    ReadCountsByLib.hpp
    ReadFlags.hpp
    Region.cpp
    Region.hpp
    namespace.hpp
    utility.hpp
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer ${Samtools_LIBRARIES} z m)
