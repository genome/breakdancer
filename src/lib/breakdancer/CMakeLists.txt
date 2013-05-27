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
    BreakDancer.cpp
    BreakDancer.hpp
    IBamReader.hpp
    LibraryInfo.hpp
    Options.hpp
    Options.cpp
    Read.cpp
    Read.hpp
    ReadCountsByLib.hpp
    ReadFlags.hpp
    namespace.hpp
    utility.hpp
)

add_library(breakdancer ${SOURCES})
target_link_libraries(breakdancer ${Samtools_LIBRARIES} z m)
