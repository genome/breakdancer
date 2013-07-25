project(breakdancer)

set(SOURCES
    BamConfig.cpp
    BamConfig.hpp
    BamSummary.cpp
    BamSummary.hpp
    LibraryConfig.hpp
    LibraryFlagDistribution.hpp
    LibraryInfo.hpp
)

add_library(config ${SOURCES})
target_link_libraries(config common io)
