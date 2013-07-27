project(breakdancer)

set(SOURCES
    BamConfig.cpp
    BamConfig.hpp
    BamConfigEntry.cpp
    BamConfigEntry.hpp
    BamSummary.cpp
    BamSummary.hpp
    ConfigLoader.cpp
    ConfigLoader.hpp
    LibraryConfig.hpp
    LibraryFlagDistribution.cpp
    LibraryFlagDistribution.hpp
    LibraryInfo.hpp
)

add_library(config ${SOURCES})
target_link_libraries(config common io)
