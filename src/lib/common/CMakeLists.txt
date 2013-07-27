project(breakdancer)

set(SOURCES
    ConfigMap.hpp
    Graph.hpp
    Options.cpp
    Options.hpp
    ReadFlags.cpp
    ReadFlags.hpp
    Timer.hpp
    namespace.hpp
    utility.hpp
)

add_library(common ${SOURCES})
