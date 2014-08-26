project(breakdancer-max)

set(SOURCES
    Params.cpp
    Params.hpp
    bam2cfg.cpp
)

set(EXECUTABLE_NAME bam2cfg)
add_executable(${EXECUTABLE_NAME} ${SOURCES})
target_link_libraries(${EXECUTABLE_NAME} common io breakdancer ${Boost_LIBRARIES})
set_target_properties(${EXECUTABLE_NAME} PROPERTIES PACKAGE_OUTPUT_NAME ${EXECUTABLE_NAME}${EXE_VERSION_SUFFIX})
install(TARGETS ${EXECUTABLE_NAME} DESTINATION bin/)
