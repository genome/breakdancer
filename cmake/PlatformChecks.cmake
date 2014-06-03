cmake_minimum_required(VERSION 2.8)

include(CheckIncludeFileCXX)
include(CheckCXXCompilerFlag)
include(CheckCXXSourceRuns)
include(CheckFunctionExists)

function(find_cxx11_flags FLAGS FOUND)
    set(CXX11_FOUND False BOOL)
    list(APPEND CANDIDATE_FLAGS "-std=c++11" "-std=c++0x" "-std=gnu++11" "-std=gnu++0x")

    set(COUNTER 0)
    set(__CXX11_FLAG_FOUND False BOOL)
    unset(__CXX11_FLAG)
    foreach(elt ${CANDIDATE_FLAGS})
        unset(CXX11_FLAG${COUNTER} CACHE)
        check_cxx_compiler_flag("${elt}" CXX11_FLAG${COUNTER})
        if (CXX11_FLAG${COUNTER})
            unset(CXX11_FLAG${COUNTER} CACHE)
            set(__CXX11_FLAG ${elt})
            set(__CXX11_FLAG_FOUND True BOOL)
            message("-- C++11 support enabled via ${elt}")
            break()
        endif()
        math(EXPR COUNTER "${COUNTER} + 1")
    endforeach()

    if(NOT __CXX11_FLAG_FOUND)
        message("Failed to find C++11 compiler flag, perhaps one is not needed.")
    endif(NOT __CXX11_FLAG_FOUND)

    set(SAVE_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    list(APPEND STDLIB_CANDIDATE_FLAGS " " "-stdlib=libc++" "-stdlib=libstdc++")

    set(COUNTER 0)
    foreach(elt ${STDLIB_CANDIDATE_FLAGS})
        unset(CXX11_STDLIB_FLAG${COUNTER} CACHE)
        set(CMAKE_CXX_FLAGS "${SAVE_CMAKE_CXX_FLAGS} ${__CXX11_FLAG} ${elt}")
        check_cxx_source_runs("
            #include <vector>
            #include <utility>
            #include <memory>

            bool check_move(std::vector<int>&& v) {
                return v.size() == 1 && v[0] == 3;
            }

            int main(int argc, char** argv) {
                std::unique_ptr<int> x(new int(5));
                std::vector<int> v;
                v.emplace_back(3);
                auto ok = check_move(std::move(v));
                return (*x == 5 && ok) ? 0 : 1;
            }
            " CXX11_STDLIB_FLAG${COUNTER})
        if (CXX11_STDLIB_FLAG${COUNTER})
            unset(CXX11_STDLIB_FLAG${COUNTER} CACHE)
            set(${FLAGS} "${__CXX11_FLAG} ${elt}" PARENT_SCOPE)
            set(${FOUND} True BOOL PARENT_SCOPE)
            message("-- Sufficient C++11 library support found with flag '${elt}'")
            break()
        endif (CXX11_STDLIB_FLAG${COUNTER})

        math(EXPR COUNTER "${COUNTER} + 1")
    endforeach()

    set(CMAKE_CXX_FLAGS ${SAVE_CMAKE_CXX_FLAGS})
endfunction(find_cxx11_flags FLAGS FOUND)

function(find_library_providing func FOUND LIB)
    set(POTENTIAL_LIBRARIES "" ${ARGN})

    set(ORIG_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    foreach(lib ${POTENTIAL_LIBRARIES})
        set(CMAKE_REQUIRED_LIBRARIES ${ORIG_CMAKE_REQUIRED_LIBRARIES} ${lib})
        check_function_exists(${func} ${FOUND})
        if (${FOUND})
            message("-- * clock_gettime found in library ${lib}")
            set(${LIB} ${lib} PARENT_SCOPE)
            break()
        endif (${FOUND})
    endforeach(lib ${POTENTIAL_LIBRARIES})
    set(CMAKE_REQUIRED_LIBRARIES ${ORIG_CMAKE_REQUIRED_LIBRARIES})
endfunction(find_library_providing func FOUND LIB)
