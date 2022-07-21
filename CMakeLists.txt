cmake_minimum_required(VERSION 3.9.0)
set(PROJECT_NAME textopt)

project(${PROJECT_NAME})

set(VERSION_MAJOR "0")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

enable_language(CXX)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(SFML COMPONENTS system window graphics REQUIRED)
find_package(OpenMP)

if(WIN32 OR CYGWIN)
    list(APPEND GMSH_LIB winmm wsock32 ws2_32 psapi)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -DNDEBUG")
endif()

# Enable link time optimization if supported
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        list(APPEND CMAKE_CXX_COMPILE_OPTIONS_IPO "-flto=auto")
    endif()
else()
    message(WARNING "IPO/LTO is not supported: ${output}")
endif()


include_directories(
    ${PROJECT_SOURCE_DIR}/include
    SYSTEM
    ${SFML_INCLUDE_DIR}
)
   
add_subdirectory(src)
