cmake_minimum_required(VERSION 3.5.0)
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

if(WIN32 OR CYGWIN)
    list(APPEND GMSH_LIB winmm wsock32 ws2_32 psapi)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
endif()


include_directories(
    ${PROJECT_SOURCE_DIR}/include
    SYSTEM
    ${SFML_INCLUDE_DIR}
)
   
add_subdirectory(src)
