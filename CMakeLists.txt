cmake_minimum_required(VERSION 3.9.0)
set(PROJECT_NAME textopt)

project(${PROJECT_NAME})

set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "1")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

enable_language(CXX)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 17)

set(NO_OPENMP OFF CACHE BOOL "Compile without OpenMP support")
set(BUILD_NATIVE ON CACHE BOOL "Compile using the architecture currently in use")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
    if(BUILD_NATIVE)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Og")
endif()

find_package(SFML COMPONENTS system window graphics REQUIRED)
find_package(YAML-CPP REQUIRED)

if(NOT NO_OPENMP)
    find_package(OpenMP)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
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
    ${YAML_CPP_INCLUDE_DIR}
)
   
add_subdirectory(src)
