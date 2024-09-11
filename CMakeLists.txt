cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(CtrModel)

# general compilation flags
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall /std:c++17")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17")
endif()

# Boost library (required by gnuplot-iostream, used for plotting)
FIND_PACKAGE( Boost 1.40 COMPONENTS iostreams system filesystem REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
LINK_LIBRARIES( ${Boost_LIBRARIES} )

# OpenMP library (used for parallel computing)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(SYSTEM Library)
include_directories(SYSTEM Library/eigen-3.4.0)
include_directories(SYSTEM src)

add_subdirectory(src)
add_subdirectory(demo)

set(CMAKE_DISABLE_SOURCE_CHANGES OFF)
set(CMAKE_DISABLE_IN_SOURCE_BUILD OFF)
