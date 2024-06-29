cmake_minimum_required(VERSION 3.14)
project(WCAPSolver VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

option(ENABLE_VERBOSE_LOG "Enable --verbose logging flag" ON)

if(ENABLE_VERBOSE_LOG)
  add_compile_definitions(DBG)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(ARGTABLE3_ENABLE_EXAMPLES OFF)
set(ARGTABLE3_ENABLE_TESTS OFF)

include(FetchContent)
FetchContent_Declare(
  pugixml
  GIT_REPOSITORY https://github.com/zeux/pugixml
  GIT_TAG master
)
FetchContent_Declare(
  argtable
  GIT_REPOSITORY https://github.com/argtable/argtable3
  GIT_TAG master
)
FetchContent_MakeAvailable(pugixml)
FetchContent_MakeAvailable(argtable)

include_directories(${argtable_SOURCE_DIR}/src)
include_directories(${pugixml_SOURCE_DIR}/src)

find_package(Boost REQUIRED)
find_package(LEMON REQUIRED)

include_directories(
  ./src
)

add_library(util src/util.cpp)
add_library(graph src/graph.cpp)
add_library(greedy src/greedy.cpp)
add_library(ilp src/ilp.cpp)
add_library(config src/config.cpp)
add_library(approx src/approx.cpp)
add_library(ring_graph src/ring_graph.cpp)
add_library(watanabe src/watanabe.cpp)
target_link_libraries(config PRIVATE argtable3 util)
target_link_libraries(graph util)
target_link_libraries(approx ring_graph util)

add_executable(Solver src/Solver.cpp)
target_link_libraries(Solver PUBLIC graph ilp greedy approx pugixml config watanabe util)

# Gurobi
find_path(GUROBI_INCLUDE_DIRS
  NAMES gurobi_c.h
  HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
  PATH_SUFFIXES include)

find_library(GUROBI_LIBRARY
  NAMES gurobi gurobi100
  HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
  PATH_SUFFIXES lib)

find_library(GUROBI_CXX_LIBRARY
  NAMES gurobi_c++
  HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
  PATH_SUFFIXES lib)
set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_CXX_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY)

include_directories(${GUROBI_INCLUDE_DIRS})

target_link_libraries(Solver PUBLIC optimized ${GUROBI_CXX_LIBRARY}
  debug ${GUROBI_CXX_DEBUG_LIBRARY})
target_link_libraries(Solver PUBLIC ${GUROBI_LIBRARY})
