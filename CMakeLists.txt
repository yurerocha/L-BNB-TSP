# To build on Linux/MacOS:
#   cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release [Debug]
#   cmake --build build

cmake_minimum_required(VERSION 3.08)

# Set the project name and version.
project(bnb VERSION 0.1 LANGUAGES CXX)

option(CXX "enable C++ compilation" ON)

# file(MAKE_DIRECTORY sol)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Specify the C++ standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler flags.
#   -Wall: Recommended compiler warnings.
#   -Wextra: Enables some extra warning flags that are not enabled by -Wall.
#   -Wpedantic: Issue all the warnings demanded by strict ISO C and ISO C++.
#   -Werror: Make all warnings into errors.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wno-deprecated-copy -pedantic-errors")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

# list all source files here
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(SRC_FILES
  ${SRC_DIR}/main.cpp
  ${SRC_DIR}/Lagrange.cpp
  ${SRC_DIR}/BNB.cpp
  ${SRC_DIR}/ext/data.cpp
  ${SRC_DIR}/ext/hungarian.cpp
  ${SRC_DIR}/ext/Kruskal.cpp
  ${SRC_DIR}/utils/Utils.cpp
)

# add the executable
include_directories(include include/ext include/utils)
add_executable(${CMAKE_PROJECT_NAME} ${SRC_FILES})