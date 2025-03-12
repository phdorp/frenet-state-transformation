cmake_minimum_required(VERSION 3.22)

project(FrenetTransform)

include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

add_subdirectory(src)

option(BUILD_TEST "Build tests for transformation." ON)
option(BUILD_BENCHMARK "Build benchmarks for transformations." ON)

if(BUILD_TEST)
    enable_testing()
    add_subdirectory(test/src/tests)
endif()

if(BUILD_BENCHMARK)
    add_subdirectory(test/src/benchmarks)
endif()