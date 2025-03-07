cmake_minimum_required(VERSION 3.22)

project(FrenetTransform)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

add_subdirectory(src)

option(BUILD_TESTING "..." OFF)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test/src/tests)
    add_subdirectory(test/src/benchmarks)
endif()