FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        c58e6d0710581e3a08d65c349664128a8d9a2461 # release-1.9.1
)
set(BENCHMARK_ENABLE_TESTING OFF)
FetchContent_MakeAvailable(googlebenchmark)

add_executable(circleBenchmark circleBenchmark.cpp)

target_include_directories(circleBenchmark PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(circleBenchmark PRIVATE
    transform
    benchmark::benchmark)