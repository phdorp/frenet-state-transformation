find_package(benchmark REQUIRED)

add_executable(circleBenchmark circleBenchmark.cpp)

target_include_directories(circleBenchmark PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(circleBenchmark PRIVATE
    transform
    benchmark::benchmark)