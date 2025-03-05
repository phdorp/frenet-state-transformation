add_library(transform path.cpp polyline.cpp transform.cpp)

target_include_directories(transform PUBLIC "${CMAKE_SOURCE_DIR}/include")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)