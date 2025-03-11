FetchContent_Declare(
  eigen
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4/eigen-3.4.zip
)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
FetchContent_MakeAvailable(eigen)

add_library(transform path.cpp polyline.cpp transform.cpp)
target_link_libraries(transform PRIVATE eigen)

target_include_directories(transform PUBLIC "${CMAKE_SOURCE_DIR}/include")