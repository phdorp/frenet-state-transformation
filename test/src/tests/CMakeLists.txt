include(GoogleTest)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.16.0.zip
)
FetchContent_MakeAvailable(googletest)

add_executable(
    mathTest
    mathTest.cpp
)

target_include_directories(mathTest PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(
    mathTest
    PRIVATE
    transform
    GTest::gtest_main
    eigen
)

gtest_discover_tests(mathTest)

add_executable(
    pathPolylineTest
    pathPolylineTest.cpp
)

target_include_directories(pathPolylineTest PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(
    pathPolylineTest
    PRIVATE
    transform
    GTest::gtest_main
    eigen
)

gtest_discover_tests(pathPolylineTest)

add_executable(
    pathCircleTest
    pathCircleTest.cpp
)

target_include_directories(pathCircleTest PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(
    pathCircleTest
    PRIVATE
    transform
    GTest::gtest_main
    eigen
)

gtest_discover_tests(pathCircleTest)

add_executable(
    pathPolylineLineTest
    pathPolylineLineTest.cpp
)

target_include_directories(pathPolylineLineTest PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(
    pathPolylineLineTest
    PRIVATE
    transform
    GTest::gtest_main
    eigen
)

gtest_discover_tests(pathPolylineLineTest)

