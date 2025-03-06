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
)

gtest_discover_tests(pathCircleTest)