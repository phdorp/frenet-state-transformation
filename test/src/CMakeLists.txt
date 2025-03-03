add_executable(
    MathTest
    MathTest.cpp
)

target_include_directories(MathTest PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(
    MathTest
    PRIVATE
    transform
    GTest::gtest_main
)

gtest_discover_tests(MathTest)

add_executable(
    PolylineTest
    PolylineTest.cpp
)

target_include_directories(PolylineTest PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(
    PolylineTest
    PRIVATE
    transform
    GTest::gtest_main
)

gtest_discover_tests(PolylineTest)

add_executable(
    CircleTest
    CircleTest.cpp
)

target_include_directories(CircleTest PRIVATE "${CMAKE_SOURCE_DIR}/test/include")

target_link_libraries(
    CircleTest
    PRIVATE
    transform
    GTest::gtest_main
)

gtest_discover_tests(CircleTest)