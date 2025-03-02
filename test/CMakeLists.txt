include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.16.0.zip
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_executable(
    MathTest
    MathTest.cpp
)
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
target_link_libraries(
  PolylineTest
  PRIVATE
  transform
  GTest::gtest_main
)
gtest_discover_tests(PolylineTest)

