include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.16.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

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

include(GoogleTest)
gtest_discover_tests(MathTest)
