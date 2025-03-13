# frenet-state-transformation
Transformation of dynamic state from Cartisian to Frenet coordinates and vice versa

```
mkdir build && cd build
cmake .. -DBUILD_TEST=ON -DBUILD_BENCHMARK=ON
cmake --build .
```

## Usage with CMake

```
FetchContent_Declare(
  transform
  GIT_REPOSITORY https://github.com/CoteGoal/frenet-state-transformation
)
set(BUILD_TEST OFF)
set(BUILD_BENCHMARK OFF)
FetchContent_MakeAvailable(transform)

...

target_link_libraries("target" PRIVATE transform)
```