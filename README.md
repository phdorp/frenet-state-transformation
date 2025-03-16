# Frenet transformation

This project implements a library to perform transformations from a 2-dimensional Cartesian to a Frenet coordinate system and vice versa.
Transformations are implemented for states up to the second time derivative.
The reference path defining the Frenet coordinate system is represented by a polychain.

## Installation

Since this is a header-only library it is sufficient to include the required header files from the directory *include/Transform*.

### Integration in CMake project

If you desire an integration in a CMake project you may utilize the *FetchContent* module.

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

### Build project

Building the project with the benchmark and tests requires enabling the respective CMake configuration flags.

```
mkdir build
cmake -Bbuild -DBUILD_TEST=ON -DBUILD_BENCHMARK=ON
cmake --build build
```

## Usage



## Introduction Frenet coordinates