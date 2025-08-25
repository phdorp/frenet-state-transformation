# Frenet transformation

![Frenet transform of position, velocity, and acceleration](docs/media/circlePolychain.png)

This project implements a library to perform transformations from a 2-dimensional Cartesian to a Frenet coordinate system and vice versa.
Transformations are implemented for states up to the second time derivative.
The reference path defining the Frenet coordinate system is represented by a polychain.

## Installation

Since this is a header-only library it is sufficient to include the required header files from the directory *include/Transform*.

### Integration in CMake project

If you desire an integration in a CMake project you may utilize the *FetchContent* module.

```bash
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

```bash
cmake -Bbuild
cmake --build build
```

The following options are avaible for configuration:
- BUILD_DOCS=ON|OFF: whether to build the [documentation](build/docs/html/index.html)
- BUILD_TEST=ON|OFF: whether to build the tests
- BUILD_BENCHMARK=ON|OFF: whether to build the benchmarks
- BUILD_EXAMPLES=ON|OFF: whether to build the examples

### Python bindings

```bash
pip install -e .
```

```bash
python bindings/python/frenetTransform/build_stubs.py $PATH_BUILD
```

## Examples

For running the examples setup the environment variables and invoke the script *examples/examples.bash* from the repository root.

```bash
source setup.bash
examples/examples.bash [options] <source-file-names>
```

Invoke `examples/examples.bash -h` for additional information.
Note that the examples require [gnuplot](http://gnuplot.info/) to plot the results.

## Benchmarks

For running the benchmarks setup the environment variables and invoke the script *benchmarks/benchmarks.bash* from the repository root.

```bash
source setup.bash
benchmarks/benchmarks.bash [options] <source-file-names>
```

Invoke `benchmarks/benchmarks.bash -h` for additional information.
