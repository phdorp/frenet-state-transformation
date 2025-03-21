# Frenet transformation

![Frenet transformation](docs/img/stateTransform.png)

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

The following example includes snippets from the file [examples/circlePolychain.cpp](examples/circlePolychain.cpp)
First, create a polychain path as a shared pointer.
In this example, the polychain approximates a circle.
Coordinates are represented with Eigen arrays.
The example uses dynamic arrays, but static arrays are also supported.
Create a separate transform object with the polychain.

```
// create circle with radius 10 m
const double radius { 10.0 };
const Eigen::ArrayXd lengthsCircle { Eigen::ArrayXd::LinSpaced(101, 0.0, 2 * M_PI) };
const Eigen::ArrayXd circlePointsX { radius * lengthsCircle.cos() };
const Eigen::ArrayXd circlePointsY { radius * lengthsCircle.sin() };

// generate polychain along circle
const auto polychain { std::make_shared<FrenetTransform::Polychain<Eigen::Dynamic>>(circlePointsX, circlePointsY) };

// instantiate transform
FrenetTransform::Transform<Eigen::Dynamic> transform { polychain };
```

For the transformation from the Cartesian to the Frenet coordinate frame define positions, velocities, and accelerations in the Cartesian frame.
In this case, a set of positions is arranged on a grid.
The velocities and accelerations are the same for each point.

```
// create point grid from -15 to 15 in x- and y-direction
const double bound { 15.0 };
auto [posMeshX, posMeshY] { matplot::meshgrid(matplot::iota(0.5 - bound, 1, 0.5 + bound), matplot::iota(0.5 - bound, 1, 0.5 + bound)) };

// copy mesh to Cartesian points
FrenetTransform::Points<Eigen::Dynamic> cartesPoints { toArray(ravel(posMeshX)), toArray(ravel(posMeshY))};

// instantiate Cartesian velocities and accelerations
const FrenetTransform::Points<Eigen::Dynamic> cartesVels { Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 2, Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 2 };
const FrenetTransform::Points<Eigen::Dynamic> cartesAccs { 3 * Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 4, -3 * Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 4 };
```

Finally, query the transform object to perform the transformation to the Frenet coordinates.
High-order derivatives require the solution from the low-order ones.

```
// transform query points to Frenet frame
const FrenetTransform::Points<Eigen::Dynamic> frenetPointsTf { transform.posFrenet(cartesPoints) };
const FrenetTransform::Points<Eigen::Dynamic> frenetVelsTf { transform.velFrenet(cartesVels, frenetPointsTf) };
const FrenetTransform::Points<Eigen::Dynamic> frenetAccsTf { transform.accFrenet(cartesAccs, frenetVelsTf, frenetPointsTf) };
```

For the transformation back to the Cartesian coordinate system query the transformation object with points the points in the Frenet frame.
Note that the tranformation from Cartesian to the Frenet frame is in general surjective when performed with a polyline i.e. there are multiple Cartesian points that map to one Frenet point if the next point on the polychain is a vertex point.
Thus the last back transformation to the Cartesian frame cannot reproduce the input.

```
// transform query points back to Cartesian frame
const FrenetTransform::Points<Eigen::Dynamic> cartesPointsTf { transform.posCartes(frenetPointsTf) };
const FrenetTransform::Points<Eigen::Dynamic> cartesVelsTf { transform.velCartes(frenetVelsTf, frenetPointsTf) };
const FrenetTransform::Points<Eigen::Dynamic> cartesAccsTf { transform.accCartes(frenetAccsTf, frenetVelsTf, frenetPointsTf) };
```