#include "frenetTransform/point.h"
#include "frenetTransform/points.h"
#include "frenetTransform/polychain.h"
#include "frenetTransform/transform.h"

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace FrenetTransform {
using namespace pybind11::literals;

PYBIND11_MODULE(_core, handle) {
  py::classh<Point>(handle, "Point")
      .def(py::init<>(), R"doc(Construct a default Point at (0, 0).)doc")
      .def(py::init<double, double>(), "x"_a, "y"_a,
           R"doc(Construct a new Point with coordinates in x- and y-direction.

Args:
    x (float): point's coordinate in x-direction.
    y (float): point's coordinate in y-direction.
Returns:
    Point: Point object.
)doc")
      .def("x", &Point::x,
           R"doc(Provides the point's coordinate in x-direction.

Returns:
    float: x-coordinate.
)doc")
      .def("y", &Point::y,
           R"doc(Provides the point's coordinate in y-direction.

Returns:
    float: y-coordinate.
)doc")
      .def("distanceSquare", &Point::distanceSquare, "point"_a,
           R"doc(Squared distance between Point and a query point.

Args:
    point (Point): Point to determine the squared distance to.
Returns:
    float: squared distance to point.
)doc")
      .def("distance", &Point::distance, "point"_a,
           R"doc(Distance between Point and a query point.

Args:
    point (Point): Point to determine the distance to.
Returns:
    float: distance to point.
)doc")
      .def("__neg__", &Point::operator-,
           R"doc(Negate x- and y-coordinates.

Returns:
    Point: Point with negated x- and y-coordinates.
)doc")
      .def(py::self + py::self,
           R"doc(Sum between coordinates of two points.

Returns:
    Point: Point with sum of coordinates of two points.
)doc")
      .def(py::self - py::self,
           R"doc(Difference between coordinates of two points.

Returns:
    Point: Point with difference of coordinates of two points.
)doc");

  using PointsD = Points<Eigen::Dynamic>;
  py::classh<PointsD>(handle, "Points")
      .def(
          py::init<Eigen::ArrayXd, Eigen::ArrayXd>(), "x"_a, "y"_a,
          R"doc(Construct a new Points object from sequences of x- and y-coordinates. The number of coordinates in both directions must be equal.

Args:
    x (np.ndarray): coordinates in x-direction of point sequence.
    y (np.ndarray): coordinates in y-direction of point sequence.
Returns:
    Points: Points object.
)doc")
      .def("__call__", &PointsD::operator(), "index"_a,
           R"doc(Get point at index along the point sequence.

Args:
    index (int): index of point along the point sequence.
Returns:
    Point: point at the index.
)doc")
      .def("numPoints", &PointsD::numPoints,
           R"doc(Provide the number of points along the sequence.

Returns:
    int: number of points.
)doc")
      .def("x", py::overload_cast<const int>(&PointsD::x, py::const_),
           "index"_a,
           R"doc(Get x-coordinate at index along the point sequence.

Args:
    index (int): index of x-coordinate along the point sequence.
Returns:
    float: x-coordinate.
)doc")
      .def("x", py::overload_cast<>(&PointsD::x, py::const_),
           R"doc(Get x-coordinates along the point sequence.

Returns:
    np.ndarray: x-coordinates.
)doc")
      .def("y", py::overload_cast<const int>(&PointsD::y, py::const_),
           "index"_a,
           R"doc(Get y-coordinate at index along the point sequence.

Args:
    index (int): index of y-coordinate along the point sequence.
Returns:
    float: y-coordinate.
)doc")
      .def("y", py::overload_cast<>(&PointsD::y, py::const_),
           R"doc(Get y-coordinates along the point sequence.

Returns:
    np.ndarray: y-coordinates.
)doc")
      .def("setX", &PointsD::setX, "x"_a,
           R"doc(Set x-coordinates along the point sequence.

Args:
    x (np.ndarray): new x-coordinates.
)doc")
      .def("setY", &PointsD::setY, "y"_a,
           R"doc(Set y-coordinates along the point sequence.

Args:
    y (np.ndarray): new y-coordinates.
)doc")
      .def("distanceSquare", &PointsD::distanceSquare, "point"_a,
           R"doc(Squared distance between Points and point.

Args:
    point (Point): point to determine the squared distances to.
Returns:
    np.ndarray: squared distances to point.
)doc")
      .def("distance", &PointsD::distance, "point"_a,
           R"doc(Distance between Points and point.

Args:
    point (Point): point to determine the distances to.
Returns:
    np.ndarray: distances to point.
)doc")
      .def("__neg__", &PointsD::operator-,
           R"doc(Negate x- and y-coordinates.

Returns:
    Points: Points with negated x- and y-coordinates.
)doc")
      .def(py::self + py::self,
           R"doc(Element-wise sum of two Points.

Returns:
    Points: Points with sum of coordinates of two Points.
)doc")
      .def(py::self - py::self,
           R"doc(Element-wise difference between two Points.

Returns:
    Points: Points with difference between coordinates of two Points.
)doc")
      .def(py::self * py::self,
           R"doc(Element-wise scalar product of two Points.

Args:
    other (Points): other Points sequence for scalar product.
Returns:
    np.ndarray: scalar products of Points.
)doc");

  using PathD = Path<Eigen::Dynamic>;
  py::classh<PathD> path(handle, "Path");
  path.def("__call__", &PathD::operator(), "lengths"_a,
           R"doc(Determines points at the given path lengths.

Args:
    lengths (np.ndarray): lengths along the path.
Returns:
    Points: points at given path lengths.
)doc")
      .def("lengths", &PathD::lengths, "points"_a,
           R"doc(Determines next points to the query points.

Args:
    points (Points): query points.
Returns:
    np.ndarray: next to query points.
)doc")
      .def("tangent", &PathD::tangent, "lengths"_a,
           R"doc(Determines tangent vectors at the given path lengths.

Args:
    lengths (np.ndarray): lengths along the path.
Returns:
    Points: tangent vectors.
)doc")
      .def("normal", &PathD::normal, "lengths"_a,
           R"doc(Determines normal vectors at the given path lengths.

Args:
    lengths (np.ndarray): lengths along the path.
Returns:
    Points: normal vectors.
)doc")
      .def("angle0", &PathD::angle0, "lengths"_a,
           R"doc(Determines path angle at the given path lengths.

Args:
    lengths (np.ndarray): lengths along the path.
Returns:
    np.ndarray: path angles.
)doc")
      .def("angle1", &PathD::angle1, "lengths"_a,
           R"doc(Determines path curvature at the given path lengths.

Args:
    lengths (np.ndarray): lengths along the path.
Returns:
    np.ndarray: path curvatures.
)doc")
      .def("angle2", &PathD::angle2, "lengths"_a,
           R"doc(Determines path curvature derivative at the given path lengths.

Args:
    lengths (np.ndarray): lengths along the path.
Returns:
    np.ndarray: path curvatures.
)doc");

  using PolychainD = Polychain<Eigen::Dynamic>;
  py::classh<PolychainD>(handle, "Polychain", path)
      .def(
          py::init<Eigen::ArrayXd, Eigen::ArrayXd>(), "x"_a, "y"_a,
          R"doc(Construct a new Polychain object from Cartesian x- and y-positions.

Args:
    x (np.ndarray): coordinates in x-direction along the path.
    y (np.ndarray): coordinates in y-direction along the path.
Returns:
    Polychain: Polychain object.
)doc")
      .def(
          "setPoints", &PolychainD::setPoints, "x"_a, "y"_a,
          R"doc(Provide new points for the polychain. Update lengths and gradient information.

Args:
    x (np.ndarray): coordinates in x-direction of new points.
    y (np.ndarray): coordinates in y-direction of new points.
)doc")
      .def(
          "__call__", &PolychainD::operator(), "lengths"_a,
          R"doc(Gets points along the polychain at the query lengths. Lengths exceeding the polychain's domain either resolve to the first or last point.

Args:
    lengths (np.ndarray): query lengths along the polychain.
Returns:
    Points: points at the query lengths.
)doc")
      .def(
          "lengths", &PolychainD::lengths, "points"_a,
          R"doc(Determines next points to the query points. Performs a linear search over all polychain segments to identify the closest one.

Args:
    points (Points): query points.
Returns:
    np.ndarray: next points to query points.
)doc");

  using TransformD = Transform<Eigen::Dynamic>;
  py::classh<TransformD>(handle, "Transform")
      .def(py::init<std::shared_ptr<PathD>>(), "path"_a,
           R"doc(Construct a new Transform object from a given Path.

Args:
    path (Path): defines the Frenet frame.
Returns:
    Transform: Transform object.
)doc")
      .def(
          "posFrenet", &TransformD::posFrenet, "posCartes"_a,
          R"doc(Transform Cartesian positions to Frenet positions. Projects the query points onto the path. Determines the signed lengths along the path from the path origin to the projections. Determines the signed shortest distances to the query point.

Args:
    posCartes (Points): query points in Cartesian coordinates.
Returns:
    Points: result points in Frenet coordinates.
)doc")
      .def("posCartes", &TransformD::posCartes, "posFrenet"_a,
           R"doc(Transform Frenet positions to Cartesian positions.

Args:
    posFrenet (Points): query points in Frenet coordinates.
Returns:
    Points: result points in Cartesian coordinates.
)doc")
      .def("velFrenet", &TransformD::velFrenet, "velCartes"_a, "posFrenet"_a,
           R"doc(Transform Cartesian velocities to Frenet velocities.

Args:
    velCartes (Points): query velocities in Cartesian coordinates.
    posFrenet (Points): positions corresponding to velocities in Frenet coordinates.
Returns:
    Points: result velocities in Frenet coordinates.
)doc")
      .def("velCartes", &TransformD::velCartes, "velFrenet"_a, "posFrenet"_a,
           R"doc(Transform Frenet velocities to Cartesian velocities.

Args:
    velFrenet (Points): query velocities in Frenet coordinates.
    posFrenet (Points): positions corresponding to velocities in Frenet coordinates.
Returns:
    Points: result velocities in Cartesian coordinates.
)doc")
      .def("accFrenet", &TransformD::accFrenet, "accCartes"_a, "velFrenet"_a,
           "posFrenet"_a,
           R"doc(Transform Cartesian accelerations to Frenet accelerations.

Args:
    accCartes (Points): query accelerations in Cartesian coordinates.
    velFrenet (Points): velocities corresponding to accelerations in Frenet coordinates.
    posFrenet (Points): positions corresponding to accelerations in Frenet coordinates.
Returns:
    Points: result accelerations in Frenet coordinates.
)doc")
      .def("accCartes", &TransformD::accCartes, "accFrenet"_a, "velFrenet"_a,
           "posFrenet"_a,
           R"doc(Transform Frenet accelerations to Cartesian accelerations.

Args:
    accFrenet (Points): query accelerations in Frenet coordinates.
    velFrenet (Points): velocities corresponding to accelerations in Frenet coordinates.
    posFrenet (Points): positions corresponding to accelerations in Frenet coordinates.
Returns:
    Points: result accelerations in Cartesian coordinates.
)doc");
}
} // namespace FrenetTransform