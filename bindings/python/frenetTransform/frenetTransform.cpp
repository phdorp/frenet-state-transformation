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
      .def(py::init<>())
      .def(py::init<double, double>(), "x"_a, "y"_a)
      .def("x", &Point::x)
      .def("y", &Point::y)
      .def("distanceSquare", &Point::distanceSquare, "point"_a)
      .def("distance", &Point::distance, "point"_a)
      .def("__neg__", &Point::operator-)
      .def(py::self + py::self)
      .def(py::self - py::self);

  using PointsD = Points<Eigen::Dynamic>;
  py::classh<PointsD>(handle, "Points")
      .def(py::init<Eigen::ArrayXd, Eigen::ArrayXd>(), "x"_a, "y"_a)
      .def("__call__", &PointsD::operator(), "index"_a)
      .def("numPoints", &PointsD::numPoints)
      .def("x", py::overload_cast<const int>(&PointsD::x, py::const_),
           "index"_a)
      .def("x", py::overload_cast<>(&PointsD::x, py::const_))
      .def("y", py::overload_cast<const int>(&PointsD::y, py::const_),
           "index"_a)
      .def("y", py::overload_cast<>(&PointsD::y, py::const_))
      .def("setX", &PointsD::setX, "x"_a)
      .def("setY", &PointsD::setY, "y"_a)
      .def("distanceSquare", &PointsD::distanceSquare, "point"_a)
      .def("distance", &PointsD::distance, "point"_a)
      .def("__neg__", &PointsD::operator-)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self * py::self);

  using PathD = Path<Eigen::Dynamic>;
  py::classh<PathD> path(handle, "Path");
  path.def("__call__", &PathD::operator(), "lengths"_a)
      .def("lengths", &PathD::lengths, "points"_a)
      .def("tangent", &PathD::tangent, "lengths"_a)
      .def("normal", &PathD::normal, "lengths"_a)
      .def("angle0", &PathD::angle0, "lengths"_a)
      .def("angle1", &PathD::angle1, "lengths"_a)
      .def("angle2", &PathD::angle2, "lengths"_a);

  using PolychainD = Polychain<Eigen::Dynamic>;
  py::classh<PolychainD>(handle, "Polychain", path)
      .def(py::init<Eigen::ArrayXd, Eigen::ArrayXd>(), "x"_a, "y"_a)
      .def("setPoints", &PolychainD::setPoints, "x"_a, "y"_a)
      .def("__call__", &PolychainD::operator(), "lengths"_a)
      .def("lengths", &PolychainD::lengths, "points"_a);

  using TransformD = Transform<Eigen::Dynamic>;
  py::classh<TransformD>(handle, "Transform")
      .def(py::init<std::shared_ptr<PathD>>(), "path"_a)
      .def("posFrenet", &TransformD::posFrenet, "posCartes"_a)
      .def("posCartes", &TransformD::posCartes, "posFrenet"_a)
      .def("velFrenet", &TransformD::velFrenet, "velCartes"_a, "posFrenet"_a)
      .def("velCartes", &TransformD::velCartes, "velFrenet"_a, "posFrenet"_a)
      .def("accFrenet", &TransformD::accFrenet, "accCartes"_a, "velFrenet"_a,
           "posFrenet"_a)
      .def("accCartes", &TransformD::accCartes, "accFrenet"_a, "velFrenet"_a,
           "posFrenet"_a);
}
} // namespace FrenetTransform