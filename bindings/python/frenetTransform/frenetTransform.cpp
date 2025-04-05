#include "frenetTransform/point.h"
#include "frenetTransform/points.h"
#include "frenetTransform/polychain.h"
#include "frenetTransform/transform.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;

namespace FrenetTransform {
PYBIND11_MODULE(_core, handle) {
  py::classh<Point>(handle, "Point")
      .def(py::init<>())
      .def(py::init<double, double>())
      .def("x", &Point::x)
      .def("y", &Point::y)
      .def("distanceSquare", &Point::distanceSquare)
      .def("distance", &Point::distance)
      .def("__neg__", &Point::operator-)
      .def(py::self + py::self)
      .def(py::self - py::self);

  using PointsD = Points<Eigen::Dynamic>;
  py::classh<PointsD>(handle, "Points")
      .def(py::init<Eigen::ArrayXd, Eigen::ArrayXd>())
      .def("x", py::overload_cast<const int>(&PointsD::x, py::const_))
      .def("x", py::overload_cast<>(&PointsD::x, py::const_))
      .def("y", py::overload_cast<const int>(&PointsD::y, py::const_))
      .def("y", py::overload_cast<>(&PointsD::y, py::const_))
      .def("setX", &PointsD::setX)
      .def("setY", &PointsD::setY);

  using PathD = Path<Eigen::Dynamic>;
  py::classh<PathD> path(handle, "Path");
  path.def("__call__", &PathD::operator())
      .def("lengths", &PathD::lengths)
      .def("tangent", &PathD::tangent)
      .def("normal", &PathD::normal)
      .def("angle0", &PathD::angle0)
      .def("angle1", &PathD::angle1)
      .def("angle2", &PathD::angle2);

  using PolychainD = Polychain<Eigen::Dynamic>;
  py::classh<PolychainD>(handle, "Polychain", path)
      .def(py::init<Eigen::ArrayXd, Eigen::ArrayXd>())
      .def("setPoints", &PolychainD::setPoints)
      .def("__call__", &PolychainD::operator())
      .def("lengths", &PolychainD::lengths);

  using TransformD = Transform<Eigen::Dynamic>;
  py::classh<TransformD>(handle, "Transform")
      .def(py::init<std::shared_ptr<PathD>>())
      .def("posFrenet", &TransformD::posFrenet)
      .def("posCartes", &TransformD::posCartes)
      .def("velFrenet", &TransformD::velFrenet)
      .def("velCartes", &TransformD::velCartes)
      .def("accFrenet", &TransformD::accFrenet)
      .def("accCartes", &TransformD::accCartes);
}
} // namespace FrenetTransform