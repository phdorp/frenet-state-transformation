#ifndef CIRCLE_H
#define CIRCLE_H

#include <Eigen/Core>
#include <math.h>

#include "frenetTransform/internal/math.h"
#include "frenetTransform/path.h"
#include "frenetTransform/point.h"

namespace FrenetTransform {
namespace Internal {
/**
 * @brief Path representation as a Circle.
 * The class is thought to provide ground truth for the Polychain
 * implementation.
 *
 * @tparam NumQueries number of query points with -1 for dynamic point number.
 */
template <int NumQueries = Eigen::Dynamic>
class Circle : public Path<NumQueries> {
public:
  using ArrayQueries = Eigen::Array<double, NumQueries, 1>;

  Circle() = delete;

  /**
   * @brief Construct a new Circle object from its radius, center and angle
   * offset.
   *
   * @param radius of the Circle.
   * @param center of the Circle.
   * @param angle0 angle offset.
   */
  Circle(double radius, Point center, double angle0 = 0)
      : m_radius{radius}, m_center{center}, m_angle0{angle0} {}

  /**
   * @brief Gets points along the Circle at the query lengths.
   * Lengths exceeding the polychain's domain are wrapped at 2PI.
   *
   * @param lengths query lengths along the Circle.
   * @return Points<NumQueries> at the query lengths.
   */
  Points<NumQueries> operator()(const ArrayQueries &lengths) const override {
    return {m_center.x() + m_radius * angle(lengths).cos(),
            m_center.y() + m_radius * angle(lengths).sin()};
  }

  /**
   * @brief Determines distances along the Circle to next points to the query
   * points.
   *
   * @param points query points.
   * @return Points<NumQueries> next points to query points.
   */
  ArrayQueries lengths(const Points<NumQueries> &points) const override {
    const ArrayQueries dirsx{points.x() - m_center.x()};
    const ArrayQueries dirsy{points.y() - m_center.y()};
    return m_radius * (angleDir(dirsx, dirsy) - m_angle0);
  }

  /**
   * @brief Determines distances along the Circle to the query angles.
   *
   * @param angles query angles.
   * @return ArrayQueries distances to next points on Circle.
   */
  ArrayQueries lengths(const ArrayQueries &angles) const {
    return angles * m_radius;
  }

  /**
   * @brief Determines angles along the Circle to the given lengths.
   *
   * @param lengths query lengths.
   * @return ArrayQueries distances to next points on Circle.
   */
  ArrayQueries angle(const ArrayQueries &lengths) const {
    return lengths / m_radius + m_angle0;
  }

  /**
   * @brief Gets the Circle radius.
   *
   * @return double Circle radius.
   */
  double radius() const { return m_radius; }

  /**
   * @brief Gets the Circle center.
   *
   * @return const Point& Circle center.
   */
  const Point &center() const { return m_center; }

  /**
   * @brief gets the angle offset.
   *
   * @return double angle offset.
   */
  double angleOffset() const { return m_angle0; }

private:
  const double m_radius{}; /*<< Circle radius */
  const Point m_center;    /*<< Circle center */
  const double m_angle0{}; /*<< Circle offset angle */

  /**
   * @brief Determines 1st order gradient at the given path lengths.
   *
   * @param lengths lengths along the path.
   * @return 1st order gradient at given path lengths.
   */
  Points<NumQueries> gradient1(const ArrayQueries &lengths) const override {
    return {-angle(lengths).sin(), angle(lengths).cos()};
  }

  /**
   * @brief Determines 2nd order gradient at the given path lengths.
   *
   * @param lengths lengths along the path.
   * @return 2nd order gradient at given path lengths.
   */
  Points<NumQueries> gradient2(const ArrayQueries &lengths) const override {
    return {-angle(lengths).cos() / m_radius, -angle(lengths).sin() / m_radius};
  }

  /**
   * @brief Determines 3rd order gradient at the given path lengths.
   *
   * @param lengths lengths along the path.
   * @return 3rd order gradient at given path lengths.
   */
  Points<NumQueries> gradient3(const ArrayQueries &lengths) const override {
    return {angle(lengths).sin() / std::pow(m_radius, 2),
            -angle(lengths).cos() / std::pow(m_radius, 2)};
  }
};
}; // namespace Internal
}; // namespace FrenetTransform

#endif