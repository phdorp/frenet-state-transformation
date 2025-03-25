#ifndef TRANSFORM_CIRCLE_H
#define TRANSFORM_CRICLE_H

#include "frenetTransform/internal/circle.h"
#include "frenetTransform/point.h"
#include "frenetTransform/transform.h"

namespace FrenetTransform {
namespace Internal {
  /**
   * @brief Point in Circle frame.
   * Used to overload the functions in Transform Circle to transform from Circle to Frenet and Cartesian frames.
   *
   */
class PointCircle : public Point {
public:
  PointCircle(double radius, double angle) : Point(radius, angle) {}
};

/**
 * @brief Transformation between Cartesian, Frenet, and Circle frame.
 * Thought to provide ground truth for the Transformation.
 *
 * @tparam NumQueries number of query points with -1 for dynamic point number.
 */
template <int NumQueries = Eigen::Dynamic>
class TransformCircle : public Transform<NumQueries> {
public:
  /**
   * @brief Construct a new Transform Circle from a given Circle.
   *
   * @param circle defines the Frenet frame.
   */
  TransformCircle(std::shared_ptr<Circle<NumQueries>> circle)
      : Transform<NumQueries>(circle), m_path{circle} {}

  /**
   * @brief Transform Circle positions to Frenet positions.
   *
   * @param posCircle Query points in Circle coordinates.
   * @return Points<NumQueries> Result points in Frenet coordinates.
   */
  Points<NumQueries>
  posFrenet(const Points<NumQueries, PointCircle> &posCircle) const {
    return {m_path->radius() * (posCircle.y() - m_path->angleOffset()),
            m_path->radius() - posCircle.x()};
  }

  /**
   * @brief Transform Circle velocities to Frenet velocities.
   *
   * @param velCircle Query velocities in Circle coordinates.
   * @return Points<NumQueries> Result velocities in Frenet coordinates.
   */
  Points<NumQueries>
  velFrenet(const Points<NumQueries, PointCircle> &velCircle) const {
    return {m_path->radius() * velCircle.y(), -velCircle.x()};
  }

  /**
   * @brief Transform Circle accelerations to Frenet accelerations.
   *
   * @param accCircle Query accelerations in Circle coordinates.
   * @return Points<NumQueries> Result accelerations in Frenet coordinates.
   */
  Points<NumQueries>
  accFrenet(const Points<NumQueries, PointCircle> &accCircle) const {
    return {m_path->radius() * accCircle.y(), -accCircle.x()};
  }

  /**
   * @brief Transform Circle positions to Cartesian positions.
   *
   * @param posCircle Query points in Circle coordinates.
   * @return Points<NumQueries> Result points in Cartesian coordinates.
   */
  Points<NumQueries>
  posCartes(const Points<NumQueries, PointCircle> &posCircle) const {
    return {posCircle.x() * posCircle.y().cos() + m_path->center().x(),
            posCircle.x() * posCircle.y().sin() + m_path->center().y()};
  }

  /**
   * @brief Transform Circle velocities to Cartesian velocities.
   *
   * @param velCircle Query velocities in Circle coordinates.
   * @param posCircle Query positions in Circle coordinates.
   * @return Points<NumQueries> Result velocities in Cartesian coordinates.
   */
  Points<NumQueries>
  velCartes(const Points<NumQueries, PointCircle> &velCircle,
            const Points<NumQueries, PointCircle> &posCircle) const {
    return {velCircle.x() * posCircle.y().cos() -
                posCircle.x() * velCircle.y() * posCircle.y().sin(),
            velCircle.x() * posCircle.y().sin() +
                posCircle.x() * velCircle.y() * posCircle.y().cos()};
  }

  /**
   * @brief Transform Circle accelerations to Cartesian accelerations.
   *
   * @param accCircle Query accelerations in Circle coordinates.
   * @param velCircle Query velocities in Circle coordinates.
   * @param posCircle Query positions in Circle coordinates.
   * @return Points<NumQueries> Result accelerations in Cartesian coordinates.
   */
  Points<NumQueries>
  accCartes(const Points<NumQueries, PointCircle> &accCircle,
            const Points<NumQueries, PointCircle> &velCircle,
            const Points<NumQueries, PointCircle> &posCircle) const {
    return {(accCircle.x() - posCircle.x() * velCircle.y().pow(2)) *
                    posCircle.y().cos() -
                (2 * velCircle.x() * velCircle.y() +
                 posCircle.x() * accCircle.y()) *
                    posCircle.y().sin(),
            (accCircle.x() - posCircle.x() * velCircle.y().pow(2)) *
                    posCircle.y().sin() +
                (2 * velCircle.x() * velCircle.y() +
                 posCircle.x() * accCircle.y()) *
                    posCircle.y().cos()};
  }

private:
  const std::shared_ptr<Circle<NumQueries>> m_path;
};
}; // namespace Internal
}; // namespace FrenetTransform

#endif