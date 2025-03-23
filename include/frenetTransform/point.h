#ifndef POINT_H
#define POINT_H

#include <math.h>

namespace FrenetTransform {
/**
 * @brief Representation of a 2-dimensional point.
 *
 */
class Point {
public:
  Point() = default;

  /**
   * @brief Construct a new Point with coordinates in x- and y-direction.
   *
   * @param x point's coordinate in x-direction.
   * @param y point's coordinate in y-direction.
   */
  Point(double x, double y) : m_x{x}, m_y{y} {}

  /**
   * @brief Provides the point's coordinate in x-direction.
   *
   * @return double x-coordinate.
   */
  double x() const { return m_x; }

  /**
   * @brief Provides the point's coordinate in y-direction.
   *
   * @return double y-coordinate.
   */
  double y() const { return m_y; }

  /**
   * @brief Squared distance between Point and a query "point".
   *
   * @param point to determine the squared distance to.
   * @return double squared distance to "point".
   */
  double distanceSquare(const Point &point) const {
    return std::pow(point.x() - m_x, 2) + std::pow(point.y() - m_y, 2);
  }

  /**
   * @brief Distance between Point and a query "point".
   *
   * @param point to determine the distance to.
   * @return double distance to "point".
   */
  double distance(const Point &point) const {
    return std::sqrt(distanceSquare(point));
  }

  /**
   * @brief Negate x- and y-coordinates.
   *
   * @return Point with negated x- and y-coordinates.
   */
  Point operator-() const { return {-m_x, -m_y}; }

  /**
   * @brief Sum between coordinates of two points.
   *
   * @param point1 first point of the sum.
   * @param point2 second point of the sum.
   * @return Point with sum of coordinates of "point1" and "point2".
   */
  friend Point operator+(const Point &point1, const Point &point2) {
    return {point1.x() + point2.x(), point1.y() + point2.y()};
  }

  /**
   * @brief Difference between coordinates of two points.
   *
   * @param point1 first point of the difference.
   * @param point2 second point of the difference.
   * @return Point with difference of coordinates of "point1" and "point2".
   */
  friend Point operator-(const Point &point1, const Point &point2) {
    return {point1 + (-point2)};
  }

private:
  double m_x{}; /**<< x-coordinate */
  double m_y{}; /**<< y-coordinate */
};
}; // namespace FrenetTransform

#endif