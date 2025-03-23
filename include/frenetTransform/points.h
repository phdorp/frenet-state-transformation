#ifndef POINTS_H
#define POINTS_H

#include <Eigen/Core>
#include <cassert>

#include "frenetTransform/point.h"

namespace FrenetTransform {
/**
 * @brief Point sequence.
 *
 * @tparam NumPoints number of points along the sequence.
 * @tparam PointType type of points.
 */
template <int NumPoints, typename PointType = Point> class Points {
public:
  using ArrayPoints = Eigen::Array<double, NumPoints, 1>;

  Points() = default;

  /**
   * @brief Construct a new Points object from sequences of x- and
   * y-coordinates. The number of coordinates in both directions must be equal.
   *
   * @param x coordinates in x-direction of point sequence.
   * @param y coordinates in y-direction of point sequence.
   */
  Points(const ArrayPoints &x, const ArrayPoints &y) : m_x{x}, m_y{y} {}

  /**
   * @brief Get point at "index" along the point sequence.
   *
   * @param index of point along the point sequence.
   * @return PointType point at the "index".
   */
  PointType operator()(const int index) const {
    return {m_x(index), m_y(index)};
  }

  /**
   * @brief Provide the number of points along the sequence.
   *
   * @return int number of points.
   */
  int numPoints() const { return NumPoints > -1 ? NumPoints : m_x.rows(); }

  /**
   * @brief Get x-coordinate at "index" along the point sequence.
   *
   * @param index of x-coordinate along the point sequence.
   * @return double x-coordiante.
   */
  double x(const int index) const { return m_x(index); }

  /**
   * @brief Get y-coordinate at "index" along the point sequence.
   *
   * @param index of y-coordinate along the point sequence.
   * @return double y-coordiante.
   */
  double y(const int index) const { return m_y(index); }

  /**
   * @brief Get x-coordinates along the point sequence.
   *
   * @return const ArrayPoints& x-coordinates.
   */
  const ArrayPoints &x() const { return m_x; }

  /**
   * @brief Get y-coordinates along the point sequence.
   *
   * @return const ArrayPoints& y-coordinates.
   */
  const ArrayPoints &y() const { return m_y; }

  void setX(const ArrayPoints &x) { m_x = x; }

  void setY(const ArrayPoints &y) { m_y = y; }

  /**
   * @brief Squared distance between Points and "point".
   *
   * @param point to determine the squared distances to.
   * @return ArrayPoints squared distances to "point".
   */
  ArrayPoints distanceSquare(const PointType &point) {
    return (m_x - point.x()).pow(2) + (m_y - point.y()).pow(2);
  }

  /**
   * @brief Distance between Points and "point".
   *
   * @param point to determine the distances to.
   * @return ArrayPoints distances to "point".
   */
  ArrayPoints distance(const PointType &point) {
    return distanceSquare(point).sqrt();
  }

  /**
   * @brief Negate x- and y-coordinates.
   *
   * @return Points with negated x- and y-coordinates.
   */
  Points operator-() const { return {-m_x, -m_y}; }

  /**
   * @brief Element-wise sum of two Points.
   *
   * @param points1 first Points sequence of the sum.
   * @param points2 second Points sequence of the sum.
   * @return Points<NumPoints, PointType> with sum of coordinates of "points1"
   * and "points2".
   */
  friend Points<NumPoints, PointType>
  operator+(const Points<NumPoints, PointType> &points1,
            const Points<NumPoints, PointType> &points2) {
    return {points1.m_x + points2.m_x, points1.m_y + points2.m_y};
  }

  /**
   * @brief Element-wise sum of Points and Point.
   *
   * @param points Points involved in the sum.
   * @param point Point involved in the sum.
   * @return Points<NumPoints, PointType> with sum of coordinates in "point" and
   * "points".
   */
  friend Points<NumPoints, PointType>
  operator+(const Points<NumPoints, PointType> &points,
            const PointType &point) {
    return {point.m_x + point.x(), points.m_y + point.y()};
  }

  /**
   * @brief Element-wise sum of Point and Points.
   *
   * @param point Point involved in the sum.
   * @param points Points involved in the sum.
   * @return Points<NumPoints, PointType> with sum of coordinates in "point" and
   * "points".
   */
  friend Points<NumPoints, PointType>
  operator+(const PointType &point,
            const Points<NumPoints, PointType> &points) {
    return points + point;
  }

  /**
   * @brief Element-wise difference between two Points.
   *
   * @param points1 first Points sequence of the difference.
   * @param points2 second Points sequence of the difference.
   * @return Points<NumPoints, PointType> with difference between coordinates of
   * "points1" and "points2".
   */
  friend Points<NumPoints, PointType>
  operator-(const Points<NumPoints, PointType> &points1,
            const Points<NumPoints, PointType> &points2) {
    return -points2 + points1;
  }

  /**
   * @brief Element-wise difference between Points and Point.
   *
   * @param points Points involved in the difference.
   * @param point Point involved in the difference.
   * @return Points<NumPoints, PointType> with difference between coordinates in
   * "points" and "point".
   */
  friend Points<NumPoints, PointType>
  operator-(const Points<NumPoints, PointType> &points,
            const PointType &point) {
    return -point + points;
  }

  /**
   * @brief Element-wise difference between Point and Points.
   *
   * @param point Point involved in the difference.
   * @param points Points involved in the difference.
   * @return Points<NumPoints, PointType> with difference between coordinates in
   * "point" and "points".
   */
  friend Points<NumPoints, PointType>
  operator-(const PointType &point,
            const Points<NumPoints, PointType> &points) {
    return -points + point;
  }

  /**
   * @brief Element-wise scalar product of two Points.
   *
   * @param points1 first Points involved in the product.
   * @param points2 second Points involved in the product.
   * @return ArrayPoints scalar products of "points1" and "points2".
   */
  friend ArrayPoints operator*(const Points<NumPoints, PointType> &points1,
                               const Points<NumPoints, PointType> &points2) {
    return points1.m_x * points2.m_x + points1.m_y * points2.m_y;
  }

  /**
   * @brief Element-wise product of Points and scalars.
   *
   * @param points Points involved in the product.
   * @param nums scalars invovled in the product.
   * @return Points<NumPoints, PointType> products of "points" and "nums".
   */
  friend Points<NumPoints, PointType>
  operator*(const Points<NumPoints, PointType> &points,
            const ArrayPoints &nums) {
    return {points.m_x * nums, points.m_y * nums};
  }

  /**
   * @brief Element-wise product of scalars and Points.
   *
   * @param nums scalars invovled in the product.
   * @param points Points involved in the product.
   * @return Points<NumPoints, PointType> products of "points" and "nums".
   */
  friend Points<NumPoints, PointType>
  operator*(const ArrayPoints &nums,
            const Points<NumPoints, PointType> &points) {
    return points * nums;
  }

private:
  ArrayPoints m_x{}; /**<< x-coordinates of the Point sequence */
  ArrayPoints m_y{}; /**<< y-coordinates of the Point sequence */
};
}; // namespace FrenetTransform

#endif