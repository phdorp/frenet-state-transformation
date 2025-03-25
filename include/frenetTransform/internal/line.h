#ifndef LINE_H
#define LINE_H

#include <Eigen/Core>
#include <math.h>

#include "frenetTransform/internal/math.h"
#include "frenetTransform/path.h"
#include "frenetTransform/point.h"

namespace FrenetTransform {
namespace Internal {
/**
 * @brief Path representation as as line.
 * The class is thought to provide ground truth for the Polychain
 * implementation.
 *
 * @tparam NumQueries number of query points with -1 for dynamic point number.
 */
template <int NumQueries = Eigen::Dynamic>
class Line : public Path<NumQueries> {
public:
  using ArrayQueries = Eigen::Array<double, NumQueries, 1>;

  Line() = delete;

  /**
   * @brief Construct a new Line object from its radius, center and angle
   * offset.
   *
   * @param start start point of the Line.
   * @param end end point of the Line.
   */
  Line(const Point &start, const Point &end) : m_start{start}, m_end{end} {}

  /**
   * @brief Gets points along the Line at the query lengths.
   * Lengths exceeding the Line's domain either resolve to the first or
   * last point.
   *
   * @param lengths query lengths along the Circle.
   * @return Points<NumQueries> at the query lengths.
   */
  Points<NumQueries> operator()(const ArrayQueries &lengths) const override {
    // determine relative lengths along line
    const double lineLength{m_end.distance(m_start)};
    ArrayQueries relLengths{lengths / lineLength};

    // clamp lengths to line length
    for (auto &relLength : relLengths)
      relLength = std::clamp(relLength, 0.0, 1.0);

    return {m_start.x() * relLengths + (1 - relLengths) * m_end.x(),
            m_start.y() * relLengths + (1 - relLengths) * m_end.y()};
  }

  /**
   * @brief Determines distances along the Line to next points to the query
   * points.
   *
   * @param points query points.
   * @return Points<NumQueries> next points to query points.
   */
  ArrayQueries lengths(const Points<NumQueries> &points) const override {
    // determine relative lengths of points along the line
    const auto xDiffPnt{m_end.x() - points.x()};
    const auto yDiffPnt{m_end.y() - points.y()};
    const Point pointDiff{m_end - m_start};

    ArrayQueries relLengths{
        (xDiffPnt * pointDiff.x() + yDiffPnt * pointDiff.y()) /
        (std::pow(pointDiff.x(), 2) + std::pow(pointDiff.y(), 2))};

    // clamp lengths to line length
    for (auto &relLength : relLengths)
      relLength = std::clamp(relLength, 0.0, 1.0);

    return relLengths * m_end.distance(m_start);
  }

private:
  const Point m_start; /*<< start point of the line */
  const Point m_end;   /*<< end point of the line */

  /**
   * @brief Determines 1st order gradient at the given path lengths.
   *
   * @param lengths lengths along the path.
   * @return 1st order gradient at given path lengths.
   */
  Points<NumQueries> gradient1(const ArrayQueries &lengths) const override {
    ArrayQueries gradX(lengths.rows());
    gradX += m_end.x() - m_start.x() / m_end.distance(m_start);
    ArrayQueries gradY(lengths.rows());
    gradY += m_end.y() - m_start.y() / m_end.distance(m_start);
    return {gradX, gradY};
  }

  /**
   * @brief Determines 2nd order gradient at the given path lengths.
   *
   * @param lengths lengths along the path.
   * @return 2nd order gradient at given path lengths.
   */
  Points<NumQueries> gradient2(const ArrayQueries &lengths) const override {
    return {ArrayQueries::Zero(lengths.rows()),
            ArrayQueries::Zero(lengths.rows())};
  }

  /**
   * @brief Determines 3rd order gradient at the given path lengths.
   *
   * @param lengths lengths along the path.
   * @return 3rd order gradient at given path lengths.
   */
  Points<NumQueries> gradient3(const ArrayQueries &lengths) const override {
    return {ArrayQueries::Zero(lengths.rows()),
            ArrayQueries::Zero(lengths.rows())};
  }
};
}; // namespace Internal
}; // namespace FrenetTransform

#endif