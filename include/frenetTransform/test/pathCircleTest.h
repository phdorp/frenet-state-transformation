#include "frenetTransform/internal/circle.h"
#include "frenetTransform/internal/transformCircle.h"
#include "frenetTransform/points.h"
#include "frenetTransform/test/testBase.h"

#include <Eigen/Core>
#include <gtest/gtest.h>
#include <math.h>
#include <memory>

namespace FrenetTransform {
namespace Internal {
/**
 * @brief Class for testing polyline and transformation functions against circle
 * as ground truth.
 *
 * @tparam NumQueries number of query points to transform, -1 is dynamic.
 */
template <typename NumQueries> class PathCircleTest : public TestBase {
public:
  static constexpr int s_numQueries{
      NumQueries::s_val == -1 ? 100 : NumQueries::s_val}; /**<< number of query points for point generation */

protected:
  const Circle<NumQueries::s_val> m_circle{5.0, {0.0, 0.0}, -M_PI}; /**<< Circle of radius 5, starting at -PI */
  const TransformCircle<NumQueries::s_val> m_transform{
      std::make_shared<Circle<NumQueries::s_val>>(m_circle)}; /**<< Circle Transform as transformation ground truth */

  const Points<NumQueries::s_val, PointCircle> m_posCircle{
      m_circle.radius() *
          (1 +
           Eigen::Array<double, NumQueries::s_val, 1>::Random(s_numQueries) *
               0.95),
      Eigen::Array<double, NumQueries::s_val, 1>::Random(s_numQueries) * M_PI *
          0.95}; /**<< "s_numQueries" number of points in Circle frame */
  const Points<NumQueries::s_val> m_posFrenet{
      m_transform.posFrenet(m_posCircle)}; /**<< "m_posCircle" in Frenet frame */
};
}; // namespace Internal
}; // namespace FrenetTransform