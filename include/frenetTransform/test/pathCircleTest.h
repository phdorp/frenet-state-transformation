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
template <typename NumQueries> class PathCircleTest : public TestBase {
public:
  static constexpr int s_numQueries{
      NumQueries::s_val == -1 ? 100 : NumQueries::s_val};

protected:
  const Circle<NumQueries::s_val> m_circle{5.0, {0.0, 0.0}, -M_PI};
  const TransformCircle<NumQueries::s_val> m_transform{
      std::make_shared<Circle<NumQueries::s_val>>(m_circle)};

  const Points<NumQueries::s_val, PointCircle> m_posCircle{
      m_circle.radius() *
          (1 +
           Eigen::Array<double, NumQueries::s_val, 1>::Random(s_numQueries) *
               0.95),
      Eigen::Array<double, NumQueries::s_val, 1>::Random(s_numQueries) * M_PI *
          0.95};
  const Points<NumQueries::s_val> m_posFrenet{
      m_transform.posFrenet(m_posCircle)};
};
}; // namespace Internal
}; // namespace FrenetTransform