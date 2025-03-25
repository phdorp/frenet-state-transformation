#include "frenetTransform/internal/line.h"
#include "frenetTransform/points.h"
#include "frenetTransform/polychain.h"
#include "frenetTransform/test/testBase.h"
#include "frenetTransform/transform.h"

#include <Eigen/Core>
#include <gtest/gtest.h>
#include <math.h>
#include <memory>

namespace FrenetTransform {
namespace Internal {
/**
 * @brief Class for testing polyline and transformation functions against line
 * as ground truth.
 *
 */
class PathPolylineLineTest : public TestBase {
protected:
  const Line<Eigen::Dynamic> m_line{
      {0.0, 0.0}, {1.0, 2.0}}; /**<< Line between (0, 0) and (1, 2) */

  const Points<Eigen::Dynamic> m_pointsCartes{
      Eigen::ArrayXd::Random(100).abs(),
      Eigen::ArrayXd::Random(100)
          .abs()}; /**<< 100 random points between in unit square */

  const Polychain<Eigen::Dynamic> m_polyline{m_line(
      Eigen::ArrayXd{{0.0, 1.0, 1.5, 2.23}})}; /**<< Polychain along the Line */
};
}; // namespace Internal
}; // namespace FrenetTransform