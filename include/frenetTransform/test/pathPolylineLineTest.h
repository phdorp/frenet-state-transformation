#include "frenetTransform/polychain.h"
#include "frenetTransform/points.h"
#include "frenetTransform/transform.h"
#include "frenetTransform/internal/line.h"
#include "frenetTransform/test/testBase.h"

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

namespace FrenetTransform
{
    namespace Internal
    {
        class PathPolylineLineTest : public TestBase
        {
        protected:
            const Line<Eigen::Dynamic> m_line { {0.0, 0.0}, {1.0, 2.0} };

            const Points<Eigen::Dynamic> m_pointsCartes { Eigen::ArrayXd::Random(100).abs(), Eigen::ArrayXd::Random(100).abs() };

            const Polychain<Eigen::Dynamic> m_polyline { m_line(Eigen::ArrayXd {{0.0, 1.0, 1.5, 2.23}}) };
        };
    };
};