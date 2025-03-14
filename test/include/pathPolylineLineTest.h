#include "frenetTransform/polyline.h"
#include "frenetTransform/points.h"
#include "frenetTransform/transform.h"
#include "line.h"
#include "testBase.h"

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

namespace FrenetTransform
{
    namespace Testing
    {
        class PathPolylineLineTest : public TestBase
        {
        protected:
            const Line m_line { {0.0, 0.0}, {1.0, 2.0} };

            const Points<Eigen::Dynamic, PointCartes> m_pointsCartes { Eigen::ArrayXd::Random(100).abs(), Eigen::ArrayXd::Random(100).abs() };

            const Polyline<-1> m_polyline { m_line(Eigen::ArrayXd {{0.0, 1.0, 1.5, 2.23}}) };
        };
    };
};