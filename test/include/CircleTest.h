#include "points.h"
#include "circle.h"

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include <math.h>
#include <memory>

using FrenetTransform::Points;

namespace Testing
{
    class CircleTest : public testing::Test
    {
    protected:
        const Circle m_circle { 5.0, {0.0, 0.0}, -M_PI };
        const Points<Eigen::Dynamic> m_pointsCircle { Eigen::ArrayXd::Zero(100), Eigen::ArrayXd::LinSpaced(100, 0, 2 * M_PI) };
    };
};