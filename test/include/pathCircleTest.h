#include "points.h"
#include "circle.h"
#include "transformCircle.h"
#include "testBase.h"

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

namespace FrenetTransform
{
    namespace Testing
    {
        class PathCircleTest : public TestBase
        {
        protected:
            const Circle m_circle { 5.0, {0.0, 0.0}, -M_PI };
            const TransformCircle m_transform { std::make_shared<Circle>(m_circle) };

            const int numQuery { 100 };
            const Points<Eigen::Dynamic, PointCircle> m_posCircle { m_circle.radius() * (1 + Eigen::ArrayXd::Random(numQuery) * 0.95) , Eigen::ArrayXd::Random(numQuery) * M_PI * 0.95 };
            const Points<Eigen::Dynamic, PointFrenet> m_posFrenet { m_transform.posFrenet(m_posCircle) };
        };
    };
};