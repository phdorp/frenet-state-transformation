#include "polyline.h"
#include "points.h"
#include "transform.h"
#include "circle.h"
#include "pathCircleTest.h"

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include <math.h>
#include <memory>

namespace FrenetTransform
{
    namespace Testing
    {
        class PathPolylineTest : public PathCircleTest
        {
        protected:
            const Polyline<Eigen::Dynamic> m_circlePoly { m_circle(Eigen::ArrayXd::LinSpaced(1079, 0.0, 2 * M_PI) * m_circle.radius()) };
            const Transform m_circleTransform { std::make_shared<Polyline<Eigen::Dynamic>>(m_circlePoly) };

            const Points<Eigen::Dynamic, PointCartes> m_posCartes { m_transform.posCartes(m_posCircle) };

            const Points<Eigen::Dynamic, PointCircle> m_velCircle { Eigen::ArrayXd::Random(100).abs() * m_circle.radius(), Eigen::ArrayXd::Random(100) * M_PI / 4 };
            const Points<Eigen::Dynamic, PointCartes> m_velCartes { m_transform.velCartes(m_velCircle, m_posCircle) };
            const Points<Eigen::Dynamic, PointFrenet> m_velFrenet { m_transform.velFrenet(m_velCircle) };

            const Points<Eigen::Dynamic, PointCircle> m_accCircle { Eigen::ArrayXd::Random(100).abs() * m_circle.radius(), Eigen::ArrayXd::Random(100) * M_PI / 4 };
            const Points<Eigen::Dynamic, PointCartes> m_accCartes { m_transform.accCartes(m_accCircle, m_velCircle, m_posCircle) };
            const Points<Eigen::Dynamic, PointFrenet> m_accFrenet { m_transform.accFrenet(m_accCircle) };
        };

        TEST_F(PathPolylineTest, GetPointsCircle)
        {
            const auto pointsCircleEst { m_circlePoly(m_posFrenet.x()) };
            const auto pointsCircleGtr { m_circle(m_posFrenet.x()) };

            expectAllClose(pointsCircleEst.x(), pointsCircleGtr.x(), 1e-2);
            expectAllClose(pointsCircleEst.y(), pointsCircleGtr.y(), 1e-2);
        }

        TEST_F(PathPolylineTest, GetTangentsCircle)
        {
            const auto tangentsEst { m_circlePoly.tangent(m_posFrenet.x()) };
            const auto tangentsGtr { m_circle.tangent(m_posFrenet.x()) };

            expectAllClose(tangentsEst.x(), tangentsGtr.x(), 1e-2);
            expectAllClose(tangentsEst.y(), tangentsGtr.y(), 1e-2);
        }

        TEST_F(PathPolylineTest, GetNormalsCircle)
        {
            const auto normalsEst { m_circlePoly.normal(m_posFrenet.x()) };
            const auto normalsGtr { m_circle.normal(m_posFrenet.x()) };

            expectAllClose(normalsEst.x(), normalsGtr.x(), 1e-2);
            expectAllClose(normalsEst.y(), normalsGtr.y(), 1e-2);
        }

        TEST_F(PathPolylineTest, GetAnglesCircle)
        {
            const auto angles0Est { m_circlePoly.angle0(m_posFrenet.x()) };
            const auto angles0Gtr { m_circle.angle0(m_posFrenet.x()) };

            expectAllClose(angles0Est, angles0Gtr, 1e-2);
        }

        TEST_F(PathPolylineTest, GetAngles1Circle)
        {
            const auto angles1Est { m_circlePoly.angle1(m_posFrenet.x()) };
            const auto angles1Gtr { m_circle.angle1(m_posFrenet.x()) };

            expectAllClose(angles1Est, angles1Gtr, 1e-2);
        }

        TEST_F(PathPolylineTest, GetAngles2Circle)
        {
            const auto angles2Est { m_circlePoly.angle2(m_posFrenet.x()) };
            const auto angles2Gtr { m_circle.angle2(m_posFrenet.x()) };

            expectAllClose(angles2Est,angles2Gtr, 1e-2);
        }

        TEST_F(PathPolylineTest, NextPointsCircle)
        {
            const auto lengthsEst { m_circlePoly.lengths(m_posCartes) };
            const auto lengthsGtr { m_circle.lengths(m_posCartes) };

            expectAllClose(lengthsEst, lengthsGtr, 2e-2);
        }

        TEST_F(PathPolylineTest, PosFrenetCircle)
        {
            // determine frenet positions
            const auto posFrenet { m_circleTransform.posFrenet(m_posCartes) };

            expectAllClose(posFrenet.x(), m_posFrenet.x(), 1e-1, 2e-2);
            expectAllClose(posFrenet.y(), m_posFrenet.y(), 1e-1, 1e-2);
        }

        TEST_F(PathPolylineTest, PosCartCircle)
        {
            // determine frenet positions
            const auto posCartes { m_circleTransform.posCartes(m_posFrenet) };

            expectAllClose(posCartes.x(), m_posCartes.x(), 2e-2, 0.7);
            expectAllClose(posCartes.y(), m_posCartes.y(), 2e-2, 0.7);
        }

        TEST_F(PathPolylineTest, VelFrenetCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto velFrenet { m_circleTransform.velFrenet(m_velCartes, m_posFrenet) };

            expectAllClose(velFrenet.x(), m_velFrenet.x(), 2e-1, 0.55);
            expectAllClose(velFrenet.y(), m_velFrenet.y(), 2e-2, 3e-2);
        }

        TEST_F(PathPolylineTest, VelCartCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto velCartes { m_circleTransform.velCartes(m_velFrenet, m_posFrenet) };

            expectAllClose(velCartes.x(), m_velCartes.x(), 6e-2, 0.2);
            expectAllClose(velCartes.y(), m_velCartes.y(), 2e-2, 0.15);
        }

        TEST_F(PathPolylineTest, AccFrenetCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto accFrenet { m_circleTransform.accFrenet(m_accCartes, m_velFrenet, m_posFrenet) };

            expectAllClose(accFrenet.x(), m_accFrenet.x(), 7e-2, 4.5e-1);
            expectAllClose(accFrenet.y(), m_accFrenet.y(), 7e-2, 4.5e-1);
        }

        TEST_F(PathPolylineTest, AccCartCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto accCartes { m_circleTransform.accCartes(m_accFrenet, m_velFrenet, m_posFrenet) };

            expectAllClose(accCartes.x(), m_accCartes.x(), 2e-2, 2e-1);
            expectAllClose(accCartes.y(), m_accCartes.y(), 2e-2, 2e-1);
        }
    };
};

int main(int argc, char **argv)
{
    std::srand(0);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
