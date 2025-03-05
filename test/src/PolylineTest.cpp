#include "polyline.h"
#include "points.h"
#include "transform.h"
#include "circle.h"
#include "CircleTest.h"

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include <math.h>
#include <memory>

namespace FrenetTransform
{
    namespace Testing
    {
        class PolylineTest : public CircleTest
        {
        protected:
            const Polyline<5> m_straight { {0.0, 1.0, 2.0, 4.0, 7.0},  {0.0, 1.0, 2.0, 4.0, 7.0} };

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

        TEST_F(PolylineTest, GetPointsStraight)
        {
            Eigen::ArrayXd input { {0.0, std::sqrt(8), std::sqrt(98)} };

            const auto result { m_straight(input) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth { {{0.0, 2.0, 7.0}}, {{0.0, 2.0, 7.0}} };

            for(int index {}; index < input.rows(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-10);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-10);
            }
        }

        // TEST_F(PolylineTest, GetPointsCircle)
        // {
        //     const auto result { m_circlePoly(m_posFrenet.x()) };

        //     const Points<Eigen::Dynamic, PointCartes> groundTruth { m_circle(m_posFrenet.x()) };

        //     for(int index {}; index < m_posCartes.numPoints(); ++index)
        //     {
        //         EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-2);
        //         EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-2);
        //     }
        // }

        TEST_F(PolylineTest, GetTangentsCircle)
        {
            const auto result { m_circlePoly.tangent(m_posFrenet.x()) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth { m_circle.tangent(m_posFrenet.x()) };

            for(int index {}; index < groundTruth.numPoints(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-2);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, GetNormalsCircle)
        {
            const auto result { m_circlePoly.normal(m_posFrenet.x()) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth { m_circle.normal(m_posFrenet.x()) };

            for(int index {}; index < groundTruth.numPoints(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-2);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, GetAnglesCircle)
        {
            const auto result { m_circlePoly.angle0(m_posFrenet.x()) };

            const Eigen::ArrayXd groundTruth { m_circle.angle0(m_posFrenet.x()) };

            for(int index {}; index < groundTruth.rows(); ++index)
                EXPECT_NEAR(groundTruth(index), result(index), 1e-2);
        }

        TEST_F(PolylineTest, GetAngles1Circle)
        {
            const auto result { m_circlePoly.angle1(m_posFrenet.x()) };

            const Eigen::ArrayXd groundTruth { m_circle.angle1(m_posFrenet.x()) };

            for(int index {}; index < groundTruth.rows(); ++index)
            {
                EXPECT_NEAR(groundTruth(index), result(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, GetAngles2Circle)
        {
            const auto result { m_circlePoly.angle2(m_posFrenet.x()) };

            const Eigen::ArrayXd groundTruth { m_circle.angle2(m_posFrenet.x()) };

            for(int index {}; index < groundTruth.rows(); ++index)
            {
                EXPECT_NEAR(groundTruth(index), result(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, NextPointsStraight)
        {
            Points<Eigen::Dynamic, PointCartes> input {
                {{0.0, 3.0, 2.0}},
                {{0.0, 3.0, 4.0}}
            };

            const auto lengths { m_straight.lengths(input) };
            const auto result { m_straight(lengths) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth {
                {{0.0, 3.0, 3.0}},
                {{0.0, 3.0, 3.0}}
            };

            for(int index {}; index < input.x().rows(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-10);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-10);
            }
        }

        TEST_F(PolylineTest, NextPointsCircle)
        {
            const auto lengths { m_circlePoly.lengths(m_posCartes) };
            const auto result { m_circlePoly(lengths) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth { m_circle(m_posFrenet.x()) };

            for(int index {}; index < groundTruth.x().rows(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 2e-2);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 2e-2);
            }
        }

        TEST_F(PolylineTest, PosFrenetCircle)
        {
            // determine frenet positions
            const auto result { m_circleTransform.posFrenet(m_posCartes) };

            // determine ground truth frenet positions
            const Points<Eigen::Dynamic, PointFrenet> groundTruth  { m_posFrenet };

            for(int index {}; index < groundTruth.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 2e-2);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 1e-2);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 1e-1);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 1e-1);
            }
        }

        TEST_F(PolylineTest, PosCartCircle)
        {
            // determine frenet positions
            const auto result { m_circleTransform.posCartes(m_posFrenet) };

            // determine ground truth cartesian positions
            const Points<Eigen::Dynamic, PointCartes> groundTruth { m_posCartes };

            for(int index {}; index < groundTruth.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 0.7);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 0.7);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 2e-2);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 2e-2);
            }
        }

        TEST_F(PolylineTest, VelFrenetCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto result { m_circleTransform.velFrenet(m_velCartes, m_posFrenet) };

            const Points<Eigen::Dynamic, PointFrenet> groundTruth { m_velFrenet };

            for(int index {}; index < groundTruth.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 0.55);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 3e-2);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 0.2);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 2e-2);
            }
        }

        TEST_F(PolylineTest, VelCartCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto result { m_circleTransform.velCartes(m_velFrenet, m_posFrenet) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth { m_velCartes };

            for(int index {}; index < groundTruth.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 0.2);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 0.15);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 6e-2);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 2e-2);
            }
        }

        TEST_F(PolylineTest, AccFrenetCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto result { m_circleTransform.accFrenet(m_accCartes, m_velFrenet, m_posFrenet) };

            const Points<Eigen::Dynamic, PointFrenet> groundTruth { m_accFrenet };

            for(int index {}; index < groundTruth.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 4.5e-1);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 4.5e-1);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 7e-2);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 7e-2);
            }
        }

        TEST_F(PolylineTest, AccCartCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto result { m_circleTransform.accCartes(m_accFrenet, m_velFrenet, m_posFrenet) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth { m_accCartes };

            for(int index {}; index < groundTruth.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 2e-1);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 2e-1);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 2e-2);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 2e-2);
            }
        }
    };
};

int main(int argc, char **argv)
{
    std::srand(0);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
