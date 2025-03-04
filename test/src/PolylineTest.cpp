#include "polyline.h"
#include "points.h"
#include "transform.h"
#include "circle.h"

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include <math.h>
#include <memory>

namespace FrenetTransform
{
    namespace Testing
    {
        class PolylineTest : public testing::Test
        {
            protected:
            PolylineTest()
                : m_straight { {0.0, 1.0, 2.0, 4.0, 7.0},  {0.0, 1.0, 2.0, 4.0, 7.0} }
                , m_circlePoly { m_radius * m_lengthsCircle.cos(), m_radius * m_lengthsCircle.sin() }
                , m_circleTransform { std::make_shared<FrenetTransform::Polyline<1079>>(m_circlePoly) }
                , m_circle { m_radius, PointCartes { 0.0, 0.0 }, -M_PI }
            {
            }

            const FrenetTransform::Polyline<5> m_straight;

            const double m_radius { 5.0 };
            const Eigen::Array<double, 1079, 1> m_lengthsCircle { Eigen::Array<double, 1079, 1>::LinSpaced(-M_PI, M_PI)} ;
            const FrenetTransform::Polyline<1079> m_circlePoly {};
            const FrenetTransform::Transform m_circleTransform;
            const Testing::Circle m_circle;
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

        TEST_F(PolylineTest, GetPointsCircle)
        {
            Eigen::ArrayXd input { {M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
            input *= m_radius;

            const auto result { m_circlePoly(input) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth
            {
                {{m_radius, m_radius / std::sqrt(2),      0.0}},
                {{ 0.0,     m_radius / std::sqrt(2), m_radius}}
            };

            for(int index {}; index < input.rows(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-2);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, GetTangentsCircle)
        {
            Eigen::ArrayXd input { {M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
            input *= m_radius;

            const auto result { m_circlePoly.tangent(input) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth
            {
                {{0.0, -std::sqrt(50) / 10, -1.0}},
                {{1.0,  std::sqrt(50) / 10,  0.0}}
            };

            for(int index {}; index < input.rows(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-2);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, GetNormalsCircle)
        {
            Eigen::ArrayXd input { {M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
            input *= m_radius;

            const auto result { m_circlePoly.normal(input) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth
            {
                {{-1.0, -std::sqrt(50) / 10,  0.0}},
                {{ 0.0, -std::sqrt(50) / 10, -1.0}}
            };

            for(int index {}; index < input.rows(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-2);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, GetAnglesCircle)
        {
            const Eigen::ArrayXd angles { Eigen::ArrayXd::Random(20).abs() * M_PI * 2 };
            const Eigen::ArrayXd input { angles * m_radius };

            const auto result { m_circlePoly.angle0(input) };

            const Eigen::ArrayXd groundTruth { m_circle.angle0(input) };

            for(int index {}; index < input.rows(); ++index)
                EXPECT_NEAR(groundTruth(index), result(index), 1e-2);
        }

        TEST_F(PolylineTest, GetAngles1Circle)
        {
            Eigen::ArrayXd input { {M_PI / 2, M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
            input *= m_radius;

            const auto result { m_circlePoly.angle1(input) };

            const Eigen::ArrayXd groundTruth { Eigen::ArrayXd::Ones(input.size()) / m_radius };

            for(int index {}; index < input.rows(); ++index)
            {
                EXPECT_NEAR(groundTruth(index), result(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, GetAngles2Circle)
        {
            Eigen::ArrayXd input { {M_PI / 2, M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
            input *= m_radius;

            const auto result { m_circlePoly.angle2(input) };

            const Eigen::ArrayXd groundTruth { Eigen::ArrayXd::Zero(input.size()) };

            for(int index {}; index < input.rows(); ++index)
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
            Points<Eigen::Dynamic, PointCartes> input {
                {{  0.0,  0.0, 1.0}},
                {{-m_radius, -2.0, 1.0}}
            };

            const auto lengths { m_circlePoly.lengths(input) };
            const auto result { m_circlePoly(lengths) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth {
                {{  0.0,   0.0, m_radius / std::sqrt(2)}},
                {{-m_radius, -m_radius, m_radius / std::sqrt(2)}}
            };

            for(int index {}; index < input.x().rows(); ++index)
            {
                EXPECT_NEAR(groundTruth.x()(index), result.x()(index), 1e-2);
                EXPECT_NEAR(groundTruth.y()(index), result.y()(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, PosFrenetCircle)
        {
            // sample radii and angles
            const Eigen::Array<double, Eigen::Dynamic, 1> radii { Eigen::Array<double, Eigen::Dynamic, 1>::Random(20).abs() * m_radius };
            const Eigen::Array<double, Eigen::Dynamic, 1> angles { Eigen::Array<double, Eigen::Dynamic, 1>::Random(20) * M_PI };

            // determine cartesian positions
            const Points<Eigen::Dynamic, PointCartes> input { radii * angles.cos(), radii * angles.sin() };

            // determine frenet positions
            const auto result { m_circleTransform.posFrenet(input) };

            // determine ground truth frenet positions
            const Points<Eigen::Dynamic, PointFrenet> groundTruth  {m_radius * (angles + M_PI), m_radius - radii };

            for(int index {}; index < input.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 1e-2);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 1e-2);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 1e-1);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 1e-1);
            }
        }

        TEST_F(PolylineTest, PosCartCircle)
        {
            // sample radii and angles
            const Eigen::Array<double, Eigen::Dynamic, 1> radii { Eigen::Array<double, Eigen::Dynamic, 1>::Random(20).abs() * m_radius };
            const Eigen::Array<double, Eigen::Dynamic, 1> angles { Eigen::Array<double, Eigen::Dynamic, 1>::Random(20) * M_PI };

            // determine frenet positions
            const Points<Eigen::Dynamic, PointFrenet> input  {m_radius * (angles + M_PI), m_radius - radii };

            // determine frenet positions
            const auto result { m_circleTransform.posCartes(input) };

            // determine ground truth cartesian positions
            const Points<Eigen::Dynamic, PointCartes> groundTruth { radii * angles.cos(), radii * angles.sin() };

            for(int index {}; index < input.numPoints(); ++index)
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
            // sample radii and angles
            const Eigen::Array<double, Eigen::Dynamic, 1> radii { Eigen::Array<double, Eigen::Dynamic, 1>::Random(20).abs() * m_radius };
            const Eigen::Array<double, Eigen::Dynamic, 1> angles { Eigen::Array<double, Eigen::Dynamic, 1>::Random(20) * M_PI };

            // sample radial and angular velocities
            const Eigen::Array<double, Eigen::Dynamic, 1> radiiDt { Eigen::Array<double, Eigen::Dynamic, 1>::Random(20) * m_radius };
            const Eigen::Array<double, Eigen::Dynamic, 1> anglesDt { Eigen::Array<double, Eigen::Dynamic, 1>::Random(20) * M_PI / 8 };

            // determine cartesian positions
            const Points<Eigen::Dynamic, PointFrenet> posFrenet {m_radius * (angles + M_PI), m_radius - radii };

            // determine cartesian velocities
            const Points<Eigen::Dynamic, PointCartes> input { radiiDt * angles.cos() - radii * anglesDt * angles.sin(), radiiDt * angles.sin() + radii * anglesDt * angles.cos() };

            // test frenet frame results since error grows with distance from path
            const auto result { m_circleTransform.velFrenet(input, posFrenet) };

            const Points<Eigen::Dynamic, PointFrenet> groundTruth { m_radius * anglesDt, - radiiDt };

            for(int index {}; index < posFrenet.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 0.4);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 1e-2);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 0.2);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 1e-2);
            }
        }

        TEST_F(PolylineTest, VelCartCircle)
        {
            // sample radii and angles
            const Eigen::Array<double, Eigen::Dynamic, 1> radii { Eigen::Array<double, Eigen::Dynamic, 1>::Random(1).abs() * m_radius };
            const Eigen::Array<double, Eigen::Dynamic, 1> angles { Eigen::Array<double, Eigen::Dynamic, 1>::Random(1) * M_PI };

            // sample radial and angular velocities
            const Eigen::Array<double, Eigen::Dynamic, 1> radiiDt { Eigen::Array<double, Eigen::Dynamic, 1>::Random(1) * m_radius };
            const Eigen::Array<double, Eigen::Dynamic, 1> anglesDt { Eigen::Array<double, Eigen::Dynamic, 1>::Random(1) * M_PI / 8 };

            // determine cartesian positions
            const Points<Eigen::Dynamic, PointFrenet> posFrenet {m_radius * (angles + M_PI), m_radius - radii };

            // determine cartesian velocities
            const Points<Eigen::Dynamic, PointFrenet> input { m_radius * anglesDt, - radiiDt };

            // test frenet frame results since error grows with distance from path
            const auto result { m_circleTransform.velCartes(input, posFrenet) };

            const Points<Eigen::Dynamic, PointCartes> groundTruth { radiiDt * angles.cos() - radii * anglesDt * angles.sin(), radiiDt * angles.sin() + radii * anglesDt * angles.cos() };

            for(int index {}; index < posFrenet.numPoints(); ++index)
            {
                // relative error
                EXPECT_NEAR(result.x(index) / groundTruth.x(index), 1.0, 1e-3);
                EXPECT_NEAR(result.y(index) / groundTruth.y(index), 1.0, 6e-3);
                // absolute error
                EXPECT_NEAR(result.x(index), groundTruth.x(index), 3e-3);
                EXPECT_NEAR(result.y(index), groundTruth.y(index), 6e-3);
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
