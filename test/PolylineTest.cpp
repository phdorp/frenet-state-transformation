#include "polyline.h"
#include "points.h"
#include "transform.h"

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include <math.h>
#include <memory>

using FrenetTransform::Points;

// The fixture for testing class Foo.
class PolylineTest : public testing::Test
{
    protected:
    PolylineTest()
        : m_straight { {0, 1, 2, 4, 7},  {0, 1, 2, 4, 7} }
        , m_circle { m_radius * m_lengthsCircle.cos(), m_radius * m_lengthsCircle.sin() }
        , m_circleTransform { std::make_shared<FrenetTransform::Polyline<1079>>(m_circle) }
    {
    }

    const FrenetTransform::Polyline<5> m_straight;

    const double m_radius { 5.0 };
    const Eigen::Array<double, 1079, 1> m_lengthsCircle { Eigen::Array<double, 1079, 1>::LinSpaced(-M_PI, M_PI)} ;
    const FrenetTransform::Polyline<1079> m_circle {};
    const FrenetTransform::Transform m_circleTransform;
};

TEST_F(PolylineTest, GetPointsStraight)
{
    Eigen::ArrayXd input { {0.0, std::sqrt(8), std::sqrt(98)} };

    const auto result { m_straight(input) };

    const Points<Eigen::Dynamic> groundTruth { {{0.0, 2.0, 7.0}}, {{0.0, 2.0, 7.0}} };

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

    const auto result { m_circle(input) };

    const Points<Eigen::Dynamic> groundTruth
    {
        {{m_radius, std::sqrt(50),  0.0}},
        {{ 0.0, std::sqrt(50), m_radius}}
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

    const auto result { m_circle.tangent(input) };

    const Points<Eigen::Dynamic> groundTruth
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

    const auto result { m_circle.normal(input) };

    const Points<Eigen::Dynamic> groundTruth
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
    Eigen::ArrayXd input { {M_PI / 2, M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
    input *= m_radius;

    const auto result { m_circle.angle0(input) };

    const Eigen::ArrayXd groundTruth {{0.0, M_PI / 2, 3 * M_PI / 4, 0.0}};

    for(int index {}; index < input.rows(); ++index)
    {
        EXPECT_NEAR(groundTruth(index), result(index), 1e-2);
    }
}

TEST_F(PolylineTest, GetAngles1Circle)
{
    Eigen::ArrayXd input { {M_PI / 2, M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
    input *= m_radius;

    const auto result { m_circle.angle1(input) };

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

    const auto result { m_circle.angle2(input) };

    const Eigen::ArrayXd groundTruth { Eigen::ArrayXd::Zero(input.size()) };

    for(int index {}; index < input.rows(); ++index)
    {
        EXPECT_NEAR(groundTruth(index), result(index), 1e-2);
    }
}

TEST_F(PolylineTest, NextPointsStraight)
{
    Points<Eigen::Dynamic> input {
        {{0.0, 3.0, 2.0}},
        {{0.0, 3.0, 4.0}}
    };

    const auto lengths { m_straight.lengths(input) };
    const auto result { m_straight(lengths) };

    const Points<Eigen::Dynamic> groundTruth {
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
    Points<Eigen::Dynamic> input {
        {{  0.0,  0.0, 1.0}},
        {{-m_radius, -2.0, 1.0}}
    };

    const auto lengths { m_circle.lengths(input) };
    const auto result { m_circle(lengths) };

    const Points<Eigen::Dynamic> groundTruth {
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
    const Points<Eigen::Dynamic> input {
        {{      0.0,  0.0, 1.0}},
        {{-m_radius, -2.0, 1.0}}
    };

    const auto result { m_circleTransform.posFrenet(input) };

    const Points<Eigen::Dynamic> groundTruth {
        {{ M_PI / 2 * m_radius, M_PI / 2 * m_radius, 5 * M_PI / 4 * m_radius }},
        {{                 0.0,      m_radius - 2.0, m_radius - std::sqrt(2) }}
    };

    for(int index {}; index < input.numPoints(); ++index)
    {
        EXPECT_NEAR(groundTruth.x(index), result.x(index), 1e-2);
        EXPECT_NEAR(groundTruth.y(index), result.y(index), 1e-2);
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
