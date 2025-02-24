#include "polyline.h"

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include <math.h>

// The fixture for testing class Foo.
class PolylineTest : public testing::Test
{
    protected:
    PolylineTest()
        : m_straight { {0, 1, 2, 4, 7},  {0, 1, 2, 4, 7} }
        , m_circle { m_radius * m_lengthsCircle.cos(), m_radius * m_lengthsCircle.sin() }
    {
    }

    const FrenetTransform::Polyline<5> m_straight;

    const double m_radius { 10.0 };
    const Eigen::Array<double, 400, 1> m_lengthsCircle { Eigen::Array<double, 400, 1>::LinSpaced(-M_PI, M_PI)} ;
    const FrenetTransform::Polyline<400> m_circle {};
};

TEST_F(PolylineTest, GetPointsStraight)
{
    Eigen::ArrayXd input { {0.0, std::sqrt(8), std::sqrt(98)} };

    const auto result { m_straight(input) };

    const FrenetTransform::Path::Points groundTruth { {{0.0, 2.0, 7.0}}, {{0.0, 2.0, 7.0}} };

    for(int index {}; index < input.rows(); ++index)
    {
        EXPECT_NEAR(groundTruth.x(index), result.x(index), 1e-10);
        EXPECT_NEAR(groundTruth.y(index), result.y(index), 1e-10);
    }
}

TEST_F(PolylineTest, GetPointsCircle)
{
    Eigen::ArrayXd input { {M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
    input *= m_radius;

    const auto result { m_circle(input) };

    const FrenetTransform::Path::Points groundTruth
    {
        {{10.0, std::sqrt(50),  0.0}},
        {{ 0.0, std::sqrt(50), 10.0}}
    };

    for(int index {}; index < input.rows(); ++index)
    {
        EXPECT_NEAR(groundTruth.x(index), result.x(index), 1e-2);
        EXPECT_NEAR(groundTruth.y(index), result.y(index), 1e-2);
    }
}

TEST_F(PolylineTest, GetTangentsCircle)
{
    Eigen::ArrayXd input { {M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
    input *= m_radius;

    const auto result { m_circle.tangent(input) };

    const FrenetTransform::Path::Points groundTruth
    {
        {{0.0, -std::sqrt(50) / 10, -1.0}},
        {{1.0,  std::sqrt(50) / 10,  0.0}}
    };

    for(int index {}; index < input.rows(); ++index)
    {
        EXPECT_NEAR(groundTruth.x(index), result.x(index), 1e-2);
        EXPECT_NEAR(groundTruth.y(index), result.y(index), 1e-2);
    }
}

TEST_F(PolylineTest, GetNormalsCircle)
{
    Eigen::ArrayXd input { {M_PI, M_PI + M_PI / 4, M_PI + M_PI / 2} };
    input *= m_radius;

    const auto result { m_circle.normal(input) };

    const FrenetTransform::Path::Points groundTruth
    {
        {{1.0, std::sqrt(50) / 10, 0.0}},
        {{0.0, std::sqrt(50) / 10, 1.0}}
    };

    for(int index {}; index < input.rows(); ++index)
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
