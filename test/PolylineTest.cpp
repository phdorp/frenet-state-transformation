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
    const Eigen::Array<double, 200, 1> m_lengthsCircle { Eigen::Array<double, 200, 1>::LinSpaced(0.0, M_PI)} ;
    const FrenetTransform::Polyline<200> m_circle {};
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
    Eigen::ArrayXd input { {0.0, M_PI/8, M_PI/4} };
    input *= m_radius * 2;

    const auto result { m_circle(input) };

    const FrenetTransform::Path::Points groundTruth
    {
        {{10.0, std::sqrt(50),  0.0}},
        {{ 0.0, std::sqrt(50), 10.0}}
    };

    for(int index {}; index < input.rows(); ++index)
    {
        EXPECT_NEAR(groundTruth.x(index), result.x(index), 1e-3);
        EXPECT_NEAR(groundTruth.y(index), result.y(index), 1e-3);
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
