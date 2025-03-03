#include "CircleTest.h"

namespace Testing
{
    TEST_F(CircleTest, asdf)
    {
        const auto result { m_circle.angle0(m_circle.lengths(m_pointsCircle.y())) };

        Eigen::ArrayXd groundTruth { m_pointsCircle.y() - M_PI / 2};
        for(int index {}; index < groundTruth.rows(); ++index)
            if(groundTruth(index) > M_PI)
                groundTruth(index) += -2*M_PI;

        for(int index {}; index < groundTruth.rows(); ++index)
            EXPECT_NEAR(groundTruth(index), result(index), 1e-2);
    }
};

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}