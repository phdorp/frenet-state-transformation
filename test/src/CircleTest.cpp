#include "CircleTest.h"

namespace FrenetTransform
{
    namespace Testing
    {
        TEST_F(CircleTest, asdf)
        {
            const auto anglesEst { m_circle.angle0(m_posFrenet.x()) };

            Eigen::ArrayXd anglesGtr { m_posCircle.y() + M_PI / 2};
            // map angle in -PI to PI
            std::for_each(anglesGtr.begin(), anglesGtr.end(),
                [](double& angle) { angle = std::fmod(angle + M_PI, 2 * M_PI) - M_PI; });

            expectAllClose(anglesEst, anglesGtr, 1e-2);
        }
    };
};

int main(int argc, char **argv)
{
    std::srand(0);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}