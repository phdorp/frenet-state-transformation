#include "frenetTransform/test/pathCircleTest.h"
#include "frenetTransform/internal/constexprTypes.h"

namespace FrenetTransform
{
    namespace Internal
    {
        // using Queries = testing::Types<Integral<Eigen::Dynamic>, Integral<100>>;
        using Queries = testing::Types<Integral<Eigen::Dynamic>, Integral<100>>;
        TYPED_TEST_SUITE(PathCircleTest, Queries);

        TYPED_TEST(PathCircleTest, DirectionAngle)
        {
            const auto anglesEst { this->m_circle.angle0(this->m_posFrenet.x()) };

            Eigen::Array<double, TypeParam::s_val, 1> anglesGtr { this->m_posCircle.y() + M_PI / 2};
            // map angle in -PI to PI
            std::for_each(anglesGtr.begin(), anglesGtr.end(),
                [](double& angle) { angle = std::fmod(angle + M_PI, 2 * M_PI) - M_PI; });

            this->expectAllClose(anglesEst, anglesGtr, 1e-2);
        }
    };
};

int main(int argc, char **argv)
{
    std::srand(0);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}