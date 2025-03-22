#include "frenetTransform/test/pathPolylineLineTest.h"

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

namespace FrenetTransform
{
    namespace Testing
    {
        TEST_F(PathPolylineLineTest, Lengths)
        {
            const auto lengthsGtr { m_line.lengths(m_pointsCartes) };
            const auto lengthsEst { m_polyline.lengths(m_pointsCartes) };

            expectAllClose(lengthsGtr, lengthsEst, 1e-10);
        }
    };
};

int main(int argc, char **argv)
{
    std::srand(0);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
