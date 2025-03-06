#ifndef TEST_BASE_H
#define TEST_BASE_H

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>

namespace FrenetTransform
{
    namespace Testing
    {
        using limits = std::numeric_limits<double>;

        class TestBase : public testing::Test
        {
        protected:
            template <typename Tarray>
            void expectAllClose(const Eigen::ArrayBase<Tarray>& estimate, const Eigen::ArrayBase<Tarray>& groundTruth, double errAbs = limits::infinity(), double errRel = limits::infinity())
            {
                for(int row {}; row < estimate.rows(); ++row)
                {
                    for(int col {}; col < estimate.cols(); ++col)
                    {
                        EXPECT_NEAR(estimate(row, col), groundTruth(row, col), errAbs);
                        EXPECT_NEAR(estimate(row, col) / groundTruth(row, col), 1.0, errRel);
                    }
                }
            }
        };
    };
};

#endif