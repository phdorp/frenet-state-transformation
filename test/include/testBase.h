#ifndef TEST_BASE_H
#define TEST_BASE_H

#include <gtest/gtest.h>
#include <Eigen/Core>

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
                const auto errAbsVals { (estimate - groundTruth).abs() };
                EXPECT_LT(*std::max_element(errAbsVals.begin(), errAbsVals.end()), errAbs);
            }
        };
    };
};

#endif