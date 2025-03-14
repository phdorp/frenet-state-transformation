#ifndef TEST_BASE_H
#define TEST_BASE_H

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <array>

namespace FrenetTransform
{
    namespace Testing
    {
        template <int Value>
        struct Integral { static constexpr int s_val { Value }; };

        template <typename ValType, int NumVals, std::array<ValType, NumVals> Vals>
        struct ConstVals { static constexpr std::array<ValType, NumVals> s_vals { Vals }; };

        using limits = std::numeric_limits<double>;

        class TestBase : public testing::Test
        {
        protected:
            template <typename Tarray>
            void expectAllClose(const Eigen::ArrayBase<Tarray>& estimate, const Eigen::ArrayBase<Tarray>& groundTruth, double errAbs = limits::infinity(), double errRel = limits::infinity())
            {
                const auto errAbsVals { (estimate - groundTruth).abs() };
                const auto errRelVals { errAbsVals / groundTruth };

                for(int idx {}; idx < estimate.size(); ++idx)
                {
                    EXPECT_LT(errAbsVals(idx), errAbs) << "Error bound violated at index " << idx << '!';
                    if(!std::isinf(errRelVals(idx)))
                        EXPECT_LT(errRelVals(idx), errRel) << "Error bound violated at index " << idx << '!';
                }
            }
        };
    };
};

#endif