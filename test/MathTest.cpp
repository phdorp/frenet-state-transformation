#include "math.h"

#include <eigen3/Eigen/Core>
#include <gtest/gtest.h>


TEST(diff, RandomVector)
{
    Eigen::Array<double, 5, 1> vector {Eigen::Array<double, 5, 1>::Random()};
    auto result {FrenetTransform::diff(vector)};
    EXPECT_EQ(result.rows(), vector.rows());

    EXPECT_NEAR(result(0), 0.0, 0.000001);

    for(int cNumber {1}; cNumber < result.rows(); ++cNumber)
    {
        EXPECT_DOUBLE_EQ(result(cNumber), vector(cNumber) - vector(cNumber-1));
        ++cNumber;
    }
}

TEST(partialLength, StraightLine)
{
    Eigen::Array<double, 5, 1> x {0.0, 1.0, 2.0, 3.0, 4.0};
    Eigen::Array<double, 5, 1> y {0.0, 0.0, 0.0, 0.0, 0.0};

    Eigen::Array<double, 5, 1> result {FrenetTransform::partialLength(x, y)};
    EXPECT_EQ(result.rows(), x.rows());

    Eigen::Array<double, 5, 1> groundTruth {0.0, 1.0, 2.0, 3.0, 4.0};

    int cNumber {};
    for(auto number : result)
    {
        EXPECT_NEAR(number, groundTruth(cNumber), 0.000001);
        ++cNumber;
    }
}