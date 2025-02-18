#include "math.h"

#include <eigen3/Eigen/Core>
#include <gtest/gtest.h>

/**
 * @brief Test function diff with vector of length 5.
 *
 */
TEST(diff, Vector5)
{
    Eigen::Array<double, 5, 1> input {1, 2, 5, 1, 3}; // input vector

    auto result {FrenetTransform::diff(input)}; // function result

    EXPECT_EQ(result.rows(), input.rows()); // test equality of input and result rows

    Eigen::Array<double, 5, 1> groundTruth {0, 1, 3, -4, 2}; // ground truth result

    // test for equality between ground truth and result
    for(int index {}; index < result.rows(); ++index)
        EXPECT_NEAR(result(index), groundTruth(index), 1.0e-10);
}

TEST(diff, Matrix52)
{
    // input matrix
    Eigen::Array<double, 5, 2> input {};
    input.col(0) << 1, 2, 5, 1, 3;
    input.col(1) << 2, 2, 5, 1, 3;

    auto result {FrenetTransform::diff(input)}; // function result

    EXPECT_EQ(result.rows(), input.rows()); // test equality of input and result rows
    EXPECT_EQ(result.cols(), input.cols()); // test equality of input and result columns

    // ground truth result
    Eigen::Array<double, 5, 2> groundTruth {};
    groundTruth.col(0) << 0, 1, 3, -4, 2;
    groundTruth.col(1) << 0, 0, 3, -4, 2;

    // test for equality between ground truth and result
    for(int row {}; row < result.rows(); ++row)
        for(int col {}; col < result.cols(); ++col)
            EXPECT_NEAR(result(row, col), groundTruth(row, col), 1.0e-10);
}

/**
 * @brief Test partial length with straigth line.
 *
 */
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