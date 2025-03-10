#include "math.h"

#include <eigen3/Eigen/Core>
#include <gtest/gtest.h>

namespace FrenetTransform
{
    namespace Testing
    {
        /**
         * @brief Test function diffBackward with vector of length 5.
         *
         */
        TEST(diffBackward, Vector5)
        {
            Eigen::Array<double, 5, 1> input {1, 2, 5, 1, 3}; // input vector

            auto result {FrenetTransform::diffBackward(input)}; // function result

            EXPECT_EQ(result.rows(), input.rows()); // test equality of input and result rows

            Eigen::Array<double, 5, 1> groundTruth {0, 1, 3, -4, 2}; // ground truth result

            // test for equality between ground truth and result
            for(int index {}; index < result.rows(); ++index)
                EXPECT_NEAR(result(index), groundTruth(index), 1.0e-10);
        }

        TEST(diffBackward, Matrix52)
        {
            // input matrix
            Eigen::Array<double, 5, 2> input {};
            input.col(0) << 1, 2, 5, 1, 3;
            input.col(1) << 2, 2, 5, 1, 3;

            auto result {FrenetTransform::diffBackward(input)}; // function result

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

        TEST(diffForward, Matrix52)
        {
            // input matrix
            Eigen::Array<double, 5, 2> input {};
            input.col(0) << 1, 2, 5, 1, 3;
            input.col(1) << 2, 2, 5, 1, 3;

            auto result {FrenetTransform::diffForward(input)}; // function result

            EXPECT_EQ(result.rows(), input.rows()); // test equality of input and result rows
            EXPECT_EQ(result.cols(), input.cols()); // test equality of input and result columns

            // ground truth result
            Eigen::Array<double, 5, 2> groundTruth {};
            groundTruth.col(0) << 1, 3, -4, 2, 0;
            groundTruth.col(1) << 0, 3, -4, 2, 0;

            // test for equality between ground truth and result
            for(int row {}; row < result.rows(); ++row)
                for(int col {}; col < result.cols(); ++col)
                    EXPECT_NEAR(result(row, col), groundTruth(row, col), 1.0e-10);
        }

        TEST(gradient, Matrix52)
        {
            // dependent matrix
            Eigen::Array<double, 5, 2> depents {};
            depents.col(0) << 1.0, 2.0, 5.0, 1.0, 3.0;
            depents.col(1) << 2.0, 2.0, 5.0, 1.0, 3.0;

            // independent vector
            Eigen::Array<double, 5, 1> indepents { 1.0, 2.0, 3.0, 5.0, 0.0};

            auto result {FrenetTransform::gradient(depents, indepents)}; // function result

            EXPECT_EQ(result.rows(), depents.rows()); // test equality of input and result rows
            EXPECT_EQ(result.cols(), depents.cols()); // test equality of input and result columns

            // ground truth result
            Eigen::Array<double, 5, 2> groundTruth {};
            groundTruth.col(0) << 0.0, 1.0, 3.0, -2.0, 2.0 / -5.0;
            groundTruth.col(1) << 0.0, 0.0, 3.0, -2.0, 2.0 / -5.0;

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
            Eigen::Array<double, 5, 1> x {0.0, 1.0, 2.0, 3.0, 4.0}; // input in x-direction
            Eigen::Array<double, 5, 1> y {0.0, 0.0, 0.0, 0.0, 0.0}; // input in y-direction

            Eigen::Array<double, 5, 1> result {FrenetTransform::partialLength(x, y)}; // partial length result

            EXPECT_EQ(result.rows(), x.rows()); // test equality of input and result rows

            Eigen::Array<double, 5, 1> groundTruth {0.0, 1.0, 2.0, 3.0, 4.0}; // ground truth result

            // test for equality between ground truth and result
            for(int row {}; row < result.rows(); ++row)
                EXPECT_NEAR(result(row), groundTruth(row), 1e-10);
        }

        /**
         * @brief Test correct index from increasing sequence with positive elements.
         *
         */
        TEST(first, IncreasingSequencePositive)
        {
            Eigen::Array<double, 5, 1> input { -3, -1, 3, 5, 6 }; // increasing input sequence

            int result { FrenetTransform::first(input) }; // get first index

            int groundTruth { 2 }; // ground truth index

            EXPECT_EQ(result, groundTruth);
        }

        /**
         * @brief Test negative index from sequence with no positive element.
         *
         */
        TEST(first, IncreasingSequenceNegative)
        {
            Eigen::Array<double, 5, 1> input { -3, -1, -2, -6, -2 }; // negative input sequence

            int result { FrenetTransform::first(input) }; // get negative index

            int groundTruth { 4 }; // ground truth index

            EXPECT_EQ(result, groundTruth);
        }
    };
};