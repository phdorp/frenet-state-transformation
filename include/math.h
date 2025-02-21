#ifndef MATH_H
#define MATH_H

#include <eigen3/Eigen/Core>

namespace FrenetTransform
{
    template <int T>
    using ArrayT1 = Eigen::Array<double, T, 1>;

    /**
     * @brief Perform row-wise backward differences on the input array.
     *
     * Due to backward differences, the first row of the return array is zero.
     *
     * @tparam T array type.
     * @param numbers input array.
     * @return T array of row-wise differences.
     */
    template <typename T>
    T diffBackward(const Eigen::ArrayBase<T>& numbers)
    {
        const auto colsAll {Eigen::seq(0, numbers.cols() - 1)}; // sequence over all column indices
        const auto rowsExLast {Eigen::seq(0, numbers.rows()-2)}; // sequence over all row indices except last one
        const auto rowsExFirst {Eigen::seq(1, numbers.rows()-1)}; // sequence over all row indices except first one

        T result {}; // result array

        // perform differences
        result(rowsExFirst, colsAll) = numbers(rowsExFirst, colsAll) - numbers(rowsExLast, colsAll);

        return result;
    }

    /**
     * @brief Perform row-wise forward differences on the input array.
     *
     * Due to forward differences, the first row of the return array is zero.
     *
     * @tparam T array type.
     * @param numbers input array.
     * @return T array of row-wise differences.
     */
    template <typename T>
    T diffForward(const Eigen::ArrayBase<T>& numbers)
    {
        T result { diffBackward(numbers) }; // get backward differences

        const auto colsAll {Eigen::seq(0, numbers.cols() - 1)}; // sequence over all column indices
        const auto rowsExLast {Eigen::seq(0, numbers.rows() - 2)}; // sequence over all row indices except last one
        const auto rowsExFirst {Eigen::seq(1, numbers.rows() - 1)}; // sequence over all row indices except first one

        // shift numbers by one element to vector beginning
        result(rowsExLast, colsAll) = result(rowsExFirst, colsAll);
        result(numbers.rows() - 1, colsAll) = 0.0;

        return result;
    }

    /**
     * @brief Provide cumulative lengths input given points.
     *
     * @tparam T row of input vector.
     * @param x coordinates in x-direction.
     * @param y coordinates in y-direction.
     * @return ArrayT1<T> cumulative lengths input given points.
     */
    template <int T>
    ArrayT1<T> partialLength(const ArrayT1<T>& x, const ArrayT1<T>& y)
    {
        ArrayT1<T> result {FrenetTransform::diffBackward(x).pow(2) + FrenetTransform::diffBackward(y).pow(2)}; // sum of squared distances
        result = result.sqrt(); // determine root of sum of squared distances

        // accumulate distances
        for(int cNumber {}; cNumber < result.rows() - 1; ++cNumber)
            result(cNumber+1) += result(cNumber);

        return result;
    }

    /**
     * @brief Provide index of minimum positive element in an ordered sequence.
     *
     * @tparam T rows of input vector.
     * @param sequence input to search.
     * @return int index of minimum positive element. -1 if no positive element.
     */
    template <typename T>
    int first(const Eigen::ArrayBase<T>& sequence)
    {
        for(int index {}; index < sequence.rows(); ++index)
        {
            if(sequence(index) >=0)
                return index;
        }

        return -1;
    }

    // template <typename T>
    // T
};

#endif