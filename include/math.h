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
    T diff(const Eigen::ArrayBase<T>& numbers)
    {
        T result {};
        result(Eigen::seq(1,result.rows()-1)) = numbers(Eigen::seq(1,numbers.rows()-1)) - numbers(Eigen::seq(0,numbers.rows()-2));
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
        ArrayT1<T> length { };
        length = FrenetTransform::diff(x).pow(2) + FrenetTransform::diff(y).pow(2);
        length = length.sqrt();

        for(int cNumber {}; cNumber < length.rows() - 1; ++cNumber)
            length(cNumber+1) += length(cNumber);

        return length;
    }
};

#endif