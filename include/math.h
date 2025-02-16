#ifndef MATH_H
#define MATH_H

#include <eigen3/Eigen/Core>

namespace FrenetTransform
{
    template <int T>
    using ArrayT1 = Eigen::Array<double, T, 1>;

    template <typename T>
    T diff(const Eigen::ArrayBase<T>& numbers)
    {
        T result {};
        result(Eigen::seq(1,result.rows()-1)) = numbers(Eigen::seq(1,numbers.rows()-1)) - numbers(Eigen::seq(0,numbers.rows()-2));
        return result;
    }

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