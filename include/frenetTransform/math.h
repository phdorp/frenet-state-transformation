#ifndef MATH_H
#define MATH_H

#include <Eigen/Core>

namespace FrenetTransform
{
    template <int NumRows>
    using ArrayRows = Eigen::Array<double, NumRows, 1>;

    /**
     * @brief Perform row-wise backward differences on the input array.
     *
     * Due to backward differences, the first row of the return array is zero.
     *
     * @tparam ArrayType array type.
     * @param numbers input array.
     * @return ArrayType array of row-wise differences.
     */
    template <typename ArrayType>
    ArrayType diffBackward(const Eigen::ArrayBase<ArrayType>& numbers)
    {
        const auto colsAll {Eigen::seq(0, numbers.cols() - 1)}; // sequence over all column indices
        const auto rowsExLast {Eigen::seq(0, numbers.rows()-2)}; // sequence over all row indices except last one
        const auto rowsExFirst {Eigen::seq(1, numbers.rows()-1)}; // sequence over all row indices except first one

        ArrayType result { Eigen::ArrayBase<ArrayType>::Zero(numbers.rows(), numbers.cols()) }; // result array

        // perform differences
        result(rowsExFirst, colsAll) = numbers(rowsExFirst, colsAll) - numbers(rowsExLast, colsAll);

        return result;
    }

    /**
     * @brief Perform row-wise forward differences on the input array.
     *
     * Somehow induces wrong results on repeated use in TEST_F preparation.
     *
     * Due to forward differences, the first row of the return array is zero.
     *
     * @tparam ArrayType array type.
     * @param numbers input array.
     * @return ArrayType array of row-wise differences.
     */
    template <typename ArrayType>
    ArrayType diffForward(const Eigen::ArrayBase<ArrayType>& numbers)
    {
        ArrayType result { diffBackward(numbers) }; // get backward differences

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
     * @tparam NumPoints row of input vector.
     * @param x coordinates in x-direction.
     * @param y coordinates in y-direction.
     * @return ArrayRows<NumPoints> cumulative lengths input given points.
     */
    template <int NumPoints>
    ArrayRows<NumPoints> partialLength(const ArrayRows<NumPoints>& x, const ArrayRows<NumPoints>& y)
    {
        ArrayRows<NumPoints> result {FrenetTransform::diffBackward(x).pow(2) + FrenetTransform::diffBackward(y).pow(2)}; // sum of squared distances
        result = result.sqrt(); // determine root of sum of squared distances

        // accumulate distances
        for(int cNumber {}; cNumber < result.rows() - 1; ++cNumber)
            result(cNumber+1) += result(cNumber);

        return result;
    }

    /**
     * @brief Provide index of minimum positive element in an ordered sequence.
     *
     * @tparam ArrayType rows of input vector.
     * @param sequence input to search.
     * @return int index of minimum positive element. -1 if no positive element.
     */
    template <typename ArrayType>
    int first(const Eigen::ArrayBase<ArrayType>& sequence)
    {
        int index {};

        for(auto val : sequence)
        {
            if(val >=0)
                return index;
            ++index;
        }

        return index - 1;
    }

    /**
     * @brief Determine gradient from finite differences.
     *
     * @tparam NumRows number of samples.
     * @tparam NumCols number of dependent variables.
     * @param depents dependent variables.
     * @param indepents independent variables.
     * @return Eigen::Array<double, NumRows, NumCols> finite difference gradient.
     */
    template <int NumRows, int NumCols>
    Eigen::Array<double, NumRows, NumCols> gradient(const Eigen::Array<double, NumRows, NumCols>& depents, const Eigen::Array<double, NumRows, 1>& indepents)
    {
        const Eigen::Array<double, NumRows, NumCols> diffDepents { diffBackward(depents) }; // finite differences dependent variables
        const Eigen::Array<double, NumRows, 1> diffIndepents { diffBackward(indepents) }; // finite differences independent variables

        Eigen::Array<double, NumRows, NumCols> result { Eigen::Array<double, NumRows, NumCols>::Zero(depents.rows(), depents.cols()) }; // instantiate result array

        // perform column-wise normalization of dependent differences
        for(int col {}; col < NumCols; ++col)
            result(Eigen::all, col) = diffDepents(Eigen::all, col) / diffIndepents;

        // results from 0 / 0 are set to 0
        for(int row {}; row < NumRows; ++row)
            if(std::abs(diffIndepents(row)) < 1e-10)
                for(int col {}; col < NumCols; ++col)
                    result(row, col) = std::abs(diffDepents(row, col)) < 1e-10 ? 0 : result(row, col);

        return result;
    }

    template <typename ArrayType>
    ArrayType angleDir(const Eigen::ArrayBase<ArrayType>& dirxs, const Eigen::ArrayBase<ArrayType>& dirys)
    {
        ArrayType angles { Eigen::ArrayBase<ArrayType>::Zero(dirxs.rows(), dirxs.cols()) };

        for(int row {}; row < dirxs.rows(); ++row)
        {
            for(int col{}; col < dirxs.cols(); ++col)
            {
                if(dirxs(row, col) > 0 && dirxs(row, col) > std::abs(dirys(row, col)))
                    angles(row, col) = std::atan(dirys(row, col) / dirxs(row, col));
                else if(dirys(row, col) > 0 && std::abs(dirxs(row, col)) < dirys(row, col))
                    angles(row, col) = M_PI / 2 - std::atan(dirxs(row, col) / dirys(row, col));
                else if(dirxs(row, col) < 0 && dirys(row, col) > 0 && -dirxs(row, col) > dirys(row, col))
                    angles(row, col) = M_PI + std::atan(dirys(row, col) / dirxs(row, col));
                else if(dirys(row, col) < 0 && std::abs(dirxs(row, col)) < -dirys(row, col))
                    angles(row, col) = -M_PI / 2 - std::atan(dirxs(row, col) / dirys(row, col));
                else
                    angles(row, col) = -M_PI + std::atan(dirys(row, col) / dirxs(row, col));
            }
        }

        return angles;
    }
};

#endif