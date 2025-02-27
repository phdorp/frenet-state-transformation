#ifndef POLYLINE_H
#define POLYLINE_H

#include <math.h>

#include "path.h"
#include "math.h"
#include "points.h"

using FrenetTransform::Points;

namespace FrenetTransform
{
    /**
     * @brief Path representation as polyline.
     *
     * Represents a 2-dimensional path as a polyline.
     * Provide path properties based on finite differences at query points.
     *
     * @tparam T determines the number points along the path.
     */
    template <int T>
    class Polyline : public Path
    {
    public:
        using MatrixT1 = Eigen::Matrix<double, T, 1>;
        using MatrixT2 = Eigen::Matrix<double, T, 2>;

        using ArrayT1 = Eigen::Array<double, T, 1>;

        Polyline() = default;

        /**
         * @brief Construct a new Path object from Cartesian x- and y-positions.
         *
         * @param x coordinates in x-direction.
         * @param y coordinates in y-direction.
         */
        Polyline(const ArrayT1& x, const ArrayT1& y) { setPoints(x, y); }

        Points<Eigen::Dynamic> operator()(const Eigen::ArrayXd& lengths) const override
        {
            // indices of corresponding polyline segments
            const Eigen::ArrayXi indicesLengths { indices(lengths) };

            // relative position along the linear segment
            const auto pathLengthsdiff { FrenetTransform::diffBackward(m_lengths) };
            const auto relativePos { (lengths - m_lengths(indicesLengths)) / pathLengthsdiff(indicesLengths + 1) };

            // absolute position along path
            const auto x { m_x[0](indicesLengths + 1) * relativePos + m_x[0](indicesLengths) * (1 - relativePos) };
            const auto y { m_y[0](indicesLengths + 1) * relativePos + m_y[0](indicesLengths) * (1 - relativePos) };

            return { x, y };
        }

        /**
         * @brief Determines next points to the query points.
         *
         * @param points query points.
         * @return Points<Eigen::Dynamic> next to query points.
         */
        Points<Eigen::Dynamic> nextPoints(const Points<Eigen::Dynamic>& points) const override
        {
            Eigen::ArrayXd xNext (points.numPoints());
            Eigen::ArrayXd yNext (points.numPoints());

            for(int cPoints {}; cPoints < points.numPoints(); ++cPoints)
            {
                double distSq { -1.0 };

                for(int cThis {}; cThis < T; ++cThis)
                {
                    const double xDiffPnt { m_x[0].data()[cThis] - points.x(cPoints) };
                    const double yDiffPnt { m_y[0].data()[cThis] - points.y(cPoints) };
                    const double param { (xDiffPnt * m_xDiff(cThis) + yDiffPnt * m_yDiff(cThis))
                        / (std::pow(m_xDiff(cThis), 2) + std::pow(m_yDiff(cThis), 2)) };

                    if(param >= 0.0 && param <= 1.0)
                    {
                        const double xCand { m_x[0].data()[cThis - 1] * param + (1 - param) * m_x[0].data()[cThis] };
                        const double yCand { m_y[0].data()[cThis - 1] * param + (1 - param) * m_y[0].data()[cThis] };
                        const double distSqCand { std::pow(xCand - points.x(cPoints), 2) + std::pow(yCand - points.y(cPoints), 2) };

                        if(distSq < 0.0 || distSq > distSqCand)
                        {
                            xNext(cPoints) = xCand;
                            yNext(cPoints) = yCand;
                            distSq = distSqCand;
                        }
                    }
                }
            }

            return { xNext, yNext };
        }

        void setPoints(const ArrayT1& x, const ArrayT1& y)
        {
            m_x[0] = x;
            m_xDiff = FrenetTransform::diffBackward(x);

            m_y[0] = y;
            m_yDiff = FrenetTransform::diffBackward(y);

            m_lengths = FrenetTransform::partialLength(x, y);

            for(unsigned int orderGrad { 1 }; orderGrad < s_numGrad; ++orderGrad)
            {
                m_x[orderGrad] = FrenetTransform::gradient(m_x[orderGrad - 1], m_lengths);
                m_y[orderGrad] = FrenetTransform::gradient(m_y[orderGrad - 1], m_lengths);
            }
        }

    private:
        static constexpr int s_numGrad { 4 };
        ArrayT1 m_xDiff {};
        ArrayT1 m_yDiff {};
        std::array<ArrayT1, s_numGrad> m_x {}; /*<< coordinates and gradients in x-direction*/
        std::array<ArrayT1, s_numGrad> m_y {}; /*<< coordinates and gradients in y-direction*/
        ArrayT1 m_lengths {}; /*<< partial lengths along polyline*/

        /**
         * @brief Determines 1st order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 1st order gradient at given path lengths.
         */
        Points<Eigen::Dynamic> gradient1 (const Eigen::MatrixXd& lengths) const override
        {
            const auto indicesGrad { indices(lengths) };
            return { m_x[1](indicesGrad), m_y[1](indicesGrad) };
        }

        /**
         * @brief Determines 2nd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 2nd order gradient at given path lengths.
         */
        Points<Eigen::Dynamic> gradient2 (const Eigen::MatrixXd& lengths) const override
        {
            const auto indicesGrad { indices(lengths) };
            return { m_x[2](indicesGrad), m_y[2](indicesGrad) };
        }

        /**
         * @brief Determines 3rd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 3rd order gradient at given path lengths.
         */
        Points<Eigen::Dynamic> gradient3 (const Eigen::MatrixXd& lengths) const override
        {
            const auto indicesGrad { indices(lengths) };
            return { m_x[3](indicesGrad), m_y[3](indicesGrad) };
        }

        /**
         * @brief Determines indices of polyline segment corresponding to the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return indices corresponding to given path lengths.
         */
        Eigen::ArrayXi indices(const Eigen::ArrayXd& lengths) const
        {
            Eigen::ArrayXi result(lengths.size()); // vector of segment indices

            // get indices of next segments
            for(int row {}; row < lengths.rows(); ++row)
                result(row) = FrenetTransform::first(m_lengths - lengths(row));

            return result;
        }
    };
};

#endif