#ifndef POLYLINE_H
#define POLYLINE_H

#include <math.h>

#include "frenetTransform/path.h"
#include "frenetTransform/math.h"
#include "frenetTransform/points.h"
#include "frenetTransform/point.h"

namespace FrenetTransform
{
    /**
     * @brief Path representation as polyline.
     *
     * Represents a 2-dimensional path as a polyline.
     * Provide path properties based on finite differences at query points.
     *
     * @tparam NumPoints determines the number points along the path.
     */
    template <int NumPoints>
    class Polyline : public Path
    {
    public:
        using ArrayPoints = Eigen::Array<double, NumPoints, 1>;

        Polyline() = default;

        /**
         * @brief Construct a new Path object from Cartesian x- and y-positions.
         *
         * @param x coordinates in x-direction.
         * @param y coordinates in y-direction.
         */
        Polyline(const ArrayPoints& x, const ArrayPoints& y) { setPoints(x, y); }

        Polyline(const Points<NumPoints>& points) { setPoints(points.x(), points.y()); }

        Points<Eigen::Dynamic> operator()(const Eigen::ArrayXd& lengths) const override
        {
            // indices of corresponding polyline segments
            const Eigen::ArrayXi indicesLengths { indices(lengths) };

            // relative position along the linear segment
            const Eigen::ArrayXd pathLengthsdiff { FrenetTransform::diffBackward(m_lengths) };
            Eigen::ArrayXd relativePos { (lengths - m_lengths(indicesLengths)) / pathLengthsdiff(indicesLengths + 1) };

            for(double& pos : relativePos)
                pos = std::clamp(pos, 0.0, 1.0);

            // absolute position along path
            Eigen::ArrayXd x { m_x[0](indicesLengths + 1) * relativePos + m_x[0](indicesLengths) * (1 - relativePos) };
            Eigen::ArrayXd y { m_y[0](indicesLengths + 1) * relativePos + m_y[0](indicesLengths) * (1 - relativePos) };

            for(int row {}; row < indicesLengths.rows(); ++row)
            {
                if(indicesLengths(row) >= m_lengths.rows() - 1)
                {
                    x(row) = m_x[0](m_lengths.rows() - 1);
                    y(row) = m_y[0](m_lengths.rows() - 1);
                }
            }

            return { x, y };
        }

        /**
         * @brief Determines next points to the query points.
         *
         * @param points query points.
         * @return Points<Eigen::Dynamic> next to query points.
         */
        Eigen::ArrayXd lengths(const Points<Eigen::Dynamic>& points) const override
        {
            Eigen::ArrayXd lLenghts (points.numPoints());

            for(int cPoints {}; cPoints < points.numPoints(); ++cPoints)
            {
                double distSq { -1.0 };
                for(int cThis {1}; cThis < m_lengths.rows(); ++cThis)
                {
                    const Point pathNext { m_x[0].data()[cThis], m_y[0].data()[cThis] };
                    const double xDiffPnt { pathNext.x() - points.x(cPoints) };
                    const double yDiffPnt { pathNext.y() - points.y(cPoints) };
                    const double param { (xDiffPnt * m_xDiff(cThis) + yDiffPnt * m_yDiff(cThis))
                        / (std::pow(m_xDiff(cThis), 2) + std::pow(m_yDiff(cThis), 2)) };

                    double distSqCand {};
                    double candLength {};
                    if(param >= 1.0)
                    {
                        distSqCand = std::pow( m_x[0].data()[cThis - 1] - points.x(cPoints), 2) + std::pow(m_y[0].data()[cThis - 1] - points.y(cPoints), 2);
                        candLength = m_lengths.data()[cThis - 1];
                    }
                    else if(param <= 0.0)
                    {
                        distSqCand = std::pow( m_x[0].data()[cThis] - points.x(cPoints), 2) + std::pow(m_y[0].data()[cThis] - points.y(cPoints), 2);
                        candLength = m_lengths.data()[cThis];
                    }
                    else
                    {
                        const Point pathPrev { m_x[0].data()[cThis - 1], m_y[0].data()[cThis - 1] };
                        const double xCand { pathPrev.x() * param + (1 - param) * pathNext.x() };
                        const double yCand { pathPrev.y() * param + (1 - param) * pathNext.y() };
                        distSqCand = std::pow(xCand - points.x(cPoints), 2) + std::pow(yCand - points.y(cPoints), 2);
                        candLength = m_lengths(cThis - 1) +  (Point { xCand, yCand }).distance(pathPrev);
                    }

                    if(distSqCand < distSq || distSq < 0)
                    {
                        distSq = distSqCand;
                        lLenghts(cPoints) = candLength;
                    }
                }
            }

            return lLenghts;
        }

        void setPoints(const ArrayPoints& x, const ArrayPoints& y)
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
        ArrayPoints m_xDiff {};
        ArrayPoints m_yDiff {};
        std::array<ArrayPoints, s_numGrad> m_x {}; /*<< coordinates and gradients in x-direction*/
        std::array<ArrayPoints, s_numGrad> m_y {}; /*<< coordinates and gradients in y-direction*/
        ArrayPoints m_lengths {}; /*<< partial lengths along polyline*/

        /**
         * @brief Determines 1st order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 1st order gradient at given path lengths.
         */
        Points<Eigen::Dynamic> gradient1 (const Eigen::ArrayXd& lengths) const override
        {
            auto indicesGrad { indices(lengths) };
            for(int& idx : indicesGrad)
                idx = idx + 1 > m_lengths.rows() - 1 ? m_lengths.rows() - 2 : idx;
            return { m_x[1](indicesGrad + 1), m_y[1](indicesGrad + 1) };
        }

        /**
         * @brief Determines 2nd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 2nd order gradient at given path lengths.
         */
        Points<Eigen::Dynamic> gradient2 (const Eigen::ArrayXd& lengths) const override
        {
            auto indicesGrad { indices(lengths) };
            for(int& idx : indicesGrad)
                idx = idx + 2 > m_lengths.rows() - 1 ? m_lengths.rows() - 3 : idx;
            return { m_x[2](indicesGrad + 2), m_y[2](indicesGrad + 2) };
        }

        /**
         * @brief Determines 3rd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 3rd order gradient at given path lengths.
         */
        Points<Eigen::Dynamic> gradient3 (const Eigen::ArrayXd& lengths) const override
        {
            auto indicesGrad { indices(lengths) };
            for(int& idx : indicesGrad)
                idx = idx + 3 > m_lengths.rows() - 1 ? m_lengths.rows() - 4 : idx;
            return { m_x[3](indicesGrad + 3), m_y[3](indicesGrad + 3) };
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
                result(row) = FrenetTransform::first(m_lengths(Eigen::seqN(0, m_lengths.rows() - 2)) - lengths(row));

            return result;
        }
    };
};

#endif