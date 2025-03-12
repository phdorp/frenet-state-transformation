#ifndef POLYLINE_H
#define POLYLINE_H

#include <math.h>

#include "path.h"
#include "math.h"
#include "points.h"
#include "point.h"

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
        using ArrayT1 = Eigen::Array<double, T, 1>;
        using ArrayT2 = Eigen::Array<double, T, 2>;

        Polyline() = default;

        /**
         * @brief Construct a new Path object from Cartesian x- and y-positions.
         *
         * @param x coordinates in x-direction.
         * @param y coordinates in y-direction.
         */
        Polyline(const ArrayT1& x, const ArrayT1& y) { setPoints(x, y); }

        Polyline(const Points<T, PointCartes>& points) { setPoints(points.x(), points.y()); }

        Points<Eigen::Dynamic, PointCartes> operator()(const Eigen::ArrayXd& lengths) const override
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
        Eigen::ArrayXd lengths(const Points<Eigen::Dynamic, PointCartes>& points) const override
        {
            Eigen::ArrayXd lLenghts (points.numPoints());

            for(int cPoints {}; cPoints < points.numPoints(); ++cPoints)
            {
                const auto distSqCorners { (points(cPoints).x() - m_x[0]).pow(2) + (points(cPoints).y() - m_y[0]).pow(2) };

                double distSq { -1.0 };

                for(int cCorder {}; cCorder < distSqCorners.rows(); ++cCorder)
                {
                    if(distSq < 0.0 ||  distSqCorners(cCorder) < distSq)
                    {
                        distSq = distSqCorners(cCorder);
                        lLenghts(cPoints) = m_lengths(cCorder);
                    }
                }

                for(int cThis {}; cThis < T; ++cThis)
                {
                    const PointCartes pathNext { m_x[0].data()[cThis], m_y[0].data()[cThis] };
                    const double xDiffPnt { pathNext.x() - points.x(cPoints) };
                    const double yDiffPnt { pathNext.y() - points.y(cPoints) };
                    const double param { (xDiffPnt * m_xDiff(cThis) + yDiffPnt * m_yDiff(cThis))
                        / (std::pow(m_xDiff(cThis), 2) + std::pow(m_yDiff(cThis), 2)) };

                    if(param >= 0.0 && param <= 1.0)
                    {
                        const PointCartes pathPrev { m_x[0].data()[cThis - 1], m_y[0].data()[cThis - 1] };
                        const double xCand { pathPrev.x() * param + (1 - param) * pathNext.x() };
                        const double yCand { pathPrev.y() * param + (1 - param) * pathNext.y() };
                        const PointCartes cand { xCand, yCand };
                        const double distSqCand { std::pow(xCand - points.x(cPoints), 2) + std::pow(yCand - points.y(cPoints), 2) };

                        if(distSq < 0.0 || distSq > distSqCand)
                        {
                            distSq = distSqCand;
                            lLenghts(cPoints) = m_lengths(cThis - 1) + cand.distance(pathPrev);
                        }
                    }
                }
            }

            return lLenghts;
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
        Points<Eigen::Dynamic, PointCartes> gradient1 (const Eigen::ArrayXd& lengths) const override
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
        Points<Eigen::Dynamic, PointCartes> gradient2 (const Eigen::ArrayXd& lengths) const override
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
        Points<Eigen::Dynamic, PointCartes> gradient3 (const Eigen::ArrayXd& lengths) const override
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