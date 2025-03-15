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
     * @brief Path representation as polychain.
     *
     * Represents a 2-dimensional path as a polychain.
     * Provide path properties based on finite differences at query points.
     *
     * @tparam NumPoints determines the number points along the path.
     */
    template <int NumPoints, int NumQueries=Eigen::Dynamic>
    class Polychain : public Path<NumQueries>
    {
    public:
        using ArrayPoints = Eigen::Array<double, NumPoints, 1>;
        using ArrayQueries = Eigen::Array<double, NumQueries, 1>;

        Polychain() = default;

        /**
         * @brief Construct a new Path object from Cartesian x- and y-positions.
         *
         * @param x coordinates in x-direction.
         * @param y coordinates in y-direction.
         */
        Polychain(const ArrayPoints& x, const ArrayPoints& y) { setPoints(x, y); }

        Polychain(const Points<NumPoints>& points) { setPoints(points.x(), points.y()); }

        Points<NumQueries> operator()(const ArrayQueries& lengths) const override
        {
            // indices of corresponding polychain segments
            const Eigen::Array<int, NumQueries, 1> indicesLengths { indices(lengths) };

            // relative position along the linear segment
            const ArrayPoints pathLengthsdiff { diffBackward(m_lengths) };
            ArrayQueries relativePos { (lengths - m_lengths(indicesLengths)) / pathLengthsdiff(indicesLengths + 1) };

            for(double& pos : relativePos)
                pos = std::clamp(pos, 0.0, 1.0);

            // absolute position along path
            ArrayQueries x { m_points[0].x()(indicesLengths + 1) * relativePos + m_points[0].x()(indicesLengths) * (1 - relativePos) };
            ArrayQueries y { m_points[0].y()(indicesLengths + 1) * relativePos + m_points[0].y()(indicesLengths) * (1 - relativePos) };

            for(int row {}; row < indicesLengths.rows(); ++row)
            {
                if(indicesLengths(row) >= m_lengths.rows() - 1)
                {
                    x(row) = m_points[0].x(m_lengths.rows() - 1);
                    y(row) = m_points[0].y(m_lengths.rows() - 1);
                }
            }

            return { x, y };
        }

        /**
         * @brief Determines next points to the query points.
         *
         * @param points query points.
         * @return Points<NumQueries> next to query points.
         */
        ArrayQueries lengths(const Points<NumQueries>& points) const override
        {
            ArrayQueries lLenghts (points.numPoints());

            for(int cQuery {}; cQuery < points.numPoints(); ++cQuery)
            {
                double distSq { -1.0 };

                for(int cPoints {1}; cPoints < m_numPoints; ++cPoints)
                {
                    const Point pathNext { m_points[0](cPoints) };
                    const Point diffPoint { pathNext - points(cQuery) };
                    const double segmentPart { (diffPoint.x() * m_xDiff(cPoints) + diffPoint.y() * m_yDiff(cPoints)) / m_diffSq(cPoints) };

                    double distSqCand {};
                    double candLength {};

                    if(segmentPart >= 1.0)
                    {
                        distSqCand = m_points[0](cPoints - 1).distanceSquare(points(cQuery));
                        candLength = m_lengths.data()[cPoints - 1];
                    }
                    else if(segmentPart <= 0.0)
                    {
                        distSqCand = m_points[0](cPoints).distanceSquare(points(cQuery));
                        candLength = m_lengths.data()[cPoints];
                    }
                    else
                    {
                        const Point pathPrev { m_points[0](cPoints - 1) };
                        const Point pointCand {
                            pathPrev.x() * segmentPart + (1 - segmentPart) * pathNext.x(),
                            pathPrev.y() * segmentPart + (1 - segmentPart) * pathNext.y()
                        };
                        distSqCand = points(cQuery).distanceSquare(pointCand);
                        candLength = m_lengths(cPoints - 1) +  pointCand.distance(pathPrev);
                    }

                    if(distSqCand < distSq || distSq < 0)
                    {
                        distSq = distSqCand;
                        lLenghts(cQuery) = candLength;
                    }
                }
            }

            return lLenghts;
        }

        void setPoints(const ArrayPoints& x, const ArrayPoints& y)
        {
            m_numPoints = x.rows();

            m_points[0] = Points<NumPoints> { x, y };

            m_xDiff = FrenetTransform::diffBackward(x);
            m_yDiff = FrenetTransform::diffBackward(y);

            m_diffSq = m_xDiff.pow(2) + m_yDiff.pow(2);
            m_lengths = FrenetTransform::partialLength(x, y);

            for(unsigned int orderGrad { 1 }; orderGrad < s_numGrad; ++orderGrad)
                m_points[orderGrad] = Points<NumPoints> {
                    FrenetTransform::gradient(m_points[orderGrad - 1].x(), m_lengths),
                    FrenetTransform::gradient(m_points[orderGrad - 1].y(), m_lengths)
                };
        }

        int numPoints() { return m_numPoints; }

    private:
        static constexpr int s_numGrad { 4 };
        ArrayPoints m_xDiff {};
        ArrayPoints m_yDiff {};
        ArrayPoints m_diffSq {};
        std::array<Points<NumPoints>, s_numGrad> m_points {};
        ArrayPoints m_lengths {}; /*<< partial lengths along polychain*/
        int m_numPoints {};

        /**
         * @brief Determines 1st order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 1st order gradient at given path lengths.
         */
        Points<NumQueries> gradient1 (const ArrayQueries& lengths) const override
        {
            auto indicesGrad { indices(lengths) };
            for(int& idx : indicesGrad)
                idx = idx + 1 > m_lengths.rows() - 1 ? m_lengths.rows() - 2 : idx;
            return { m_points[1].x()(indicesGrad + 1),m_points[1].y()(indicesGrad + 1) };
        }

        /**
         * @brief Determines 2nd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 2nd order gradient at given path lengths.
         */
        Points<NumQueries> gradient2 (const ArrayQueries& lengths) const override
        {
            auto indicesGrad { indices(lengths) };
            for(int& idx : indicesGrad)
                idx = idx + 2 > m_lengths.rows() - 1 ? m_lengths.rows() - 3 : idx;
            return { m_points[2].x()(indicesGrad + 1),m_points[2].y()(indicesGrad + 1) };
        }

        /**
         * @brief Determines 3rd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 3rd order gradient at given path lengths.
         */
        Points<NumQueries> gradient3 (const ArrayQueries& lengths) const override
        {
            auto indicesGrad { indices(lengths) };
            for(int& idx : indicesGrad)
                idx = idx + 3 > m_lengths.rows() - 1 ? m_lengths.rows() - 4 : idx;
            return { m_points[3].x()(indicesGrad + 1),m_points[3].y()(indicesGrad + 1) };
        }

        /**
         * @brief Determines indices of polychain segment corresponding to the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return indices corresponding to given path lengths.
         */
        Eigen::Array<int, NumQueries, 1> indices(const ArrayQueries& lengths) const
        {
            Eigen::Array<int, NumQueries, 1> result(lengths.rows()); // vector of segment indices

            // get indices of next segments
            for(int row {}; row < lengths.rows(); ++row)
                result(row) = FrenetTransform::first(m_lengths(Eigen::seqN(0, m_lengths.rows() - 2)) - lengths(row));

            return result;
        }
    };
};

#endif