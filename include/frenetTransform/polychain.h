#ifndef polychain_H
#define polychain_H

#include <math.h>

#include "frenetTransform/path.h"
#include "frenetTransform/math.h"
#include "frenetTransform/points.h"
#include "frenetTransform/point.h"

namespace FrenetTransform
{
    /**
     * @brief Path representation as polychain.
     * Represents a 2-dimensional path as a polychain.
     * Provide path properties based on finite differences at query points.
     *
     *
     * @tparam NumPoints number of points along the path with -1 for dynamic point number.
     * @tparam NumQueries number of query points with -1 for dynamic point number.
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
         * @param x coordinates in x-direction along the path.
         * @param y coordinates in y-direction along the path.
         */
        Polychain(const ArrayPoints& x, const ArrayPoints& y) { setPoints(x, y); }

        /**
         * @brief Construct a new Polychain object from points in Cartesian coordinates.
         *
         * @param points points along the path.
         */
        Polychain(const Points<NumPoints>& points) { setPoints(points.x(), points.y()); }

        /**
         * @brief Gets points along the polychain at the query lengths.
         * Lengths exceeding the polychain's domain either resolve to the first or last point.
         *
         * @param lengths query lengths along the polychain.
         * @return Points<NumQueries> at the query lengths.
         */
        Points<NumQueries> operator()(const ArrayQueries& lengths) const override
        {
            // indices of corresponding polychain segments
            Eigen::Array<int, NumQueries, 1> indicesLengths { indices(lengths) };
            for(int& idx : indicesLengths)
                idx = idx > 0 ? idx - 1 : idx;

            // relative position along the linear segment
            const ArrayPoints pathLengthsdiff { diffBackward(m_lengths) };
            ArrayQueries segmentPart { (lengths - m_lengths(indicesLengths)) / pathLengthsdiff(indicesLengths + 1) };
            for(double& part : segmentPart)
                part = std::clamp(part, 0.0, 1.0);

            // absolute position along path
            ArrayQueries x { m_points[0].x()(indicesLengths + 1) * segmentPart + m_points[0].x()(indicesLengths) * (1 - segmentPart) };
            ArrayQueries y { m_points[0].y()(indicesLengths + 1) * segmentPart + m_points[0].y()(indicesLengths) * (1 - segmentPart) };

            return { x, y };
        }

        /**
         * @brief Determines next points to the query points.
         * Performs a linear search over all polychain segments to identify the closest one.
         *
         * @param points query points.
         * @return Points<NumQueries> next points to query points.
         */
        ArrayQueries lengths(const Points<NumQueries>& points) const override
        {
            ArrayQueries lengthsPoints (points.numPoints());

            // determine lengths for all query points
            for(int cQuery {}; cQuery < points.numPoints(); ++cQuery)
            {
                // current minimum squared distance between query point and polychain
                double distanceSquare { -1.0 };

                // determine distance between query point and path segment
                for(int cPoints {1}; cPoints < m_numPoints; ++cPoints)
                {
                    // point at end of current segment
                    const Point nextPoint { m_points[0](cPoints) };
                    // difference between "nextPoint" and query point
                    const Point diffPoint { nextPoint - points(cQuery) };
                    // parameter determining next point along current linear segment
                    const double segmentPart { (diffPoint.x() * m_xDiff(cPoints) + diffPoint.y() * m_yDiff(cPoints)) / m_diffSquare(cPoints) };

                    // squared distance between current segment and query point
                    double distanceSquareCand {};
                    // length along polychain of shortest distance point on segment to query point
                    double lengthCand {};

                    // determine squared distance and length if shortest distance point is beyond segment end
                    if(segmentPart >= 1.0)
                    {
                        distanceSquareCand = m_points[0](cPoints - 1).distanceSquare(points(cQuery));
                        lengthCand = m_lengths.data()[cPoints - 1];
                    }
                    // determine squared distance and length if shortest distance point is beyond segment start
                    else if(segmentPart <= 0.0)
                    {
                        distanceSquareCand = m_points[0](cPoints).distanceSquare(points(cQuery));
                        lengthCand = m_lengths.data()[cPoints];
                    }
                    // determine squared distance and length if shortest distance point is on segment
                    else
                    {
                        // point at start of current segment
                        const Point prevPoint { m_points[0](cPoints - 1) };
                        // shortest distance point on current segment
                        const Point pointCand {
                            prevPoint.x() * segmentPart + (1 - segmentPart) * nextPoint.x(),
                            prevPoint.y() * segmentPart + (1 - segmentPart) * nextPoint.y()
                        };
                        distanceSquareCand = points(cQuery).distanceSquare(pointCand);
                        lengthCand = m_lengths(cPoints - 1) +  pointCand.distance(prevPoint);
                    }

                    // update length and squared distance if squared distance is smaller than current minimum
                    if(distanceSquareCand < distanceSquare || distanceSquare < 0)
                    {
                        distanceSquare = distanceSquareCand;
                        lengthsPoints(cQuery) = lengthCand;
                    }
                }
            }

            return lengthsPoints;
        }

        /**
         * @brief Provide new points for the polychain.
         * Update lengths and gradient information.
         *
         * @param x coordinates in x-direction of new points.
         * @param y coordinates in y-direction of new points.
         */
        void setPoints(const ArrayPoints& x, const ArrayPoints& y)
        {
            // number of points along the polychain
            m_numPoints = x.rows();

            // update points along the polychain
            m_points[0] = Points<NumPoints> { x, y };

            // point differences
            m_xDiff = FrenetTransform::diffBackward(x);
            m_yDiff = FrenetTransform::diffBackward(y);
            // point squared differences
            m_diffSquare = m_xDiff.pow(2) + m_yDiff.pow(2);

            // accumulated lengths along the new points
            m_lengths = FrenetTransform::partialLength(x, y);

            // update gradients
            for(unsigned int orderGrad { 1 }; orderGrad < s_numGrad; ++orderGrad)
                m_points[orderGrad] = Points<NumPoints> {
                    FrenetTransform::gradient(m_points[orderGrad - 1].x(), m_lengths),
                    FrenetTransform::gradient(m_points[orderGrad - 1].y(), m_lengths)
                };
        }

    private:
        int m_numPoints {}; /*<< number of points along the polychain*/

        ArrayPoints m_lengths {}; /*<< partial lengths along polychain*/

        ArrayPoints m_xDiff {}; /*<< differences between polychain points in x-direction*/
        ArrayPoints m_yDiff {}; /*<< differences between polychain points in y-direction*/
        ArrayPoints m_diffSquare {}; /*<< differences between polychain points squared*/

        static constexpr int s_numGrad { 4 }; /*<< number of time derivatives provided*/
        std::array<Points<NumPoints>, s_numGrad> m_points {}; /*<< points and gradients at polychain points*/

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
                idx = std::clamp(idx, 1, m_numPoints);
            return { m_points[1].x()(indicesGrad),m_points[1].y()(indicesGrad) };
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
                idx = std::clamp(idx, 2, m_numPoints);
            return { m_points[2].x()(indicesGrad),m_points[2].y()(indicesGrad) };
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
                idx = std::clamp(idx, 3, m_numPoints);
            return { m_points[3].x()(indicesGrad),m_points[3].y()(indicesGrad) };
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
                result(row) = FrenetTransform::first(m_lengths - lengths(row));

            return result;
        }
    };
};

#endif