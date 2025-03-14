#ifndef POINTS_H
#define POINTS_H

#include <Eigen/Core>
#include <cassert>

#include "frenetTransform/point.h"

namespace FrenetTransform{
    template <int NumPoints, typename PointType=Point>
    class Points
    {
    public:
        Points () {}

        Points(const Eigen::Array<double, NumPoints, 1>& x, const Eigen::Array<double, NumPoints, 1>& y)
            : m_x { x }, m_y { y }
        {
            assert(x.cols() == y.cols());
        }

        PointType operator()(const int index) const { return { m_x(index), m_y(index) }; }

        int numPoints() const { return NumPoints > -1 ? NumPoints : m_x.rows(); }

        template<typename P>
        Eigen::Array<double, NumPoints, 1> distance(const PointType& point) { return ((m_x - point.x()).pow(2) + (m_y - point.y()).pow(2)).sqrt(); }

        const Eigen::Array<double, NumPoints, 1>& x() const { return m_x; }

        double x(const int index) const { return m_x(index); }

        const Eigen::Array<double, NumPoints, 1>& y() const { return m_y; }

        double y(const int index) const { return m_y(index); }

        Points operator-() const { return { -m_x, -m_y }; }

        friend Points<NumPoints, PointType> operator+(const Points<NumPoints, PointType>& points1, const Points<NumPoints, PointType>& points2) { return { points1.m_x + points2.m_x, points1.m_y + points2.m_y }; }

        friend Points<NumPoints, PointType> operator+(const Points<NumPoints, PointType>& points, const PointType& point) { return { point.m_x + point.x(), points.m_y + point.y() }; }

        friend Points<NumPoints, PointType> operator+(const PointType& point, const Points<NumPoints, PointType>& points) { return points + point; }

        friend Points<NumPoints, PointType> operator-(const Points<NumPoints, PointType>& points1, const Points<NumPoints, PointType>& points2) { return -points2 + points1; }

        friend Points<NumPoints, PointType> operator-(const Points<NumPoints, PointType>& points, const PointType& point) { return -point + points; }

        friend Points<NumPoints, PointType> operator-(const PointType& point, const Points<NumPoints, PointType>& points) { return -points + point; }

        friend Eigen::Array<double, NumPoints, 1> operator*(const Points<NumPoints, PointType>& points1, const Points<NumPoints, PointType>& points2) { return points1.m_x * points2.m_x + points1.m_y * points2.m_y; }

        friend Points<NumPoints, PointType> operator*(const Points<NumPoints, PointType>& points, const Eigen::Array<double, NumPoints, 1>& nums) { return { points.m_x * nums, points.m_y * nums }; }

        friend Points<NumPoints, PointType> operator*(const Eigen::Array<double, NumPoints, 1>& nums, const Points<NumPoints, PointType>& points) { return points * nums; }

    private:
        Eigen::Array<double, NumPoints, 1> m_x {};
        Eigen::Array<double, NumPoints, 1> m_y {};
    };
};

#endif