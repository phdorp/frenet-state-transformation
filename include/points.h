#ifndef POINTS_H
#define POINTS_H

#include <Eigen/Core>
#include <cassert>

#include "point.h"

namespace FrenetTransform{
    template <int NumPoints, typename PointType>
    class Points
    {
    public:
        Points () {}

        Points(const Eigen::Array<double, NumPoints, 1>& x, const Eigen::Array<double, NumPoints, 1>& y)
            : m_x { x }, m_y { y }
        {
            assert(x.cols() == y.cols());
        }

        Point<PointType> operator()(const int index) const { return { m_x(index), m_y(index) }; }

        int numPoints() const { return NumPoints > -1 ? NumPoints : m_x.rows(); }

        template<typename P>
        Eigen::Array<double, NumPoints, 1> distance(const Point<PointType>& point) { return ((m_x - point.x()).pow(2) + (m_y - point.y()).pow(2)).sqrt(); }

        const Eigen::Array<double, NumPoints, 1>& x() const { return m_x; }

        double x(const int index) const { return m_x(index); }

        const Eigen::Array<double, NumPoints, 1>& y() const { return m_y; }

        double y(const int index) const { return m_y(index); }

        Points operator-() const { return { -m_x, -m_y }; }

        template<int T>
        friend Points<T, PointType> operator+(const Points<T, PointType>& points1, const Points<T, PointType>& points2) { return { points1.m_x + points2.m_x, points1.m_y + points2.m_y }; }

        template<int T>
        friend Points<T, PointType> operator+(const Points<T, PointType>& points, const Point<PointType>& point) { return { point.m_x + point.x(), points.m_y + point.y() }; }

        template<int T>
        friend Points<T, PointType> operator+(const Point<PointType>& point, const Points<T, PointType>& points) { return points + point; }

        template<int T>
        friend Points<T, PointType> operator-(const Points<T, PointType>& points1, const Points<T, PointType>& points2) { return -points2 + points1; }

        template<int T>
        friend Points<T, PointType> operator-(const Points<T, PointType>& points, const Point<PointType>& point) { return -point + points; }

        template<int T>
        friend Points<T, PointType> operator-(const Point<PointType>& point, const Points<T, PointType>& points) { return -points + point; }

        template<int T>
        friend Eigen::Array<double, T, 1> operator*(const Points<T, PointType>& points1, const Points<T, PointType>& points2) { return points1.m_x * points2.m_x + points1.m_y * points2.m_y; }

        template<int T>
        friend Points<T, PointType> operator*(const Points<T, PointType>& points, const Eigen::Array<double, T, 1>& nums) { return { points.m_x * nums, points.m_y * nums }; }

        template<int T>
        friend Points<T, PointType> operator*(const Eigen::Array<double, T, 1>& nums, const Points<T, PointType>& points) { return points * nums; }

    private:
        Eigen::Array<double, NumPoints, 1> m_x {};
        Eigen::Array<double, NumPoints, 1> m_y {};
    };
};

#endif