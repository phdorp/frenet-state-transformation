#ifndef POINTS_H
#define POINTS_H

#include <eigen3/Eigen/Core>
#include <cassert>

#include "point.h"

namespace FrenetTransform{
    template <int NumPoints>
    class Points
    {
    public:
        Points() = delete;

        Points(const Eigen::Array<double, NumPoints, 1>& x, const Eigen::Array<double, NumPoints, 1>& y)
            : m_x { x }, m_y { y }
        {
            assert(x.cols() == y.cols());
        }

        Point operator()(const int index) const { return { m_x(index), m_y(index) }; }

        int numPoints() const { return NumPoints > -1 ? NumPoints : m_x.rows(); }

        const Eigen::Array<double, NumPoints, 1>& x() const { return m_x; }

        double x(const int index) const { return m_x(index); }

        const Eigen::Array<double, NumPoints, 1>& y() const { return m_y; }

        double y(const int index) const { return m_y(index); }

    private:
        const Eigen::Array<double, NumPoints, 1> m_x {};
        const Eigen::Array<double, NumPoints, 1> m_y {};
    };
};

#endif