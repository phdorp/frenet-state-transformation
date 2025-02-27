#ifndef POINT_H
#define POINT_H

#include <math.h>

namespace FrenetTransform{
    class Point
    {
    public:
        Point() = delete;

        Point(double x, double y)
            : m_x { x }, m_y { y }
        {
        }

        double distance(const Point& point) const { return std::sqrt(std::pow(point.x() - m_x, 2) + std::pow(point.y() - m_y, 2)); }

        double x() const { return m_x; }

        double y() const { return m_y; }

    private:
        double m_x {};
        double m_y {};
    };
};

#endif