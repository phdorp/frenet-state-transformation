#ifndef CIRCLE_H
#define CIRCLE_H

#include <eigen3/Eigen/Core>
#include <math.h>

#include "path.h"
#include "point.h"

using FrenetTransform::Path;
using FrenetTransform::Point;

namespace Testing
{
    class Circle : public Path
    {
    public:
        Circle() = delete;

        Circle(double radius, Point center, double angle0 = 0)
            : m_radius { radius }
            , m_center { center }
            , m_angle0 { angle0 }
        {
        }

        Points<Eigen::Dynamic> operator()(const Eigen::ArrayXd& lengths) const override { return { m_center.x() + m_radius * angle(lengths).cos(), m_center.y() + m_radius * angle(lengths).sin() }; }

        Eigen::ArrayXd lengths(const Points<Eigen::Dynamic>& points) const override { return m_radius * ((points.y() - m_center.y()) / (points.y() - m_center.y())).atan(); }

        Eigen::ArrayXd angle(const Eigen::ArrayXd& lengths) const { return lengths / m_radius + m_angle0; }

    private:
        const double m_radius {};
        const Point m_center;
        const double m_angle0 {};

        Points<Eigen::Dynamic> gradient1(const Eigen::ArrayXd& lengths) const override { return { -angle(lengths).sin(), angle(lengths).cos() }; }

        Points<Eigen::Dynamic> gradient2(const Eigen::ArrayXd& lengths) const override { return { -angle(lengths).cos() / m_radius, -angle(lengths).sin() / m_radius }; }

        Points<Eigen::Dynamic> gradient3(const Eigen::ArrayXd& lengths) const override { return { angle(lengths).sin() / std::pow(m_radius, 2), -angle(lengths).cos() / std::pow(m_radius, 2) }; }
    };
};

#endif