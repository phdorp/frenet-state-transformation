#ifndef CIRCLE_H
#define CIRCLE_H

#include <Eigen/Core>
#include <math.h>

#include "frenetTransform/path.h"
#include "frenetTransform/point.h"
#include "frenetTransform/math.h"

namespace FrenetTransform
{
    namespace Testing
    {
        class Circle : public Path
        {
        public:
            Circle() = delete;

            Circle(double radius, PointCartes center, double angle0 = 0)
                : m_radius { radius }
                , m_center { center }
                , m_angle0 { angle0 }
            {
            }

            Points<Eigen::Dynamic, PointCartes> operator()(const Eigen::ArrayXd& lengths) const override { return { m_center.x() + m_radius * angle(lengths).cos(), m_center.y() + m_radius * angle(lengths).sin() }; }

            Eigen::ArrayXd lengths(const Points<Eigen::Dynamic, PointCartes>& points) const override
            {
                const Eigen::ArrayXd dirsx { points.x() - m_center.x() };
                const Eigen::ArrayXd dirsy { points.y() - m_center.y() };
                return m_radius * (angleDir(dirsx, dirsy) - m_angle0);
            }

            Eigen::ArrayXd lengths(const Eigen::ArrayXd& angles) const { return angles * m_radius; }

            Eigen::ArrayXd angle(const Eigen::ArrayXd& lengths) const { return lengths / m_radius + m_angle0; }

            double radius() const { return m_radius; }

            const PointCartes& center() const { return m_center; }

            double angleOffset() const { return m_angle0; }

        private:
            const double m_radius {};
            const PointCartes m_center;
            const double m_angle0 {};

            Points<Eigen::Dynamic, PointCartes> gradient1(const Eigen::ArrayXd& lengths) const override { return { -angle(lengths).sin(), angle(lengths).cos() }; }

            Points<Eigen::Dynamic, PointCartes> gradient2(const Eigen::ArrayXd& lengths) const override { return { -angle(lengths).cos() / m_radius, -angle(lengths).sin() / m_radius }; }

            Points<Eigen::Dynamic, PointCartes> gradient3(const Eigen::ArrayXd& lengths) const override { return { angle(lengths).sin() / std::pow(m_radius, 2), -angle(lengths).cos() / std::pow(m_radius, 2) }; }
        };
    };
};

#endif