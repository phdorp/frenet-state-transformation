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
        template<int NumQueries=Eigen::Dynamic>
        class Circle : public Path<NumQueries>
        {
        public:
            using ArrayQueries = Eigen::Array<double, NumQueries, 1>;

            Circle() = delete;

            Circle(double radius, Point center, double angle0 = 0)
                : m_radius { radius }
                , m_center { center }
                , m_angle0 { angle0 }
            {
            }

            Points<NumQueries> operator()(const ArrayQueries& lengths) const override { return { m_center.x() + m_radius * angle(lengths).cos(), m_center.y() + m_radius * angle(lengths).sin() }; }

            ArrayQueries lengths(const Points<NumQueries>& points) const override
            {
                const ArrayQueries dirsx { points.x() - m_center.x() };
                const ArrayQueries dirsy { points.y() - m_center.y() };
                return m_radius * (angleDir(dirsx, dirsy) - m_angle0);
            }

            ArrayQueries lengths(const ArrayQueries& angles) const { return angles * m_radius; }

            ArrayQueries angle(const ArrayQueries& lengths) const { return lengths / m_radius + m_angle0; }

            double radius() const { return m_radius; }

            const Point& center() const { return m_center; }

            double angleOffset() const { return m_angle0; }

        private:
            const double m_radius {};
            const Point m_center;
            const double m_angle0 {};

            Points<NumQueries> gradient1(const ArrayQueries& lengths) const override { return { -angle(lengths).sin(), angle(lengths).cos() }; }

            Points<NumQueries> gradient2(const ArrayQueries& lengths) const override { return { -angle(lengths).cos() / m_radius, -angle(lengths).sin() / m_radius }; }

            Points<NumQueries> gradient3(const ArrayQueries& lengths) const override { return { angle(lengths).sin() / std::pow(m_radius, 2), -angle(lengths).cos() / std::pow(m_radius, 2) }; }
        };
    };
};

#endif