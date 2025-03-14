#ifndef TRANSFORM_CIRCLE_H
#define TRANSFORM_CRICLE_H

#include "circle.h"
#include "frenetTransform/transform.h"
#include "frenetTransform/point.h"

namespace FrenetTransform
{
    namespace Testing
    {
        class PointCircle : public Point<PointCircle>
        {
        public:
            PointCircle(double radius, double angle)
                : Point(radius, angle)
            {
            }
        };

        class TransformCircle : public FrenetTransform::Transform
        {
        public:
            TransformCircle(std::shared_ptr<Circle> circle)
                : Transform(circle)
                , m_path { circle }
            {
            }

            Points<Eigen::Dynamic, PointFrenet> posFrenet(const Points<Eigen::Dynamic, PointCircle>& posCircle) const { return { m_path->radius() * (posCircle.y() - m_path->angleOffset()), m_path->radius() - posCircle.x() }; }

            Points<Eigen::Dynamic, PointFrenet> velFrenet(const Points<Eigen::Dynamic, PointCircle>& velCircle) const { return { m_path->radius() * velCircle.y(), -velCircle.x() }; }

            Points<Eigen::Dynamic, PointFrenet> accFrenet(const Points<Eigen::Dynamic, PointCircle>& accCircle) const { return { m_path->radius() * accCircle.y(), -accCircle.x() }; }

            Points<Eigen::Dynamic, PointCartes> posCartes(const Points<Eigen::Dynamic, PointCircle>& posCircle) const { return { posCircle.x() * posCircle.y().cos() + m_path->center().x(), posCircle.x() * posCircle.y().sin() + m_path->center().y() }; }

            Points<Eigen::Dynamic, PointCartes> velCartes(const Points<Eigen::Dynamic, PointCircle>& velCircle, const Points<Eigen::Dynamic, PointCircle>& posCircle) const
            {
                return { velCircle.x() * posCircle.y().cos() - posCircle.x() * velCircle.y() * posCircle.y().sin(),  velCircle.x() * posCircle.y().sin() + posCircle.x() * velCircle.y() * posCircle.y().cos() };
            }

            Points<Eigen::Dynamic, PointCartes> accCartes(const Points<Eigen::Dynamic, PointCircle>& accCircle, const Points<Eigen::Dynamic, PointCircle>& velCircle, const Points<Eigen::Dynamic, PointCircle>& posCircle) const
            {
                return { (accCircle.x() - posCircle.x() * velCircle.y().pow(2)) * posCircle.y().cos() - (2 * velCircle.x() * velCircle.y() + posCircle.x() * accCircle.y()) * posCircle.y().sin(),
                    (accCircle.x() - posCircle.x() * velCircle.y().pow(2)) * posCircle.y().sin() + (2 * velCircle.x() * velCircle.y() + posCircle.x() * accCircle.y()) * posCircle.y().cos() };
            }

        private:
            const std::shared_ptr<Circle> m_path;
        };
    };
};

#endif