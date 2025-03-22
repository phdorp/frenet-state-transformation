#ifndef TRANSFORM_CIRCLE_H
#define TRANSFORM_CRICLE_H

#include "frenetTransform/internal/circle.h"
#include "frenetTransform/transform.h"
#include "frenetTransform/point.h"

namespace FrenetTransform
{
    namespace Testing
    {
        class PointCircle : public Point
        {
        public:
            PointCircle(double radius, double angle)
                : Point(radius, angle)
            {
            }
        };

        template <int NumQueries=Eigen::Dynamic>
        class TransformCircle : public Transform<NumQueries>
        {
        public:
            TransformCircle(std::shared_ptr<Circle<NumQueries>> circle)
                : Transform<NumQueries>(circle)
                , m_path { circle }
            {
            }

            Points<NumQueries> posFrenet(const Points<NumQueries, PointCircle>& posCircle) const { return { m_path->radius() * (posCircle.y() - m_path->angleOffset()), m_path->radius() - posCircle.x() }; }

            Points<NumQueries> velFrenet(const Points<NumQueries, PointCircle>& velCircle) const { return { m_path->radius() * velCircle.y(), -velCircle.x() }; }

            Points<NumQueries> accFrenet(const Points<NumQueries, PointCircle>& accCircle) const { return { m_path->radius() * accCircle.y(), -accCircle.x() }; }

            Points<NumQueries> posCartes(const Points<NumQueries, PointCircle>& posCircle) const { return { posCircle.x() * posCircle.y().cos() + m_path->center().x(), posCircle.x() * posCircle.y().sin() + m_path->center().y() }; }

            Points<NumQueries> velCartes(const Points<NumQueries, PointCircle>& velCircle, const Points<NumQueries, PointCircle>& posCircle) const
            {
                return { velCircle.x() * posCircle.y().cos() - posCircle.x() * velCircle.y() * posCircle.y().sin(),  velCircle.x() * posCircle.y().sin() + posCircle.x() * velCircle.y() * posCircle.y().cos() };
            }

            Points<NumQueries> accCartes(const Points<NumQueries, PointCircle>& accCircle, const Points<NumQueries, PointCircle>& velCircle, const Points<NumQueries, PointCircle>& posCircle) const
            {
                return { (accCircle.x() - posCircle.x() * velCircle.y().pow(2)) * posCircle.y().cos() - (2 * velCircle.x() * velCircle.y() + posCircle.x() * accCircle.y()) * posCircle.y().sin(),
                    (accCircle.x() - posCircle.x() * velCircle.y().pow(2)) * posCircle.y().sin() + (2 * velCircle.x() * velCircle.y() + posCircle.x() * accCircle.y()) * posCircle.y().cos() };
            }

        private:
            const std::shared_ptr<Circle<NumQueries>> m_path;
        };
    };
};

#endif