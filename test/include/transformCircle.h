#ifndef TRANSFORM_CIRCLE_H
#define TRANSFORM_CRICLE_H

#include "circle.h"
#include "transform.h"
#include "point.h"

namespace FrenetTransform
{
    namespace Testing
    {
        class TransformCircle : public FrenetTransform::Transform
        {
        public:
            TransformCircle(std::shared_ptr<Circle> circle)
                : Transform(circle)
                , m_path { circle }
            {
            }

        private:
            const std::shared_ptr<Circle> m_path;
        };

        class PointCircle : public Point<PointCircle>
        {
        public:
            PointCircle(double radius, double angle)
                : Point(radius, angle)
            {
            }
        };
    };
};

#endif