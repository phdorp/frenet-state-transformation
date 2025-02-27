#include "transform.h"

using FrenetTransform::Points;
using FrenetTransform::Transform;

Points<Eigen::Dynamic> Transform::posFrenet(const Points<Eigen::Dynamic>& posCartes) const
{
    const auto lengths { m_path->lengths(posCartes) };
    const auto posPath { m_path->operator()(lengths) };
    const auto posDiff { posCartes - posPath };
    const auto normals { m_path->normal(lengths) };
    return { lengths, normals * posDiff };
}