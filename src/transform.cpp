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

Points<Eigen::Dynamic> Transform::posCartes(const Points<Eigen::Dynamic>& posFrenet) const
{
    const auto posPath { m_path->operator()(posFrenet.x()) };
    const auto normals { m_path->normal(posFrenet.x()) };
    return posPath + normals * posFrenet.y();
}

Eigen::Array<Eigen::ArrayXd, 2, 2> Transform::velTransform(const Points<Eigen::Dynamic>& posFrenet) const
{
    const auto tangents { m_path->tangent(posFrenet.x()) };
    const auto normals { m_path->normal(posFrenet.x()) };
    const auto curvs { m_path->angle1(posFrenet.x()) };
    return {
        {tangents.x() * (1 - curvs * posFrenet.y()), normals.x()},
        {                              tangents.y(), normals.y()}
    };
}