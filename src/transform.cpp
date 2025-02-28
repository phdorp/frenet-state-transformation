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

Points<Eigen::Dynamic> Transform::velFrenet(const Points<Eigen::Dynamic>& velCartes, const Points<Eigen::Dynamic>& posFrenet) const
{
    const auto velTransformsInv { Transform::velTransformInv(posFrenet) };
    return {
        velTransformsInv(0, 0) * velCartes.x() + velTransformsInv(0, 1) * velCartes.y(),
        velTransformsInv(1, 0) * velCartes.x() + velTransformsInv(1, 1) * velCartes.y()
    };
}

Points<Eigen::Dynamic> Transform::velCartes(const Points<Eigen::Dynamic>& velFrenet, const Points<Eigen::Dynamic>& posFrenet) const
{
    const auto velTransforms { Transform::velTransform(posFrenet) };
    return {
        velTransforms(0, 0) * velFrenet.x() + velTransforms(0, 1) * velFrenet.y(),
        velTransforms(1, 0) * velFrenet.x() + velTransforms(1, 1) * velFrenet.y()
    };
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

Eigen::Array<Eigen::ArrayXd, 2, 2> Transform::velTransformInv(const Points<Eigen::Dynamic>& posFrenet) const
{
    const auto velTransforms { Transform::velTransform(posFrenet) };
    const auto normalization { 1 / (velTransforms(0, 0) * velTransforms(1, 1) - velTransforms(1, 0) * velTransforms(0, 1)) };
    return {
        { velTransforms(1, 1) * normalization, -velTransforms(0, 1) * normalization},
        {-velTransforms(1, 0) * normalization,  velTransforms(0, 0) * normalization}
    };
}