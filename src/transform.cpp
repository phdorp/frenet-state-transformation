#include "transform.h"

using FrenetTransform::Points;
using FrenetTransform::Transform;

namespace FrenetTransform
{
Points<Eigen::Dynamic, PointFrenet> Transform::posFrenet(const Points<Eigen::Dynamic, PointCartes>& posCartes) const
{
    const auto lengths { m_path->lengths(posCartes) };
    const auto posPath { m_path->operator()(lengths) };
    const auto posDiff { posCartes - posPath };
    const auto normals { m_path->normal(lengths) };
    return { lengths, normals * posDiff };
}

Points<Eigen::Dynamic, PointCartes> Transform::posCartes(const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const
{
    const auto posPath { m_path->operator()(posFrenet.x()) };
    const auto normals { m_path->normal(posFrenet.x()) };
    return posPath + normals * posFrenet.y();
}

Points<Eigen::Dynamic, PointFrenet> Transform::velFrenet(const Points<Eigen::Dynamic, PointCartes>& velCartes, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const
{
    const auto velTransformsInv { Transform::transformInv(Transform::velTransform(posFrenet)) };
    return {
        velTransformsInv(0, 0) * velCartes.x() + velTransformsInv(0, 1) * velCartes.y(),
        velTransformsInv(1, 0) * velCartes.x() + velTransformsInv(1, 1) * velCartes.y()
    };
}

Points<Eigen::Dynamic, PointCartes> Transform::velCartes(const Points<Eigen::Dynamic, PointFrenet>& velFrenet, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const
{
    const auto velTransforms { Transform::velTransform(posFrenet) };
    return {
        velTransforms(0, 0) * velFrenet.x() + velTransforms(0, 1) * velFrenet.y(),
        velTransforms(1, 0) * velFrenet.x() + velTransforms(1, 1) * velFrenet.y()
    };
}

Points<Eigen::Dynamic, PointFrenet> Transform::accFrenet(const Points<Eigen::Dynamic, PointCartes>& accCartes, const Points<Eigen::Dynamic, PointFrenet>& velFrenet, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const
{
    const auto accTransform { Transform::accTransform(velFrenet, posFrenet) };
    const Points<Eigen::Dynamic, PointFrenet> accRel {
        accCartes.x() - accTransform(0, 0) * velFrenet.x() - accTransform(0, 1) * velFrenet.y(),
        accCartes.y() - accTransform(1, 0) * velFrenet.x() - accTransform(1, 1) * velFrenet.y()
    };

    const auto velTransformsInv { Transform::transformInv(Transform::velTransform(posFrenet)) };

    return {
        velTransformsInv(0, 0) * accRel.x() * velTransformsInv(0, 1) * accRel.y(),
        velTransformsInv(1, 0) * accRel.x() * velTransformsInv(1, 1) * accRel.y()
    };
}


Points<Eigen::Dynamic, PointCartes> Transform::accCartes(const Points<Eigen::Dynamic, PointFrenet>& accFrenet, const Points<Eigen::Dynamic, PointFrenet>& velFrenet, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const
{
    const auto accTransform { Transform::accTransform(velFrenet, posFrenet) };
    const auto velTransform { Transform::velTransform(posFrenet) };
    return {
        velTransform(0, 0) * accFrenet.x() + velTransform(0, 1) * accFrenet.y() + accTransform(0, 0) * velFrenet.x() + accTransform(0, 1) * velFrenet.y(),
        velTransform(1, 0) * accFrenet.x() + velTransform(1, 1) * accFrenet.y() + accTransform(1, 0) * velFrenet.x() + accTransform(1, 1) * velFrenet.y()
    };
}

Eigen::Array<Eigen::ArrayXd, 2, 2> Transform::velTransform(const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const
{
    const auto tangents { m_path->tangent(posFrenet.x()) };
    const auto normals { m_path->normal(posFrenet.x()) };
    const auto curvs { m_path->angle1(posFrenet.x()) };
    return {
        {tangents.x() * (1 - curvs * posFrenet.y()), normals.x()},
        {tangents.y() * (1 - curvs * posFrenet.y()), normals.y()}
    };
}

Eigen::Array<Eigen::ArrayXd, 2, 2> Transform::accTransform(const Points<Eigen::Dynamic, PointFrenet>& velFrenet, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const
{
    const auto tangents { m_path->tangent(posFrenet.x()) };
    const auto normals { m_path->normal(posFrenet.x()) };

    const auto curvs { m_path->angle1(posFrenet.x()) };
    const auto latScale { 1 - curvs * posFrenet.y() };

    const auto curv1s { m_path->angle2(posFrenet.x()) };
    const auto latScaleDer { curv1s * velFrenet.x() * posFrenet.y() + curvs * velFrenet.y() };

    return {
        {normals.x() * curvs * latScale * velFrenet.x() - tangents.x() * latScaleDer, -curvs * tangents.x() * velFrenet.x()},
        {normals.y() * curvs * latScale * velFrenet.x() - tangents.y() * latScaleDer, -curvs * tangents.y() * velFrenet.x()}
    };
}

Eigen::Array<Eigen::ArrayXd, 2, 2> Transform::transformInv(const Eigen::Array<Eigen::ArrayXd, 2, 2>& transform)
{
    const auto normalization { 1 / (transform(0, 0) * transform(1, 1) - transform(1, 0) * transform(0, 1)) };
    return {
        { transform(1, 1) * normalization, -transform(0, 1) * normalization},
        {-transform(1, 0) * normalization,  transform(0, 0) * normalization}
    };
}
};