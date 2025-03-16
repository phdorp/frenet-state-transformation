#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <Eigen/Core>
#include <memory>

#include "frenetTransform/path.h"
#include "frenetTransform/points.h"
#include "frenetTransform/math.h"

namespace FrenetTransform
{
    /**
     * @brief Transformation between Cartesian and Frenet frame
     *
     * Use the Path properties to implement the transformation independent from the Path implementation.
     */
    template <int NumQueries=Eigen::Dynamic>
    class Transform
    {
    public:
        using ArrayQueries = Eigen::Array<double, NumQueries, 1>;

        Transform() = default;

        /**
         * @brief Construct a new Transform object from a given Path.
         *
         * @param path defines the Frenet frame.
         */
        Transform(const std::shared_ptr<Path<NumQueries>> path) : m_path {path}
        {
        }

        Points<NumQueries> posFrenet(const Points<NumQueries>& posCartes) const
        {
            const auto lengths { m_path->lengths(posCartes) };
            const auto posPath { m_path->operator()(lengths) };
            const auto posDiff { posCartes - posPath };
            const auto normals { m_path->normal(lengths) };
            return { lengths, normals * posDiff };
        }

        Points<NumQueries> posCartes(const Points<NumQueries>& posFrenet) const
        {
            const auto posPath { m_path->operator()(posFrenet.x()) };
            const auto normals { m_path->normal(posFrenet.x()) };
            return posPath + normals * posFrenet.y();
        }

        Points<NumQueries> velFrenet(const Points<NumQueries>& velCartes, const Points<NumQueries>& posFrenet) const
        {
            const auto velTransformsInv { transformInv(velTransform(posFrenet)) };
            return {
                velTransformsInv(0, 0) * velCartes.x() + velTransformsInv(0, 1) * velCartes.y(),
                velTransformsInv(1, 0) * velCartes.x() + velTransformsInv(1, 1) * velCartes.y()
            };
        }

        Points<NumQueries> velCartes(const Points<NumQueries>& velFrenet, const Points<NumQueries>& posFrenet) const
        {
            const auto velTransforms { velTransform(posFrenet) };
            return {
                velTransforms(0, 0) * velFrenet.x() + velTransforms(0, 1) * velFrenet.y(),
                velTransforms(1, 0) * velFrenet.x() + velTransforms(1, 1) * velFrenet.y()
            };
        }

        Points<NumQueries> accFrenet(const Points<NumQueries>& accCartes, const Points<NumQueries>& velFrenet, const Points<NumQueries>& posFrenet) const
        {
            const auto accTransforms { accTransform(velFrenet, posFrenet) };
            const Points<NumQueries> accRel {
                accCartes.x() - accTransforms(0, 0) * velFrenet.x() - accTransforms(0, 1) * velFrenet.y(),
                accCartes.y() - accTransforms(1, 0) * velFrenet.x() - accTransforms(1, 1) * velFrenet.y()
            };

            const auto velTransformsInv { transformInv(velTransform(posFrenet)) };

            return {
                velTransformsInv(0, 0) * accRel.x() + velTransformsInv(0, 1) * accRel.y(),
                velTransformsInv(1, 0) * accRel.x() + velTransformsInv(1, 1) * accRel.y()
            };
        }

        Points<NumQueries> accCartes(const Points<NumQueries>& accFrenet, const Points<NumQueries>& velFrenet, const Points<NumQueries>& posFrenet) const
        {
            const auto accTransforms { accTransform(velFrenet, posFrenet) };
            const auto velTransforms { velTransform(posFrenet) };
            return {
                velTransforms(0, 0) * accFrenet.x() + velTransforms(0, 1) * accFrenet.y() + accTransforms(0, 0) * velFrenet.x() + accTransforms(0, 1) * velFrenet.y(),
                velTransforms(1, 0) * accFrenet.x() + velTransforms(1, 1) * accFrenet.y() + accTransforms(1, 0) * velFrenet.x() + accTransforms(1, 1) * velFrenet.y()
            };
        }

    protected:
        std::shared_ptr<Path<NumQueries>> m_path; /**< Store path. */

    private:
        Eigen::Array<ArrayQueries, 2, 2> velTransform(const Points<NumQueries>& posFrenet) const
        {
            const auto tangents { m_path->tangent(posFrenet.x()) };
            const auto normals { m_path->normal(posFrenet.x()) };
            const auto curvs { m_path->angle1(posFrenet.x()) };
            return {
                {tangents.x() * (1 - curvs * posFrenet.y()), normals.x()},
                {tangents.y() * (1 - curvs * posFrenet.y()), normals.y()}
            };
        }

        Eigen::Array<ArrayQueries, 2, 2> accTransform(const Points<NumQueries>& velFrenet, const Points<NumQueries>& posFrenet) const
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
    };
};

#endif