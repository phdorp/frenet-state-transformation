#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <eigen3/Eigen/Core>
#include <memory>

#include "path.h"
#include "points.h"

using FrenetTransform::Points;

namespace FrenetTransform
{
    /**
     * @brief Transformation between Cartesian and Frenet frame
     *
     * Use the Path properties to implement the transformation independent from the Path implementation.
     */
    class Transform
    {
    public:
        Transform() = default;

        /**
         * @brief Construct a new Transform object from a given Path.
         *
         * @param path defines the Frenet frame.
         */
        Transform(const std::shared_ptr<Path> path) : m_path {path}
        {
        }

        Points<Eigen::Dynamic, PointFrenet> posFrenet(const Points<Eigen::Dynamic, PointCartes>& posCartes) const;

        Points<Eigen::Dynamic, PointFrenet> velFrenet(const Points<Eigen::Dynamic, PointCartes>& velCartes, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const;

        Points<Eigen::Dynamic, PointFrenet> accFrenet(const Points<Eigen::Dynamic, PointCartes>& accCartes, const Points<Eigen::Dynamic, PointFrenet>& velFrenet, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const;

        Points<Eigen::Dynamic, PointCartes> posCartes(const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const;

        Points<Eigen::Dynamic, PointCartes> velCartes(const Points<Eigen::Dynamic, PointFrenet>& velFrenet, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const;

        Points<Eigen::Dynamic, PointCartes> accCartes(const Points<Eigen::Dynamic, PointFrenet>& accFrenet, const Points<Eigen::Dynamic, PointFrenet>& velFrenet, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const;

    protected:
        std::shared_ptr<Path> m_path; /**< Store path. */

    private:
        Eigen::Array<Eigen::ArrayXd, 2, 2> velTransform(const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const;

        Eigen::Array<Eigen::ArrayXd, 2, 2> accTransform(const Points<Eigen::Dynamic, PointFrenet>& velFrenet, const Points<Eigen::Dynamic, PointFrenet>& posFrenet) const;

        static Eigen::Array<Eigen::ArrayXd, 2, 2> transformInv(const Eigen::Array<Eigen::ArrayXd, 2, 2>& transform);
    };
};

#endif