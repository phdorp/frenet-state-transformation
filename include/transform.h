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
        Transform() = delete;

        /**
         * @brief Construct a new Transform object from a given Path.
         *
         * @param path defines the Frenet frame.
         */
        Transform(const std::shared_ptr<Path> path) : m_path {path}
        {
        }

        Points<Eigen::Dynamic> posFrenet(const Points<Eigen::Dynamic>& posCartes) const;

        Points<Eigen::Dynamic> velFrenet(const Points<Eigen::Dynamic>& velCartes, const Points<Eigen::Dynamic>& posFrenet) const;

        Points<Eigen::Dynamic> accFrenet(const Points<Eigen::Dynamic>& accCartes, const Points<Eigen::Dynamic>& velFrenet, const Points<Eigen::Dynamic>& posFrenet) const;

        Points<Eigen::Dynamic> posCartes(const Points<Eigen::Dynamic>& posFrenet) const;

        Points<Eigen::Dynamic> velCartes(const Points<Eigen::Dynamic>& velFrenet, const Points<Eigen::Dynamic>& posCartes) const;

        Points<Eigen::Dynamic> accCartes(const Points<Eigen::Dynamic>& accFrenet, const Points<Eigen::Dynamic>& velCartes, const Points<Eigen::Dynamic>& posCartes) const;

    private:
        const std::shared_ptr<Path> m_path; /**< Store path. */
    };
};

#endif