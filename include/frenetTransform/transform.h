#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <Eigen/Core>
#include <memory>

#include "frenetTransform/path.h"
#include "frenetTransform/points.h"

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

        Points<Eigen::Dynamic> posFrenet(const Points<Eigen::Dynamic>& posCartes) const;

        Points<Eigen::Dynamic> velFrenet(const Points<Eigen::Dynamic>& velCartes, const Points<Eigen::Dynamic>& posFrenet) const;

        Points<Eigen::Dynamic> accFrenet(const Points<Eigen::Dynamic>& accCartes, const Points<Eigen::Dynamic>& velFrenet, const Points<Eigen::Dynamic>& posFrenet) const;

        Points<Eigen::Dynamic> posCartes(const Points<Eigen::Dynamic>& posFrenet) const;

        Points<Eigen::Dynamic> velCartes(const Points<Eigen::Dynamic>& velFrenet, const Points<Eigen::Dynamic>& posFrenet) const;

        Points<Eigen::Dynamic> accCartes(const Points<Eigen::Dynamic>& accFrenet, const Points<Eigen::Dynamic>& velFrenet, const Points<Eigen::Dynamic>& posFrenet) const;

    protected:
        std::shared_ptr<Path> m_path; /**< Store path. */

    private:
        Eigen::Array<Eigen::ArrayXd, 2, 2> velTransform(const Points<Eigen::Dynamic>& posFrenet) const;

        Eigen::Array<Eigen::ArrayXd, 2, 2> accTransform(const Points<Eigen::Dynamic>& velFrenet, const Points<Eigen::Dynamic>& posFrenet) const;

        static Eigen::Array<Eigen::ArrayXd, 2, 2> transformInv(const Eigen::Array<Eigen::ArrayXd, 2, 2>& transform);
    };
};

#endif