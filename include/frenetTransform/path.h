#ifndef PATH_H
#define PATH_H

#include <Eigen/Core>
#include <array>
#include <math.h>

#include "frenetTransform/points.h"
#include "frenetTransform/math.h"

using FrenetTransform::Points;

namespace FrenetTransform
{
    /**
     * @brief Path base class.
     *
     * The properties include orientation, curvature, cuvature change, normal and tangential.
     */
    class Path
    {
    public:

        virtual ~Path() = default;

        /**
         * @brief Determines points at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return points at given path lengths.
         */
        virtual Points<Eigen::Dynamic> operator()(const Eigen::ArrayXd& lengths) const = 0;

        /**
         * @brief Determines tangent vectors at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return tangent vectors.
         */
        Points<Eigen::Dynamic> tangent(const Eigen::ArrayXd& lengths) const { return gradient1(lengths); }

        /**
         * @brief Determines normal vectors at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return normal vectors.
         */
        Points<Eigen::Dynamic> normal(const Eigen::ArrayXd& lengths) const
        {
            const Points<Eigen::Dynamic> tangents { tangent(lengths) };
            return { -tangents.y(), tangents.x() };
        }

        /**
         * @brief Determines next points to the query points.
         *
         * @param points query points.
         * @return Points<Eigen::Dynamic> next to query points.
         */
        virtual Eigen::ArrayXd lengths(const Points<Eigen::Dynamic>& points) const = 0;

        /**
         * @brief Determines path angle at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return path angles.
         */
        Eigen::ArrayXd angle0(const Eigen::ArrayXd& lengths) const
        {
            const Points<Eigen::Dynamic> tangents { tangent(lengths) };

            return angleDir(tangents.x(), tangents.y());
        }

        /**
         * @brief Determines path curvature at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return path curvatures.
         */
        Eigen::ArrayXd angle1(const Eigen::ArrayXd& lengths) const
        {
            const auto grad1 { gradient1(lengths) };
            const auto grad2 { gradient2(lengths) };
            const auto grad1Abs { (grad1.x().square() + grad1.y().square()).sqrt() };
            return -(grad1.y() * grad2.x() - grad1.x() * grad2.y()) / grad1Abs.pow(3);
        }

        /**
         * @brief Determines path curvature derivative at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return path curvatures.
         */
        Eigen::ArrayXd angle2(const Eigen::ArrayXd& lengths) const
        {
            const auto grad1 { gradient1(lengths) };
            const auto grad1Abs { (grad1.x().square() + grad1.y().square()).sqrt() };

            const auto grad2 { gradient2(lengths) };
            const auto grad3 { gradient3(lengths) };

            return (grad1.x() * grad3.y() - grad3.x() * grad1.y()) / grad1Abs.pow(3)
                - (grad1.x() * grad2.x() + grad1.y() * grad2.y()) * (grad1.x() * grad2.y() - grad2.x() * grad1.y())
                / grad1Abs.pow(5);
        }

    private:
        /**
         * @brief Determines 1st order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 1st order gradient at given path lengths.
         */
        virtual Points<Eigen::Dynamic> gradient1(const Eigen::ArrayXd& lengths) const = 0;

        /**
         * @brief Determines 2nd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 2nd order gradient at given path lengths.
         */
        virtual Points<Eigen::Dynamic> gradient2(const Eigen::ArrayXd& lengths) const = 0;

        /**
         * @brief Determines 3rd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 3rd order gradient at given path lengths.
         */
        virtual Points<Eigen::Dynamic> gradient3(const Eigen::ArrayXd& lengths) const = 0;
    };
};

#endif