#ifndef PATH_H
#define PATH_H

#include <eigen3/Eigen/Core>
#include <array>

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

        struct Points
        {
            Eigen::ArrayXd x {};
            Eigen::ArrayXd y {};
        };

        /**
         * @brief Determines points at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return points at given path lengths.
         */
        virtual Points operator()(const Eigen::ArrayXd& lengths) const = 0;

        /**
         * @brief Determines tangent vectors at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return tangent vectors.
         */
        Points tangent(const Eigen::ArrayXd& lengths) const { return gradient1(lengths); }

        /**
         * @brief Determines normal vectors at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return normal vectors.
         */
        Points normal(const Eigen::ArrayXd& lengths) const
        {
            const Points tangents { tangent(lengths) };
            return { tangents.y, -tangents.x };
        }

        /**
         * @brief Determines path angle at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return path angles.
         */
        Eigen::MatrixXd angle0(const Eigen::MatrixXd& lengths);

        /**
         * @brief Determines path curvature at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return path curvatures.
         */
        Eigen::MatrixXd angle1(const Eigen::MatrixXd& lengths);

        /**
         * @brief Determines path curvature derivative at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return path curvatures.
         */
        Eigen::MatrixXd angle2(const Eigen::MatrixXd& lengths);

    private:
        /**
         * @brief Determines 1st order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 1st order gradient at given path lengths.
         */
        virtual Points gradient1(const Eigen::MatrixXd& lengths) const = 0;

        /**
         * @brief Determines 2nd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 2nd order gradient at given path lengths.
         */
        virtual Points gradient2(const Eigen::MatrixXd& lengths) const = 0;

        /**
         * @brief Determines 3rd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 3rd order gradient at given path lengths.
         */
        virtual Points gradient3(const Eigen::MatrixXd& lengths) const = 0;
    };
};

#endif