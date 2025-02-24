#ifndef PATH_H
#define PATH_H

#include <eigen3/Eigen/Core>
#include <array>
#include <math.h>

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
        Eigen::ArrayXd angle0(const Eigen::ArrayXd& lengths) const
        {
            const Points tangents { tangent(lengths) };

            Eigen::ArrayXd result (lengths.size());
            for(int iLength {}; iLength < lengths.size(); ++iLength)
            {
                if(std::abs(tangents.x(iLength)) < 0.5)
                    result(iLength) = M_PI / 2 - std::atan(tangents.x(iLength) / tangents.y(iLength));
                else
                    result(iLength) = std::atan(tangents.y(iLength) / tangents.x(iLength));

                if(tangents.x(iLength) < 0 && tangents.y(iLength) > 0)
                    result(iLength) += M_PI;
                else if(tangents.x(iLength) < 0 && tangents.y(iLength) < 0)
                    result(iLength) -= M_PI;
            }

            return result;
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
            const auto grad1Abs { (grad1.x.square() + grad1.y.square()).sqrt() };
            return -(grad1.y * grad2.x - grad1.x * grad2.y) / grad1Abs.pow(3);
        }

        /**
         * @brief Determines path curvature derivative at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return path curvatures.
         */
        Eigen::MatrixXd angle2(const Eigen::MatrixXd& lengths) const;

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