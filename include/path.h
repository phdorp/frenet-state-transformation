#ifndef PATH_H
#define PATH_H

#include <eigen3/Eigen/Core>
#include <array>

namespace FrenetTransform
{
    /**
     * @brief Path representation as polyline.
     *
     * Represents a 2-dimensional path as a polyline.
     * Provide path properties based on finite differences at query points.
     * The properties include orientation, curvature, cuvature change, normal and tangential.
     *
     * @tparam T determines the number points along the path.
     */
    template <int T>
    class Path
    {
    public:
        using MatrixT2 = Eigen::Matrix<double, T, 2>;

        Path() = delete;

        /**
         * @brief Construct a new Path object from Cartesian x- and y-positions.
         *
         * @param points points along the path in Cartesian frame.
         */
        Path(const MatrixT2& points);

        /**
         * @brief Determines points at the given path lengths.
         *
         * @tparam N number of query lengths.
         * @param lengths lengths along the path.
         * @return points at given path lengths.
         */
        template <int N>
        Eigen::Matrix<double, N, 2> operator()(const Eigen::Matrix<double, N, 1>& lengths);

        /**
         * @brief Determines tangent vectors at the given path lengths.
         *
         * @tparam N number of query lengths.
         * @param lengths lengths along the path.
         * @return tangent vectors.
         */
        template <int N>
        Eigen::Matrix<double, N, 2> tangent(const Eigen::Matrix<double, N, 1>& lengths);

        /**
         * @brief Determines normal vectors at the given path lengths.
         *
         * @tparam N number of query lengths.
         * @param lengths lengths along the path.
         * @return normal vectors.
         */
        template <int N>
        Eigen::Matrix<double, N, 2> normal(const Eigen::Matrix<double, N, 1>& lengths);

        /**
         * @brief Determines path angle at the given path lengths.
         *
         * @tparam N number of query lengths.
         * @param lengths lengths along the path.
         * @return path angles.
         */
        template <int N>
        Eigen::Matrix<double, N, 1> angle0(const Eigen::Matrix<double, N, 1>& lengths);

        /**
         * @brief Determines path curvature at the given path lengths.
         *
         * @tparam N number of query lengths.
         * @param lengths lengths along the path.
         * @return path curvatures.
         */
        template <int N>
        Eigen::Matrix<double, N, 1> angle1(const Eigen::Matrix<double, N, 1>& lengths);

        /**
         * @brief Determines path curvature derivative at the given path lengths.
         *
         * @tparam N number of query lengths.
         * @param lengths lengths along the path.
         * @return path curvatures.
         */
        template <int N>
        Eigen::Matrix<double, N, 1> angle2(const Eigen::Matrix<double, N, 1>& lengths);

private:
        /**
         * @brief Stores gradients from order 1 to 3.
         *
         */
        std::array<const MatrixT2>, 3> m_gradients {};
    };
};

#endif