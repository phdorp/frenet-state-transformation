#ifndef POLYLINE_H
#define POLYLINE_H

#include "path.h"
#include "math.h"

namespace FrenetTransform
{
    /**
     * @brief Path representation as polyline.
     *
     * Represents a 2-dimensional path as a polyline.
     * Provide path properties based on finite differences at query points.
     *
     * @tparam T determines the number points along the path.
     */
    template <int T>
    class Polyline : public Path
    {
    public:
        using MatrixT1 = Eigen::Matrix<double, T, 1>;
        using MatrixT2 = Eigen::Matrix<double, T, 2>;

        using ArrayT1 = Eigen::Array<double, T, 1>;

        Polyline() = default;

        /**
         * @brief Construct a new Path object from Cartesian x- and y-positions.
         *
         * @param x coordinates in x-direction.
         * @param y coordinates in y-direction.
         */
        Polyline(const ArrayT1& x, const ArrayT1& y)
            : m_x { x }
            , m_y { y }
            , m_lengths { FrenetTransform::partialLength(x, y) }
        {
        }

        Points operator()(const Eigen::ArrayXd& lengths) const override
        {
            // indices of corresponding polyline segments
            const Eigen::ArrayXi indicesLengths { indices(lengths) };

            // relative position along the linear segment
            const auto pathLengthsdiff { FrenetTransform::diffBackward(m_lengths) };
            const auto relativePos { (lengths - m_lengths(indicesLengths)) / pathLengthsdiff(indicesLengths + 1) };

            // absolute position along path
            const auto x { m_x(indicesLengths + 1) * relativePos + m_x(indicesLengths) * (1 - relativePos) };
            const auto y { m_y(indicesLengths + 1) * relativePos + m_y(indicesLengths) * (1 - relativePos) };

            return { x, y };
        }

        void setPoints(const ArrayT1& x, const ArrayT1& y)
        {
            m_x = x;
            m_y = y;
            m_lengths = FrenetTransform::partialLength(x, y);
        }

    private:
        ArrayT1 m_x {}; /*<< coordinates in x-direction*/
        ArrayT1 m_y {}; /*<< coordinates in y-direction*/
        ArrayT1 m_lengths {}; /*<< partial lengths along polyline*/

        /**
         * @brief Stores gradients from order 1 to 3.
         *
         */
        std::array<const MatrixT2, 3> m_gradients {};

        /**
         * @brief Determines 1st order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 1st order gradient at given path lengths.
         */
        // Eigen::MatrixX2d gradient1 (const Eigen::MatrixXd& lengths) override;

        /**
         * @brief Determines 2nd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 2nd order gradient at given path lengths.
         */
        // Eigen::MatrixX2d gradient2 (const Eigen::MatrixXd& lengths) override;

        /**
         * @brief Determines 3rd order gradient at the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return 3rd order gradient at given path lengths.
         */
        // Eigen::MatrixX2d gradient3 (const Eigen::MatrixXd& lengths) override;

        /**
         * @brief Determines indices of polyline segment corresponding to the given path lengths.
         *
         * @param lengths lengths along the path.
         * @return indices corresponding to given path lengths.
         */
        Eigen::ArrayXi indices(const Eigen::ArrayXd& lengths) const
        {
            Eigen::ArrayXi result(lengths.size()); // vector of segment indices

            // get indices of next segments
            for(int row {}; row < lengths.rows(); ++row)
                result(row) = FrenetTransform::first(m_lengths - lengths(row));

            return result;
        }
    };
};

#endif