#include <matplot/matplot.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

#include "frenetTransform/transform.h"
#include "frenetTransform/polychain.h"
#include "circle.h"

namespace Testing = FrenetTransform::Testing;

/**
 * @brief Helper function to determine number of elements in 2-dimensional vector.
 *
 * @param vector2d to determine the element number for.
 * @return constexpr int element number of "vector2d".
 */
constexpr int size(const matplot::vector_2d& vector2d)
{
    int size {};
    // add size of each vector
    for(const auto& vector : vector2d)
        size += vector.size();

    return size;
}

/**
 * @brief Helper function to make 1-dimensional copy of 2-dimensional vector.
 *
 * @param vector2d to reduce to a 1-dimensional vector.
 * @return constexpr matplot::vector_1d 1-dimensional copy of "vector2d".
 */
constexpr matplot::vector_1d ravel(const matplot::vector_2d& vector2d)
{
    // instantiate empty vector with size to hold all elements
    matplot::vector_1d vector1d (size(vector2d));

    // copy values to vector
    int idx {};
    for(const auto& vector : vector2d)
        for(const double& val : vector)
            vector1d[idx++] = val;

    return vector1d;
}

/**
 * @brief Helper function to to make an Eigen::ArrayXd copy of a standard vector.
 *
 * @param vector to copy into an Eigen::ArrayXd.
 * @return Eigen::ArrayXd copy of "vector".
 */
Eigen::ArrayXd toArray(const matplot::vector_1d& vector)
{
    // instantiate new array
    Eigen::ArrayXd array (vector.size());

    // copy values to array
    int idx {};
    for(const double& val : vector)
        array[idx++] = val;

    return array;
}

int main(int argc, char* argv[])
{
    // create circle with radius 10 m
    const Testing::Circle circle { 10.0, {0.0, 0.0} };

    // get points around circle
    const Eigen::ArrayXd lengthsCircle { Eigen::ArrayXd::LinSpaced(101, 0.0, 2 * M_PI) * circle.radius() };
    const FrenetTransform::Points circlePoints { circle(lengthsCircle) };

    // create point grit from -15 to 15 in x- and y-direction
    const double bound { 15.0 };
    auto [posMeshX, posMeshY] { matplot::meshgrid(matplot::iota(0.5 - bound, 1, 0.5 + bound), matplot::iota(0.5 - bound, 1, 0.5 + bound)) };

    // copy mesh to Cartesian points, remove (0, 0)
    FrenetTransform::Points<Eigen::Dynamic> cartesPoints { toArray(ravel(posMeshX)), toArray(ravel(posMeshY))};

    // instantiate cartesian velocities
    const FrenetTransform::Points<Eigen::Dynamic> cartesVels { Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 2, Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 2 };

    // instantiate cartesian acceleratoins
    const FrenetTransform::Points<Eigen::Dynamic> cartesAccs { 3 * Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 4, -3 * Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 4 };

    // plot velocity vector field in Cartesian coordinates
    matplot::subplot(2, 2, 0);

    // plot circle
    matplot::plot(
        std::vector<double> { circlePoints.x().begin(), circlePoints.x().end() },
        std::vector<double> { circlePoints.y().begin(), circlePoints.y().end() }
    );
    matplot::hold(true);
    matplot::quiver(
        std::vector<double> { cartesPoints.x().begin(), cartesPoints.x().end() },
        std::vector<double> { cartesPoints.y().begin(), cartesPoints.y().end() },
        std::vector<double> { cartesVels.x().begin(), cartesVels.x().end() },
        std::vector<double> { cartesVels.y().begin(), cartesVels.y().end() },
        0.0
    );
    matplot::xlim({-bound, 1.5 + bound});
    matplot::ylim({-bound, 1.5 + bound});
    matplot::xlabel("Cartesian x-axis/m");
    matplot::ylabel("Cartesian y-axis/m");
    matplot::hold(false);

    // plot acceleration vector field in Cartesian coordinates
    matplot::subplot(2, 2, 1);
    // plot circle
    matplot::plot(
        std::vector<double> { circlePoints.x().begin(), circlePoints.x().end() },
        std::vector<double> { circlePoints.y().begin(), circlePoints.y().end() }
    );
    matplot::hold(true);
    matplot::quiver(
        std::vector<double> { cartesPoints.x().begin(), cartesPoints.x().end() },
        std::vector<double> { cartesPoints.y().begin(), cartesPoints.y().end() },
        std::vector<double> { cartesAccs.x().begin(), cartesAccs.x().end() },
        std::vector<double> { cartesAccs.y().begin(), cartesAccs.y().end() },
        0.0
    );
    matplot::xlim({-bound, 1.5 + bound});
    matplot::ylim({-bound, 1.5 + bound});
    matplot::xlabel("Cartesian x-axis/m");
    matplot::ylabel("Cartesian y-axis/m");
    matplot::hold(false);

    // generate polyline along circle
    const auto circlePoly { std::make_shared<FrenetTransform::Polychain<Eigen::Dynamic>>(circlePoints.x(), circlePoints.y()) };

    // instantiate transform
    FrenetTransform::Transform<Eigen::Dynamic> transform { circlePoly };

    // transform query points to Frenet frame
    const FrenetTransform::Points<Eigen::Dynamic> frenetPoints { transform.posFrenet(cartesPoints) };
    const FrenetTransform::Points<Eigen::Dynamic> frenetVels { transform.velFrenet(cartesVels, frenetPoints) };
    const FrenetTransform::Points<Eigen::Dynamic> frenetAccs { transform.accFrenet(cartesAccs, frenetVels, frenetPoints) };

    // draw velocity vector field Frenet coordinates
    matplot::subplot(2, 2, 2);
    matplot::quiver(
        std::vector<double> { frenetPoints.x().begin(), frenetPoints.x().end() },
        std::vector<double> { frenetPoints.y().begin(), frenetPoints.y().end() },
        std::vector<double> { frenetVels.x().begin(), frenetVels.x().end() },
        std::vector<double> { frenetVels.y().begin(), frenetVels.y().end() },
        0.0
    );
    matplot::xlim({-1.0, 2 * M_PI * circle.radius()});
    matplot::ylim({-std::sqrt(2) * bound + circle.radius(), std::sqrt(2) * bound - circle.radius()});
    matplot::xlabel("Frenet x-axis/m");
    matplot::ylabel("Frenet y-axis/m");

    // draw acceleration field Frenet coordinates
    matplot::subplot(2, 2, 3);
    matplot::quiver(
        std::vector<double> { frenetPoints.x().begin(), frenetPoints.x().end() },
        std::vector<double> { frenetPoints.y().begin(), frenetPoints.y().end() },
        std::vector<double> { frenetAccs.x().begin(), frenetAccs.x().end() },
        std::vector<double> { frenetAccs.y().begin(), frenetAccs.y().end() },
        0.0
    );
    matplot::xlim({-1.0, 2 * M_PI * circle.radius()});
    matplot::ylim({-std::sqrt(2) * bound + circle.radius(), std::sqrt(2) * bound - circle.radius()});
    matplot::xlabel("Frenet x-axis/m");
    matplot::ylabel("Frenet y-axis/m");

    matplot::save(*(argv + 1));
    // matplot::show();

    return 0;
}