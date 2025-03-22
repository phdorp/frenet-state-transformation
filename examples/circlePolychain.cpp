#include <matplot/matplot.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

#include "frenetTransform/transform.h"
#include "frenetTransform/polychain.h"
#include "frenetTransform/internal/circle.h"

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
    const double radius { 10.0 };
    const Eigen::ArrayXd lengthsCircle { Eigen::ArrayXd::LinSpaced(101, 0.0, 2 * M_PI) };
    const Eigen::ArrayXd circlePointsX { radius * lengthsCircle.cos() };
    const Eigen::ArrayXd circlePointsY { radius * lengthsCircle.sin() };

    // create point grid from -15 to 15 in x- and y-direction
    const double bound { 15.0 };
    auto [posMeshX, posMeshY] { matplot::meshgrid(matplot::iota(0.5 - bound, 1, 0.5 + bound), matplot::iota(0.5 - bound, 1, 0.5 + bound)) };

    // copy mesh to Cartesian points, remove (0, 0)
    FrenetTransform::Points<Eigen::Dynamic> cartesPoints { toArray(ravel(posMeshX)), toArray(ravel(posMeshY))};

    // instantiate cartesian velocities and accelerations
    const FrenetTransform::Points<Eigen::Dynamic> cartesVels { Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 2, Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 2 };
    const FrenetTransform::Points<Eigen::Dynamic> cartesAccs { 3 * Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 4, -3 * Eigen::ArrayXd::Ones(cartesPoints.numPoints()) / 4 };

    // generate polyline along circle
    const auto circlePoly { std::make_shared<FrenetTransform::Polychain<Eigen::Dynamic>>(circlePointsX, circlePointsY) };

    // instantiate transform
    FrenetTransform::Transform<Eigen::Dynamic> transform { circlePoly };

    // transform query points to Frenet frame
    const FrenetTransform::Points<Eigen::Dynamic> frenetPointsTf { transform.posFrenet(cartesPoints) };
    const FrenetTransform::Points<Eigen::Dynamic> frenetVelsTf { transform.velFrenet(cartesVels, frenetPointsTf) };
    const FrenetTransform::Points<Eigen::Dynamic> frenetAccsTf { transform.accFrenet(cartesAccs, frenetVelsTf, frenetPointsTf) };

    // transform query points back to Cartesian frame
    const FrenetTransform::Points<Eigen::Dynamic> cartesPointsTf { transform.posCartes(frenetPointsTf) };
    const FrenetTransform::Points<Eigen::Dynamic> cartesVelsTf { transform.velCartes(frenetVelsTf, frenetPointsTf) };
    const FrenetTransform::Points<Eigen::Dynamic> cartesAccsTf { transform.accCartes(frenetAccsTf, frenetVelsTf, frenetPointsTf) };

    // plot velocity vector field in Cartesian coordinates
    matplot::subplot(3, 2, 0);
    // plot circle
    matplot::plot(
        std::vector<double> { circlePointsX.begin(), circlePointsX.end() },
        std::vector<double> { circlePointsY.begin(), circlePointsY.end() }
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
    matplot::subplot(3, 2, 1);
    // plot circle
    matplot::plot(
        std::vector<double> { circlePointsX.begin(), circlePointsX.end() },
        std::vector<double> { circlePointsY.begin(), circlePointsY.end() }
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

    // draw velocity vector field Frenet coordinates
    matplot::subplot(3, 2, 2);
    matplot::quiver(
        std::vector<double> { frenetPointsTf.x().begin(), frenetPointsTf.x().end() },
        std::vector<double> { frenetPointsTf.y().begin(), frenetPointsTf.y().end() },
        std::vector<double> { frenetVelsTf.x().begin(), frenetVelsTf.x().end() },
        std::vector<double> { frenetVelsTf.y().begin(), frenetVelsTf.y().end() },
        0.0
    );
    matplot::xlim({-1.0, 2 * M_PI * radius});
    matplot::ylim({-std::sqrt(2) * bound + radius, std::sqrt(2) * bound - radius});
    matplot::xlabel("Frenet x-axis/m");
    matplot::ylabel("Frenet y-axis/m");

    // draw acceleration field Frenet coordinates
    matplot::subplot(3, 2, 3);
    matplot::quiver(
        std::vector<double> { frenetPointsTf.x().begin(), frenetPointsTf.x().end() },
        std::vector<double> { frenetPointsTf.y().begin(), frenetPointsTf.y().end() },
        std::vector<double> { frenetAccsTf.x().begin(), frenetAccsTf.x().end() },
        std::vector<double> { frenetAccsTf.y().begin(), frenetAccsTf.y().end() },
        0.0
    );
    matplot::xlim({-1.0, 2 * M_PI * radius});
    matplot::ylim({-std::sqrt(2) * bound + radius, std::sqrt(2) * bound - radius});
    matplot::xlabel("Frenet x-axis/m");
    matplot::ylabel("Frenet y-axis/m");

    // draw velocity vector field Cartesian coordinates
    matplot::subplot(3, 2, 4);
    matplot::quiver(
        std::vector<double> { cartesPointsTf.x().begin(), cartesPointsTf.x().end() },
        std::vector<double> { cartesPointsTf.y().begin(), cartesPointsTf.y().end() },
        std::vector<double> { cartesVelsTf.x().begin(), cartesVelsTf.x().end() },
        std::vector<double> { cartesVelsTf.y().begin(), cartesVelsTf.y().end() },
        0.0
    );
    matplot::xlim({-bound, 1.5 + bound});
    matplot::ylim({-bound, 1.5 + bound});
    matplot::xlabel("Cartesian x-axis/m");
    matplot::ylabel("Cartesian y-axis/m");

    // draw acceleration field Cartesian coordinates
    matplot::subplot(3, 2, 5);
    matplot::quiver(
        std::vector<double> { cartesPointsTf.x().begin(), cartesPointsTf.x().end() },
        std::vector<double> { cartesPointsTf.y().begin(), cartesPointsTf.y().end() },
        std::vector<double> { cartesAccsTf.x().begin(), cartesAccsTf.x().end() },
        std::vector<double> { cartesAccsTf.y().begin(), cartesAccsTf.y().end() },
        0.0
    );
    matplot::xlim({-bound, 1.5 + bound});
    matplot::ylim({-bound, 1.5 + bound});
    matplot::xlabel("Cartesian x-axis/m");
    matplot::ylabel("Cartesian y-axis/m");

    matplot::save(*(argv + 1));
    matplot::show();

    return 0;
}