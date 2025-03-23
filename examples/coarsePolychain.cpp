#include <matplot/matplot.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

#include "frenetTransform/transform.h"
#include "frenetTransform/polychain.h"
#include "frenetTransform/internal/circle.h"

namespace Internal = FrenetTransform::Internal;

int main(int argc, char* argv[])
{
    assert(argc > 1 && "Only one argument possible!");

    // create circle with radius 10 m
    const double radius { 10.0 };
    const Eigen::ArrayXd lengthsCircle { Eigen::ArrayXd::LinSpaced(101, 0.0, 2 * M_PI) };
    const Eigen::ArrayXd circlePointsX { radius * lengthsCircle.cos() };
    const Eigen::ArrayXd circlePointsY { radius * lengthsCircle.sin() };

    // generate polyline from 4 points at 90° angles
    const auto idcsSub { Eigen::seqN(0, 5, 25) };
    const auto circlePoly {
        std::make_shared<FrenetTransform::Polychain<Eigen::Dynamic>>(circlePointsX(idcsSub), circlePointsY(idcsSub))
    };

    // get points along the polychain
    const Eigen::ArrayXd lengthsPoly { Eigen::ArrayXd::LinSpaced(500, 0.0, 2 * M_PI * radius) };
    const auto polyPoints { circlePoly->operator()(lengthsPoly) };

    // instantiate transform
    FrenetTransform::Transform<Eigen::Dynamic> transform { circlePoly };

    // setup query points in Cartesian frame
    const FrenetTransform::Points<Eigen::Dynamic> cartesPoints {
        {{ 5.0, 12.0, -2.5,  0.3}},
        {{ 0.0, 12.0,  3.0, 11.5}}
    };

    // transform query points to Frenet frame
    const FrenetTransform::Points<Eigen::Dynamic> frenetPoints { transform.posFrenet(cartesPoints) };

    // get next points along polyline to query points
    const FrenetTransform::Points<Eigen::Dynamic> projPoints { circlePoly->operator()(frenetPoints.x()) };

    // plot next points
    matplot::scatter(
        std::vector<double> { projPoints.x().begin(), projPoints.x().end() },
        std::vector<double> { projPoints.y().begin(), projPoints.y().end() }
    );

    // get vectors from next points to query points
    FrenetTransform::Points<Eigen::Dynamic> normals { circlePoly->normal(frenetPoints.x()) };
    normals = normals * frenetPoints.y();

    // plot circle
    matplot::plot(
        std::vector<double> { circlePointsX.begin(), circlePointsX.end() },
        std::vector<double> { circlePointsY.begin(), circlePointsY.end() }
    );
    matplot::hold(true);

    // plot polychain
    matplot::plot(
        std::vector<double> { polyPoints.x().begin(), polyPoints.x().end() },
        std::vector<double> { polyPoints.y().begin(), polyPoints.y().end() }
    );

    // plot query points
    matplot::scatter(
        std::vector<double> { cartesPoints.x().begin(), cartesPoints.x().end() },
        std::vector<double> { cartesPoints.y().begin(), cartesPoints.y().end() }
    );

    // draw vectors
    matplot::quiver(
        std::vector<double> { projPoints.x().begin(), projPoints.x().end() },
        std::vector<double> { projPoints.y().begin(), projPoints.y().end() },
        std::vector<double> { normals.x().begin(), normals.x().end() },
        std::vector<double> { normals.y().begin(), normals.y().end() },
        0.0
    );

    matplot::save(*(argv + 1));
    matplot::show();

    return 0;
}