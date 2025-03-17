#include <matplot/matplot.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

#include "frenetTransform/transform.h"
#include "frenetTransform/polychain.h"
#include "circle.h"

namespace Testing = FrenetTransform::Testing;

int main(int argc, char* argv[])
{
    assert(argc > 1 && "Only one argument possible!");

    // create circle with radius 10 m
    const Testing::Circle circle { 10.0, {0.0, 0.0} };

    // get points around circle
    const Eigen::ArrayXd lengthsCircle { Eigen::ArrayXd::LinSpaced(101, 0.0, 2 * M_PI) * circle.radius() };
    const FrenetTransform::Points circlePoints { circle(lengthsCircle) };

    // plot circle
    matplot::plot(
        std::vector<double> { circlePoints.x().begin(), circlePoints.x().end() },
        std::vector<double> { circlePoints.y().begin(), circlePoints.y().end() }
    );
    matplot::hold(true);

    // generate polyline from 4 points at 90° angles
    const auto idcsSub { Eigen::seqN(0, 5, 25) };
    const auto circlePoly {
        std::make_shared<FrenetTransform::Polychain<Eigen::Dynamic>>(circlePoints.x()(idcsSub), circlePoints.y()(idcsSub))
    };

    // get points along the polychain
    const Eigen::ArrayXd lengthsPoly { Eigen::ArrayXd::LinSpaced(500, 0.0, 2 * M_PI * circle.radius()) };
    const auto polyPoints { circlePoly->operator()(lengthsPoly) };

    // plot polychain
    matplot::plot(
        std::vector<double> { polyPoints.x().begin(), polyPoints.x().end() },
        std::vector<double> { polyPoints.y().begin(), polyPoints.y().end() }
    );

    // instantiate transform
    FrenetTransform::Transform<Eigen::Dynamic> transform { circlePoly };

    // setup query points in Cartesian frame
    const FrenetTransform::Points<Eigen::Dynamic> cartesPoints {
        {{ 5.0, 12.0, -2.5,  0.3}},
        {{ 0.0, 12.0,  3.0, 11.5}}
    };

    // plot query points
    matplot::scatter(
        std::vector<double> { cartesPoints.x().begin(), cartesPoints.x().end() },
        std::vector<double> { cartesPoints.y().begin(), cartesPoints.y().end() }
    );

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

    // draw vectors
    matplot::quiver(
        std::vector<double> { projPoints.x().begin(), projPoints.x().end() },
        std::vector<double> { projPoints.y().begin(), projPoints.y().end() },
        std::vector<double> { normals.x().begin(), normals.x().end() },
        std::vector<double> { normals.y().begin(), normals.y().end() },
        0.0
    );

    // show plot if no path provided, other save
    switch(argc)
    {
    case 1:
        matplot::show();
    default:
        matplot::save(*(argv + 1));
    }

    return 0;
}