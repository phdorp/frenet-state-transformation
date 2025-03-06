#include <benchmark/benchmark.h>

#include "point.h"
#include "points.h"
#include "polyline.h"
#include "circle.h"
#include "transformCircle.h"

namespace FrenetTransform
{
    namespace Testing
    {
        class PolylineBenchmark : public benchmark::Fixture
        {
        protected:
            void SetUp(::benchmark::State& state) {
                m_circlePoly = Polyline { m_circle(Eigen::ArrayXd::LinSpaced(state.range(0), 0.0, 2 * M_PI) * m_circle.radius()) };
                m_circleTransform = Transform { std::make_shared<Polyline<Eigen::Dynamic>>(m_circlePoly) };
            }

            const Circle m_circle { 5.0, {0.0, 0.0}, -M_PI };
            const TransformCircle m_transform { std::make_shared<Circle>(m_circle) };

            const Points<Eigen::Dynamic, PointCircle> m_posCircle { m_circle.radius() * (1 + Eigen::ArrayXd::Random(100) * 0.95) , Eigen::ArrayXd::Random(100) * M_PI * 0.99 };
            const Points<Eigen::Dynamic, PointFrenet> m_posFrenet { m_transform.posFrenet(m_posCircle) };
            const Points<Eigen::Dynamic, PointCartes> m_posCartes { m_transform.posCartes(m_posCircle) };

            Polyline<Eigen::Dynamic> m_circlePoly;
            Transform m_circleTransform;
        };

        BENCHMARK_DEFINE_F(PolylineBenchmark, PosFrenetCartes)(benchmark::State& state)
        {
            for(auto _ : state)
                m_circleTransform.posCartes(m_posFrenet);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, PosFrenetCartes)->Range(8, 8<<10)->ArgName("NumPoints");
    };
};

BENCHMARK_MAIN();