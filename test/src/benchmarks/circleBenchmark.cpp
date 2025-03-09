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
                std::srand(0);
                const Points<Eigen::Dynamic, PointCircle> m_posCircle {
                    m_circle.radius() * (1 + Eigen::ArrayXd::Random(state.range(0)) * 0.95) ,
                    -Eigen::ArrayXd::Random(state.range(0)).abs() * M_PI };

                m_posFrenet = m_transform.posFrenet(m_posCircle);
                m_posCartes = m_transform.posCartes(m_posCircle);

                // non-uniform length discretization
                Eigen::ArrayXd m_lengths { Eigen::ArrayXd::Zero(state.range(1)) };
                m_lengths(Eigen::seqN(1,state.range(1) - 2)) = Eigen::ArrayXd::Random(state.range(1) - 2).abs() * 3 / 2 * M_PI * m_circle.radius();
                m_lengths(state.range(1) - 1) = 3 / 2 * M_PI * m_circle.radius();
                std::sort(m_lengths.begin(), m_lengths.end());

                m_circlePoly = Polyline { m_circle(m_lengths) };
                m_circleTransform = Transform { std::make_shared<Polyline<Eigen::Dynamic>>(m_circlePoly) };
            }

            template <typename Tpoint>
            void reportError(const Points<Eigen::Dynamic, Tpoint> diff, benchmark::State& state)
            {
                Eigen::ArrayXd errX { diff.x().abs() };
                std::sort(errX.begin(), errX.end());
                state.counters["ErrMaxX"] = errX(errX.rows() - 1);
                state.counters["ErrMedX"] = errX(errX.rows() / 2);
                Eigen::ArrayXd errY { diff.y().abs() };
                std::sort(errY.begin(), errY.end());
                state.counters["ErrMaxY"] = errY(errY.rows() - 1);
                state.counters["ErrMedY"] = errY(errY.rows() / 2);
            }

            const Circle m_circle { 5.0, {0.0, 0.0}, -M_PI };
            const TransformCircle m_transform { std::make_shared<Circle>(m_circle) };

            Points<Eigen::Dynamic, PointFrenet> m_posFrenet {};
            Points<Eigen::Dynamic, PointCartes> m_posCartes {};

            Polyline<Eigen::Dynamic> m_circlePoly {};
            Transform m_circleTransform {};
        };

        BENCHMARK_DEFINE_F(PolylineBenchmark, PosFrenetCartes)(benchmark::State& state)
        {
            Points<Eigen::Dynamic, PointCartes> posCartes {};
            for(auto _ : state)
                posCartes = m_circleTransform.posCartes(m_posFrenet);
            reportError(posCartes - m_posCartes, state);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, PosFrenetCartes)->Ranges({{8, 8<<10}, {8, 8<<10}})->ArgNames({"NumQueries", "NumPoints"});

        BENCHMARK_DEFINE_F(PolylineBenchmark, PosCartesFrenet)(benchmark::State& state)
        {
            Points<Eigen::Dynamic, PointFrenet> posFrenet {};
            for(auto _ : state)
                posFrenet = m_circleTransform.posFrenet(m_posCartes);
            reportError(posFrenet - m_posFrenet, state);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, PosCartesFrenet)->Ranges({{8, 8<<10}, {8, 8<<10}})->ArgNames({"NumQueries", "NumPoints"});
    };
};

BENCHMARK_MAIN();