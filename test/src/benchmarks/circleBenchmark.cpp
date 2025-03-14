#include <benchmark/benchmark.h>

#include "frenetTransform/point.h"
#include "frenetTransform/points.h"
#include "frenetTransform/polychain.h"
#include "circle.h"
#include "transformCircle.h"

namespace FrenetTransform
{
    namespace Testing
    {
        class PolylineBenchmark : public benchmark::Fixture
        {
        protected:
            const Circle m_circle { 5.0, {0.0, 0.0}, -M_PI };
            const TransformCircle m_transform { std::make_shared<Circle>(m_circle) };

            Points<Eigen::Dynamic> m_posFrenet {};
            Points<Eigen::Dynamic> m_velFrenet {};
            Points<Eigen::Dynamic> m_accFrenet {};

            Points<Eigen::Dynamic> m_posCartes {};
            Points<Eigen::Dynamic> m_velCartes {};
            Points<Eigen::Dynamic> m_accCartes {};

            Polychain<Eigen::Dynamic> m_circlePoly {};
            Transform m_circleTransform {};

            void SetUp(::benchmark::State& state) {
                std::srand(0);
                const Eigen::ArrayXd posCircleX  { m_circle.radius() * (1 + Eigen::ArrayXd::Random(state.range(0)) * 0.95) };
                std::srand(1);
                const Eigen::ArrayXd posCircleY { -Eigen::ArrayXd::Random(state.range(0)).abs() * M_PI };

                const Points<Eigen::Dynamic, PointCircle> m_posCircle { posCircleX, posCircleY };
                m_posFrenet = m_transform.posFrenet(m_posCircle);
                m_posCartes = m_transform.posCartes(m_posCircle);

                std::srand(2);
                const Eigen::ArrayXd velCircleX  { m_circle.radius() * (1 + Eigen::ArrayXd::Random(state.range(0)) * 0.95) };
                std::srand(3);
                const Eigen::ArrayXd velCircleY { Eigen::ArrayXd::Random(state.range(0)) * M_PI / 4 };

                const Points<Eigen::Dynamic, PointCircle> m_velCircle { velCircleX, velCircleY };
                m_velFrenet = m_transform.velFrenet(m_velCircle);
                m_velCartes = m_transform.velCartes(m_velCircle, m_posCircle);

                std::srand(4);
                const Eigen::ArrayXd accCircleX  { m_circle.radius() * (1 + Eigen::ArrayXd::Random(state.range(0)) * 0.95) };
                std::srand(5);
                const Eigen::ArrayXd accCircleY { Eigen::ArrayXd::Random(state.range(0)) * M_PI / 4 };

                const Points<Eigen::Dynamic, PointCircle> m_accCircle { accCircleX, accCircleY };
                m_accFrenet = m_transform.accFrenet(m_accCircle);
                m_accCartes = m_transform.accCartes(m_accCircle, m_velCircle, m_posCircle);

                // non-uniform length discretization
                std::srand(6);
                Eigen::ArrayXd m_lengths { Eigen::ArrayXd::LinSpaced(state.range(1), 0.0, 3 / 2 * M_PI) * m_circle.radius() };

                m_circlePoly = Polychain { m_circle(m_lengths) };
                m_circleTransform = Transform { std::make_shared<Polychain<Eigen::Dynamic>>(m_circlePoly) };
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
        };

        BENCHMARK_DEFINE_F(PolylineBenchmark, PosCartes)(benchmark::State& state)
        {
            Points<Eigen::Dynamic> posCartes {};
            for(auto _ : state)
                posCartes = m_circleTransform.posCartes(m_posFrenet);
            reportError(posCartes - m_posCartes, state);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, PosCartes)->Ranges({{8, 8<<10}, {8, 8<<10}})->ArgNames({"NumQueries", "NumPoints"});

        BENCHMARK_DEFINE_F(PolylineBenchmark, VelCartes)(benchmark::State& state)
        {
            Points<Eigen::Dynamic> velCartes {};
            for(auto _ : state)
                velCartes = m_circleTransform.velCartes(m_velFrenet, m_posFrenet);
            reportError(velCartes - m_velCartes, state);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, VelCartes)->Ranges({{8, 8<<10}, {8, 8<<10}})->ArgNames({"NumQueries", "NumPoints"});

        BENCHMARK_DEFINE_F(PolylineBenchmark, AccCartes)(benchmark::State& state)
        {
            Points<Eigen::Dynamic> accCartes {};
            for(auto _ : state)
                accCartes = m_circleTransform.accCartes(m_accFrenet, m_velFrenet, m_posFrenet);
            reportError(accCartes - m_accCartes, state);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, AccCartes)->Ranges({{8, 8<<10}, {8, 8<<10}})->ArgNames({"NumQueries", "NumPoints"});

        BENCHMARK_DEFINE_F(PolylineBenchmark, PosFrenet)(benchmark::State& state)
        {
            Points<Eigen::Dynamic> posFrenet {};
            for(auto _ : state)
                posFrenet = m_circleTransform.posFrenet(m_posCartes);
            reportError(posFrenet - m_posFrenet, state);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, PosFrenet)->Ranges({{8, 8<<10}, {8, 8<<10}})->ArgNames({"NumQueries", "NumPoints"});


        BENCHMARK_DEFINE_F(PolylineBenchmark, VelFrenet)(benchmark::State& state)
        {
            Points<Eigen::Dynamic> velFrenet {};
            for(auto _ : state)
                velFrenet = m_circleTransform.velFrenet(m_velCartes, m_posFrenet);
            reportError(velFrenet - m_velFrenet, state);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, VelFrenet)->Ranges({{8, 8<<10}, {8, 8<<10}})->ArgNames({"NumQueries", "NumPoints"});

        BENCHMARK_DEFINE_F(PolylineBenchmark, AccFrenet)(benchmark::State& state)
        {
            Points<Eigen::Dynamic> accFrenet {};
            for(auto _ : state)
                accFrenet = m_circleTransform.accFrenet(m_accCartes, m_velFrenet, m_posFrenet);
            reportError(accFrenet - m_accFrenet, state);
        }

        BENCHMARK_REGISTER_F(PolylineBenchmark, AccFrenet)->Ranges({{8, 8<<10}, {8, 8<<10}})->ArgNames({"NumQueries", "NumPoints"});
    };
};

BENCHMARK_MAIN();