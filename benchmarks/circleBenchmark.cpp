#include <benchmark/benchmark.h>

#include "frenetTransform/internal/circle.h"
#include "frenetTransform/internal/constexprTypes.h"
#include "frenetTransform/internal/transformCircle.h"
#include "frenetTransform/point.h"
#include "frenetTransform/points.h"
#include "frenetTransform/polychain.h"

namespace FrenetTransform {
namespace Internal {
/**
 * @brief Class to benchmark transformations of states along a polychain.
 *
 * @tparam Params struct with benchmark parameters.
 */
template <typename Params> class PolylineBenchmark : public benchmark::Fixture {
protected:
  static constexpr int numPoints{Params::s_vals[0]}; /**<< number of polychain points */
  static constexpr int numQueries{Params::s_vals[1]}; /**<< number of query points */

  using ArrayQueries = Eigen::Array<double, Params::s_vals[1], 1>;
  using ArrayPoints = Eigen::Array<double, Params::s_vals[0], 1>;

  const Circle<numQueries> m_circle{5.0, {0.0, 0.0}, -M_PI}; /**<< Circle of radius 5 starting at -PI */
  const TransformCircle<numQueries> m_transform{
      std::make_shared<Circle<numQueries>>(m_circle)}; /**<< Circle Transform as transformation ground truth */

  const Circle<numPoints> m_circleApprox{5.0, {0.0, 0.0}, -M_PI}; /**<< Circle to approximate by Polychain with "numPoints" points */

  Points<numQueries> m_posFrenet{}; /**<< holds ground truth points in Frenet frame */
  Points<numQueries> m_velFrenet{}; /**<< holds ground truth velocities in Frenet frame */
  Points<numQueries> m_accFrenet{}; /**<< holds ground truth acceleratios in Frenet frame */

  Points<numQueries> m_posCartes{}; /**<< holds ground truth points in Cartesian frame */
  Points<numQueries> m_velCartes{}; /**<< holds ground truth velocities in Cartesian frame */
  Points<numQueries> m_accCartes{}; /**<< holds ground truth acceleratios in Cartesian frame */

  Polychain<numPoints, numQueries> m_circlePoly{}; /**<< Polychain approximating "m_circleApprox" */
  Transform<numQueries> m_circleTransform{}; /**<< Transform using "m_circlePoly" */

  /**
   * @brief Setup randomized query points, the Polychain and the Transform.
   *
   * @param state benchmark state containing benchmark parameters.
   */
  void SetUp(::benchmark::State &state) {
    std::srand(0);
    const ArrayQueries posCircleX{
        m_circle.radius() * (1 + ArrayQueries::Random(state.range(0)) * 0.95)};
    std::srand(1);
    const ArrayQueries posCircleY{-ArrayQueries::Random(state.range(0)).abs() *
                                  M_PI};

    const Points<numQueries, PointCircle> m_posCircle{posCircleX, posCircleY};
    m_posFrenet = m_transform.posFrenet(m_posCircle);
    m_posCartes = m_transform.posCartes(m_posCircle);

    std::srand(2);
    const ArrayQueries velCircleX{
        m_circle.radius() * (1 + ArrayQueries::Random(state.range(0)) * 0.95)};
    std::srand(3);
    const ArrayQueries velCircleY{ArrayQueries::Random(state.range(0)) * M_PI /
                                  4};

    const Points<numQueries, PointCircle> m_velCircle{velCircleX, velCircleY};
    m_velFrenet = m_transform.velFrenet(m_velCircle);
    m_velCartes = m_transform.velCartes(m_velCircle, m_posCircle);

    std::srand(4);
    const ArrayQueries accCircleX{
        m_circle.radius() * (1 + ArrayQueries::Random(state.range(0)) * 0.95)};
    std::srand(5);
    const ArrayQueries accCircleY{ArrayQueries::Random(state.range(0)) * M_PI /
                                  4};

    const Points<numQueries, PointCircle> m_accCircle{accCircleX, accCircleY};
    m_accFrenet = m_transform.accFrenet(m_accCircle);
    m_accCartes = m_transform.accCartes(m_accCircle, m_velCircle, m_posCircle);

    // non-uniform length discretization
    std::srand(6);
    ArrayPoints m_lengths{
        ArrayPoints::LinSpaced(state.range(1), 0.0, 3 / 2 * M_PI) *
        m_circle.radius()};

    m_circlePoly = Polychain<numPoints, numQueries>{m_circleApprox(m_lengths)};
    m_circleTransform = Transform<numQueries>{
        std::make_shared<Polychain<numPoints, numQueries>>(m_circlePoly)};
  }

  template <typename Tpoint>
  void reportError(const Points<numQueries, Tpoint> diff,
                   benchmark::State &state) {
    Eigen::ArrayXd errX{diff.x().abs()};
    std::sort(errX.begin(), errX.end());
    state.counters["ErrMaxX"] = errX(errX.rows() - 1);
    state.counters["ErrMedX"] = errX(errX.rows() / 2);
    Eigen::ArrayXd errY{diff.y().abs()};
    std::sort(errY.begin(), errY.end());
    state.counters["ErrMaxY"] = errY(errY.rows() - 1);
    state.counters["ErrMedY"] = errY(errY.rows() / 2);
  }

  void posCartes(benchmark::State &state) {
    Points<numQueries> posCartes{};
    for (auto _ : state)
      posCartes = m_circleTransform.posCartes(m_posFrenet);
    reportError(posCartes - m_posCartes, state);
  }

  void velCartes(benchmark::State &state) {
    Points<numQueries> velCartes{};
    for (auto _ : state)
      velCartes = m_circleTransform.velCartes(m_velFrenet, m_posFrenet);
    reportError(velCartes - m_velCartes, state);
  }

  void accCartes(benchmark::State &state) {
    Points<numQueries> accCartes{};
    for (auto _ : state)
      accCartes =
          m_circleTransform.accCartes(m_accFrenet, m_velFrenet, m_posFrenet);
    reportError(accCartes - m_accCartes, state);
  }

  void posFrenet(benchmark::State &state) {
    Points<numQueries> posFrenet{};
    for (auto _ : state)
      posFrenet = m_circleTransform.posFrenet(m_posCartes);
    reportError(posFrenet - m_posFrenet, state);
  }

  void velFrenet(benchmark::State &state) {
    Points<numQueries> velFrenet{};
    for (auto _ : state)
      velFrenet = m_circleTransform.velFrenet(m_velCartes, m_posFrenet);
    reportError(velFrenet - m_velFrenet, state);
  }

  void accFrenet(benchmark::State &state) {
    Points<numQueries> accFrenet{};
    for (auto _ : state)
      accFrenet =
          m_circleTransform.accFrenet(m_accCartes, m_velFrenet, m_posFrenet);
    reportError(accFrenet - m_accFrenet, state);
  }
};

using Dynamic =
    ConstVals<int, 2, std::array<int, 2>{Eigen::Dynamic, Eigen::Dynamic}>;
using Static8 = ConstVals<int, 2, std::array<int, 2>{8, 8}>;
using Static512 = ConstVals<int, 2, std::array<int, 2>{512, 512}>;
using Static4096 = ConstVals<int, 2, std::array<int, 2>{4096, 4096}>;

// Benchmark: posCartes
BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, PosCartesDynamic, Dynamic)
(benchmark::State &state) {
  posCartes(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, PosCartesDynamic)
    ->Ranges({{8, 8 << 10}, {8, 8 << 10}})
    ->ArgNames({"NumQueries", "NumPoints"});

BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, PosCartesStatic8, Static8)
(benchmark::State &state) {
  posCartes(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, PosCartesStatic8)
    ->Args({8, 8})
    ->ArgNames({"NumQueries", "NumPoints"});

BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, PosCartesStatic512, Static512)
(benchmark::State &state) {
  posCartes(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, PosCartesStatic512)
    ->Args({512, 512})
    ->ArgNames({"NumQueries", "NumPoints"});

BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, PosCartesStatic4096, Static4096)
(benchmark::State &state) {
  posCartes(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, PosCartesStatic4096)
    ->Args({4096, 4096})
    ->ArgNames({"NumQueries", "NumPoints"});

// Benchmark: velCartes
BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, VelCartesDynamic, Dynamic)
(benchmark::State &state) {
  velCartes(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, VelCartesDynamic)
    ->Ranges({{8, 8 << 10}, {8, 8 << 10}})
    ->ArgNames({"NumQueries", "NumPoints"});

// Benchmark: accCartes
BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, AccCartesDynamic, Dynamic)
(benchmark::State &state) {
  accCartes(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, AccCartesDynamic)
    ->Ranges({{8, 8 << 10}, {8, 8 << 10}})
    ->ArgNames({"NumQueries", "NumPoints"});

// Benchmark: posFrenet
BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, PosFrenetDyn, Dynamic)
(benchmark::State &state) {
  posFrenet(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, PosFrenetDyn)
    ->Ranges({{8, 8 << 10}, {8, 8 << 10}})
    ->ArgNames({"NumQueries", "NumPoints"});

BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, PosFrenetStatic8, Static8)
(benchmark::State &state) {
  posFrenet(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, PosFrenetStatic8)
    ->Args({8, 8})
    ->ArgNames({"NumQueries", "NumPoints"});

BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, PosFrenetStatic512, Static512)
(benchmark::State &state) {
  posFrenet(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, PosFrenetStatic512)
    ->Args({512, 512})
    ->ArgNames({"NumQueries", "NumPoints"});

BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, PosFrenetStatic4096, Static4096)
(benchmark::State &state) {
  posFrenet(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, PosFrenetStatic4096)
    ->Args({4096, 4096})
    ->ArgNames({"NumQueries", "NumPoints"});

// Benchmark: velFrenet
BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, VelFrenetDyn, Dynamic)
(benchmark::State &state) {
  velFrenet(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, VelFrenetDyn)
    ->Ranges({{8, 8 << 10}, {8, 8 << 10}})
    ->ArgNames({"NumQueries", "NumPoints"});

// Benchmark: accFrenet
BENCHMARK_TEMPLATE_DEFINE_F(PolylineBenchmark, AccFrenetDyn, Dynamic)
(benchmark::State &state) {
  accFrenet(state);
}
BENCHMARK_REGISTER_F(PolylineBenchmark, AccFrenetDyn)
    ->Ranges({{8, 8 << 10}, {8, 8 << 10}})
    ->ArgNames({"NumQueries", "NumPoints"});
}; // namespace Internal
}; // namespace FrenetTransform

BENCHMARK_MAIN();