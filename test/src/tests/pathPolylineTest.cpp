#include "frenetTransform/polychain.h"
#include "frenetTransform/points.h"
#include "frenetTransform/transform.h"
#include "circle.h"
#include "pathCircleTest.h"

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <math.h>
#include <memory>

namespace FrenetTransform
{
    namespace Testing
    {
        // template <typename NumPoints, typename NumQueries>
        template <typename Params>
        class PathPolylineTest : public PathCircleTest<Integral<Params::s_vals[1]>>
        {
        public:
            using ArrayQueries = Eigen::Array<double, Params::s_vals[1], 1>;
            using ArrayPoints = Eigen::Array<double, Params::s_vals[0], 1>;

            static constexpr int s_numPoints { Params::s_vals[0] == -1 ? 4096 : Params::s_vals[0] };

        protected:
            const Circle<Params::s_vals[0]> m_circleApprox { 5.0, {0.0, 0.0}, -M_PI };
            const Polychain<Params::s_vals[0], Params::s_vals[1]> m_circlePoly { m_circleApprox(ArrayPoints::LinSpaced(s_numPoints, 0.0, 2 * M_PI) * m_circleApprox.radius()) };
            const Transform<Params::s_vals[1]> m_circleTransform { std::make_shared<Polychain<Params::s_vals[0], Params::s_vals[1]>>(m_circlePoly) };

            const Points<Params::s_vals[1]> m_posCartes { this->m_transform.posCartes(this->m_posCircle) };

            const Points<Params::s_vals[1], PointCircle> m_velCircle { ArrayQueries::Random(this->s_numQueries).abs() * this->m_circle.radius(), ArrayQueries::Random(this->s_numQueries) * M_PI / 4 };
            const Points<Params::s_vals[1]> m_velCartes { this->m_transform.velCartes(this->m_velCircle, this->m_posCircle) };
            const Points<Params::s_vals[1]> m_velFrenet { this->m_transform.velFrenet(this->m_velCircle) };

            const Points<Params::s_vals[1], PointCircle> m_accCircle { ArrayQueries::Random(this->s_numQueries).abs() * this->m_circle.radius(), ArrayQueries::Random(this->s_numQueries) * M_PI / 4 };
            const Points<Params::s_vals[1]> m_accCartes { this->m_transform.accCartes(this->m_accCircle, this->m_velCircle, this->m_posCircle) };
            const Points<Params::s_vals[1]> m_accFrenet { this->m_transform.accFrenet(this->m_accCircle) };
        };

        using TestParams = testing::Types<
            ConstVals<int, 2, std::array<int, 2> {{4096, 100}}>,
            ConstVals<int, 2, std::array<int, 2> {{4096, Eigen::Dynamic}}>,
            ConstVals<int, 2, std::array<int, 2> {{Eigen::Dynamic, 100}}>,
            ConstVals<int, 2, std::array<int, 2> {{Eigen::Dynamic, Eigen::Dynamic}}>
        >;
        TYPED_TEST_SUITE(PathPolylineTest, TestParams);

        TYPED_TEST(PathPolylineTest, GetPointsCircle)
        {
            const auto pointsCircleEst { this->m_circlePoly(this->m_posFrenet.x()) };
            const auto pointsCircleGtr { this->m_circle(this->m_posFrenet.x()) };

            this->expectAllClose(pointsCircleEst.x(), pointsCircleGtr.x(), 1e-2);
            this->expectAllClose(pointsCircleEst.y(), pointsCircleGtr.y(), 1e-2);
        }

        TYPED_TEST(PathPolylineTest, GetTangentsCircle)
        {
            const auto tangentsEst { this->m_circlePoly.tangent(this->m_posFrenet.x()) };
            const auto tangentsGtr { this->m_circle.tangent(this->m_posFrenet.x()) };

            this->expectAllClose(tangentsEst.x(), tangentsGtr.x(), 1e-2);
            this->expectAllClose(tangentsEst.y(), tangentsGtr.y(), 1e-2);
        }

        TYPED_TEST(PathPolylineTest, GetNormalsCircle)
        {
            const auto normalsEst { this->m_circlePoly.normal(this->m_posFrenet.x()) };
            const auto normalsGtr { this->m_circle.normal(this->m_posFrenet.x()) };

            this->expectAllClose(normalsEst.x(), normalsGtr.x(), 1e-2);
            this->expectAllClose(normalsEst.y(), normalsGtr.y(), 1e-2);
        }

        TYPED_TEST(PathPolylineTest, GetAnglesCircle)
        {
            const auto angles0Est { this->m_circlePoly.angle0(this->m_posFrenet.x()) };
            const auto angles0Gtr { this->m_circle.angle0(this->m_posFrenet.x()) };

            this->expectAllClose(angles0Est, angles0Gtr, 1e-2);
        }

        TYPED_TEST(PathPolylineTest, GetAngles1Circle)
        {
            const auto angles1Est { this->m_circlePoly.angle1(this->m_posFrenet.x()) };
            const auto angles1Gtr { this->m_circle.angle1(this->m_posFrenet.x()) };

            this->expectAllClose(angles1Est, angles1Gtr, 1e-2);
        }

        TYPED_TEST(PathPolylineTest, GetAngles2Circle)
        {
            const auto angles2Est { this->m_circlePoly.angle2(this->m_posFrenet.x()) };
            const auto angles2Gtr { this->m_circle.angle2(this->m_posFrenet.x()) };

            this->expectAllClose(angles2Est,angles2Gtr, 1e-2);
        }

        TYPED_TEST(PathPolylineTest, NextPointsCircle)
        {
            const auto lengthsEst { this->m_circlePoly.lengths(this->m_posCartes) };
            const auto lengthsGtr { this->m_circle.lengths(this->m_posCartes) };

            this->expectAllClose(lengthsEst, lengthsGtr, 2e-2);
        }

        TYPED_TEST(PathPolylineTest, PosFrenetCircle)
        {
            // determine frenet positions
            const auto posFrenet { this->m_circleTransform.posFrenet(this->m_posCartes) };

            this->expectAllClose(posFrenet.x(), this->m_posFrenet.x(), 4e-3);
            this->expectAllClose(posFrenet.y(), this->m_posFrenet.y(), 4e-3);
        }

        TYPED_TEST(PathPolylineTest, PosCartCircle)
        {
            // determine frenet positions
            const auto posCartes { this->m_circleTransform.posCartes(this->m_posFrenet) };

            this->expectAllClose(posCartes.x(), this->m_posCartes.x(), 1.6e-2);
            this->expectAllClose(posCartes.y(), this->m_posCartes.y(), 1.6e-2);
        }

        TYPED_TEST(PathPolylineTest, VelFrenetCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto velFrenet { this->m_circleTransform.velFrenet(this->m_velCartes, this->m_posFrenet) };

            this->expectAllClose(velFrenet.x(), this->m_velFrenet.x(), 5e-1);
            this->expectAllClose(velFrenet.y(), this->m_velFrenet.y(), 5e-1);
        }

        TYPED_TEST(PathPolylineTest, VelCartCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto velCartes { this->m_circleTransform.velCartes(this->m_velFrenet, this->m_posFrenet) };

            this->expectAllClose(velCartes.x(), this->m_velCartes.x(), 6e-2);
            this->expectAllClose(velCartes.y(), this->m_velCartes.y(), 2e-2);
        }

        TYPED_TEST(PathPolylineTest, AccFrenetCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto accFrenet { this->m_circleTransform.accFrenet(this->m_accCartes, this->m_velFrenet, this->m_posFrenet) };

            this->expectAllClose(accFrenet.x(), this->m_accFrenet.x(), 7e-2);
            this->expectAllClose(accFrenet.y(), this->m_accFrenet.y(), 7e-2);
        }

        TYPED_TEST(PathPolylineTest, AccCartCircle)
        {
            // test frenet frame results since error grows with distance from path
            const auto accCartes { this->m_circleTransform.accCartes(this->m_accFrenet, this->m_velFrenet, this->m_posFrenet) };

            this->expectAllClose(accCartes.x(), this->m_accCartes.x(), 5e-2);
            this->expectAllClose(accCartes.y(), this->m_accCartes.y(), 5e-2);
        }
    };
};

int main(int argc, char **argv)
{
    std::srand(0);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
