#ifndef LINE_H
#define LINE_H

#include <Eigen/Core>
#include <math.h>

#include "frenetTransform/path.h"
#include "frenetTransform/point.h"
#include "frenetTransform/internal/math.h"

namespace FrenetTransform
{
    namespace Internal
    {
        template<int NumQueries=Eigen::Dynamic>
        class Line : public Path<NumQueries>
        {
        public:
            using ArrayQueries = Eigen::Array<double, NumQueries, 1>;

            Line() = delete;

            Line(const Point& start, const Point& end)
                : m_start { start }
                , m_end { end }
            {
            }

            Points<NumQueries> operator()(const ArrayQueries& lengths) const override
            {
                const double lineLength { m_end.distance(m_start) };
                ArrayQueries relLengths { lengths / lineLength};

                for(auto& relLength : relLengths)
                    relLength = std::clamp(relLength, 0.0, 1.0);

                return { m_start.x() * relLengths + (1 - relLengths) * m_end.x(),
                         m_start.y() * relLengths + (1 - relLengths) * m_end.y() };
            }

            ArrayQueries lengths(const Points<NumQueries>& points) const override
            {
                const auto xDiffPnt { m_end.x() - points.x() };
                const auto yDiffPnt { m_end.y() - points.y() };
                const Point pointDiff { m_end - m_start };

                ArrayQueries relLengths { (xDiffPnt * pointDiff.x() + yDiffPnt * pointDiff.y())
                    / ( std::pow(pointDiff.x(), 2) + std::pow(pointDiff.y(), 2)) };

                for(auto& relLength : relLengths)
                    relLength = std::clamp(relLength, 0.0, 1.0);

                return relLengths * m_end.distance(m_start);
            }

        private:
            const Point m_start;
            const Point m_end;

            Points<NumQueries> gradient1(const ArrayQueries& lengths) const override
            {
                ArrayQueries gradX(lengths.rows());
                gradX += m_end.x() - m_start.x() / m_end.distance(m_start);
                ArrayQueries gradY(lengths.rows());
                gradY += m_end.y() - m_start.y() / m_end.distance(m_start);
                return { gradX, gradY };
            }

            Points<NumQueries> gradient2(const ArrayQueries& lengths) const override { return { ArrayQueries::Zero(lengths.rows()), ArrayQueries::Zero(lengths.rows()) }; }

            Points<NumQueries> gradient3(const ArrayQueries& lengths) const override { return { ArrayQueries::Zero(lengths.rows()), ArrayQueries::Zero(lengths.rows()) }; }
        };
    };
};

#endif