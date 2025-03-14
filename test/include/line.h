#ifndef LINE_H
#define LINE_H

#include <Eigen/Core>
#include <math.h>

#include "frenetTransform/path.h"
#include "frenetTransform/point.h"
#include "frenetTransform/math.h"

namespace FrenetTransform
{
    namespace Testing
    {
        class Line : public Path
        {
        public:
            Line() = delete;

            Line(const Point& start, const Point& end)
                : m_start { start }
                , m_end { end }
            {
            }

            Points<Eigen::Dynamic> operator()(const Eigen::ArrayXd& lengths) const override
            {
                const double lineLength { m_end.distance(m_start) };
                Eigen::ArrayXd relLengths { lengths / lineLength};

                for(auto& relLength : relLengths)
                    relLength = std::clamp(relLength, 0.0, 1.0);

                return { m_start.x() * relLengths + (1 - relLengths) * m_end.x(),
                         m_start.y() * relLengths + (1 - relLengths) * m_end.y() };
            }

            Eigen::ArrayXd lengths(const Points<Eigen::Dynamic>& points) const override
            {
                const auto xDiffPnt { m_end.x() - points.x() };
                const auto yDiffPnt { m_end.y() - points.y() };
                const Point pointDiff { m_end - m_start };

                Eigen::ArrayXd relLengths { (xDiffPnt * pointDiff.x() + yDiffPnt * pointDiff.y())
                    / ( std::pow(pointDiff.x(), 2) + std::pow(pointDiff.y(), 2)) };

                for(auto& relLength : relLengths)
                    relLength = std::clamp(relLength, 0.0, 1.0);

                return relLengths * m_end.distance(m_start);
            }

        private:
            const Point m_start;
            const Point m_end;

            Points<Eigen::Dynamic> gradient1(const Eigen::ArrayXd& lengths) const override
            {
                Eigen::ArrayXd gradX(lengths.rows());
                gradX += m_end.x() - m_start.x() / m_end.distance(m_start);
                Eigen::ArrayXd gradY(lengths.rows());
                gradY += m_end.y() - m_start.y() / m_end.distance(m_start);
                return { gradX, gradY };
            }

            Points<Eigen::Dynamic> gradient2(const Eigen::ArrayXd& lengths) const override { return { Eigen::ArrayXd::Zero(lengths.rows()), Eigen::ArrayXd::Zero(lengths.rows()) }; }

            Points<Eigen::Dynamic> gradient3(const Eigen::ArrayXd& lengths) const override { return { Eigen::ArrayXd::Zero(lengths.rows()), Eigen::ArrayXd::Zero(lengths.rows()) }; }
        };
    };
};

#endif