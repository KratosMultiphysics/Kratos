#pragma once

#include "CurveBase.h"
#include "CurveGeometry.h"
#include "CurveTessellation.h"
#include "SurfaceGeometry.h"
#include "Pointer.h"

#include <array>
#include <utility>

namespace ANurbs {

template <typename TVector>
class CurveSpanIntersection
{
public:
    using VectorType = TVector;
    using ScalarType = ScalarTypeOf<VectorType>;

    using CurveType = CurveBase<VectorType>;

private:
    using ParameterPoint = std::pair<ScalarType, VectorType>;

private:
    std::vector<ScalarType> m_intersectionParameters;

public:
    CurveSpanIntersection()
    {
        static_assert(CurveType::Dimension() == 2,
            "Only planar curves are supported");
    }

    int
    NbIntersections() const
    {
        return static_cast<int>(m_intersectionParameters.size());
    }

    ScalarType
    Parameter(
        const int& index) const
    {
        return m_intersectionParameters[index];
    }

private:
    struct Axis
    {
        int m_index;
        std::vector<ScalarType> m_values;
        ScalarType m_tolerance;

        inline ScalarType
        GetValue(
            const VectorType& point) const
        {
            return point[m_index];
        }

        void Initialize(
            const int axis,
            const std::vector<ScalarType>& knots,
            ScalarType tolerance)
        {
            m_index = axis;

            m_values.clear();

            ScalarType mInf = std::numeric_limits<ScalarType>::lowest();
            ScalarType pInf = std::numeric_limits<ScalarType>::max();

            m_values.push_back(mInf);

            for (ScalarType knot : knots) {
                if (std::abs(m_values.back() - knot) > tolerance) {
                    m_values.push_back(knot);
                }
            }

            m_values.push_back(pInf);

            m_tolerance = tolerance;
        }

        void
        Intersect(
            const CurveType& curve,
            const ParameterPoint& a,
            const ParameterPoint& b,
            std::vector<ScalarType>& parameters)
        {
            ScalarType tA = std::get<0>(a);
            ScalarType valueA = GetValue(std::get<1>(a));

            ScalarType tB = std::get<0>(b);
            ScalarType valueB = GetValue(std::get<1>(b));
            
            // make sure that valueA <= valueB
            if (valueA > valueB) {
                std::swap(valueA, valueB);
                std::swap(tA, tB);
            }
            
            // index of the first intersect value
            std::size_t indexA;
            {
                auto it = std::lower_bound(std::begin(m_values),
                    std::end(m_values), valueA - m_tolerance);
                indexA = std::distance(std::begin(m_values), it);
            }

            // index of the first non intersect value
            std::size_t indexB;
            {
                auto it = std::upper_bound(std::begin(m_values),
                    std::end(m_values), valueB + m_tolerance);
                indexB = std::distance(std::begin(m_values), it);
            }

            // find intersections
            for (std::size_t i = indexA; i < indexB; i++) {
                ScalarType target = m_values[i];
            
                ScalarType t = tA;

                ScalarType delta = valueA - valueB;

                if (delta != 0) {
                    t += (valueA - target) / delta * (tB - tA);
                }

                for (int j = 0; j < 100; j++) {
                    auto c = curve.DerivativesAt(t, 1);

                    ScalarType f = GetValue(c[0]) - target;
                    
                    if (std::abs(f) < m_tolerance) {
                        break;
                    }

                    ScalarType df = GetValue(c[1]);

                    if (df != 0) {
                        t -= f / df;
                    } else {
                        break;
                    }
                }

                // FIXME: check for convergency

                if (curve.Domain().Contains(t)) {
                parameters.push_back(t);
            }
        }
        }
    };

    template <typename TContainer>
    static void
    UniqueSorted(
        TContainer& container,
        const ScalarType& tolerance = 1e-7)
    {
        std::sort(std::begin(container), std::end(container));

        auto last = std::unique(std::begin(container), std::end(container),
            [=](ScalarType a, ScalarType b) -> bool {
                return b - a < tolerance;
            });

        auto nbUnique = std::distance(std::begin(container), last);

        container.resize(nbUnique);
    }

public:
    void
    Compute(
        const CurveType& curve,
        const std::vector<ScalarType>& knotsU,
        const std::vector<ScalarType>& knotsV,
        const ScalarType tolerance,
        const bool includeCurveKnots)
    {
        // approximate curve with a polyline

        CurveTessellation<VectorType> tessellation;

        tessellation.Compute(curve, tolerance);

        // initialize axes

        std::array<Axis, 2> axes;

        axes[0].Initialize(0, knotsU, tolerance);
        axes[1].Initialize(1, knotsV, tolerance);

        // add curve knots

        if (includeCurveKnots) {
            ScalarType t0 = tessellation.Parameter(0);

            m_intersectionParameters.push_back(t0);

            for (const auto& span : curve.Spans()) {
                m_intersectionParameters.push_back(span.T1());
            }
        }

        // check line segments

        for (int i = 1; i < tessellation.NbPoints(); i++) {
            ParameterPoint a = {
                tessellation.Parameter(i - 1),
                tessellation.Point(i - 1)
            };

            ParameterPoint b = {
                tessellation.Parameter(i - 0),
                tessellation.Point(i - 0)
            };

            for (auto& axis : axes) {
                axis.Intersect(curve, a, b, m_intersectionParameters);
            }
        }

        UniqueSorted(m_intersectionParameters);
    }
};

using CurveSpanIntersection2D = CurveSpanIntersection<Point2D>;

} // namespace ANurbs
