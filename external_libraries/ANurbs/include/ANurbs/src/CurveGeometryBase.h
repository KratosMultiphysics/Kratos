#pragma once

#include "CurveShapeEvaluator.h"
#include "Interval.h"
#include "Knots.h"
#include "VectorTraits.h"

#include <vector>

namespace ANurbs {

template <typename TVector>
class CurveGeometryBase
{
public:
    using VectorType = TVector;
    using ScalarType = ScalarTypeOf<VectorType>;
    using KnotsType = std::vector<ScalarType>;
    using IntervalType = Interval<ScalarType>;

protected:
    int m_degree;
    KnotsType m_knots;

public:
    CurveGeometryBase(
        const int degree,
        const int nbPoles)
        : m_degree(degree)
        , m_knots(nbPoles + degree - 1)
    {
    }

    static constexpr int
    Dimension()
    {
        return DimensionOf<VectorType>();
    }

    int
    Degree() const
    {
        return m_degree;
    }

    int
    NbKnots() const
    {
        return static_cast<int>(m_knots.size());
    }

    ScalarType
    Knot(
        const int index) const
    {
        return m_knots[index];
    }

    void
    SetKnot(
        const int index,
        const ScalarType value)
    {
        m_knots[index] = value;
    }

    const KnotsType&
    Knots() const
    {
        return m_knots;
    }

    int
    NbPoles() const
    {
        return NbKnots() - Degree() + 1;
        ;
    }

    virtual VectorType
    Pole(
        const int index) const = 0;

    virtual void
    SetPole(
        const int index,
        const VectorType& value)
        = 0;

    std::vector<VectorType>
    Poles() const
    {
        std::vector<VectorType> poles(NbPoles());

        for (int i = 0; i < poles.size(); i++) {
            poles[i] = Pole(i);
        }

        return poles;
    }

    virtual bool
    IsRational() const = 0;

    virtual ScalarType
    Weight(
        const int index) const = 0;

    virtual void
    SetWeight(
        const int index,
        const ScalarType value)
        = 0;

    std::vector<ScalarType>
    Weights() const
    {
        std::vector<ScalarType> weights(NbPoles());

        for (int i = 0; i < weights.size(); i++) {
            weights[i] = Weight(i);
        }

        return weights;
    }

    VectorType
    WeightedPole(
        const int index) const
    {
        return Pole(index) * Weight(index);
    }

    IntervalType
    Domain()
    {
        ScalarType t0 = Knot(Degree() - 1);
        ScalarType t1 = Knot(NbKnots() - Degree());

        return IntervalType(t0, t1);
    }

    int
    SpanAt(
        const ScalarType t) const
    {
        return Knots::UpperSpan(Degree(), Knots(), t);
    }

    template <typename TValue, typename TValues>
    TValue
    EvaluateAt(
        TValues values,
        const ScalarType t) const
    {
        // compute shape functions

        CurveShapeEvaluator<ScalarType> shape(Degree(), 0);

        if (IsRational()) {
            shape.Compute(Knots(), [&](int i) -> ScalarType {
                return Weight(i);
            }, t);
        } else {
            shape.Compute(Knots(), t);
        }

        // compute point

        TValue value;

        for (int i = 0; i < shape.NbNonzeroPoles(); i++) {
            int index = shape.FirstNonzeroPole() + i;

            if (i == 0) {
                value = values(index) * shape(0, i);
            } else {
                value += values(index) * shape(0, i);
            }
        }

        return value;
    }

    template <typename TValue, typename TValues>
    std::vector<TValue>
    EvaluateAt(
        TValues values,
        const ScalarType t,
        const int order) const
    {
        // evaluate shape functions

        CurveShapeEvaluator<ScalarType> shape(Degree(), order);

        if (IsRational()) {
            shape.Compute(Knots(), [&](int i) -> ScalarType {
                return Weight(i);
            }, t);
        } else {
            shape.Compute(Knots(), t);
        }

        // calculate derivatives

        std::vector<TValue> derivatives(shape.NbShapes());

        for (int order = 0; order < shape.NbShapes(); order++) {
            for (int i = 0; i < shape.NbNonzeroPoles(); i++) {
                int index = shape.FirstNonzeroPole() + i;

                if (i == 0) {
                    derivatives[order] = values(index) * shape(order, i);
                } else {
                    derivatives[order] += values(index) * shape(order, i);
                }
            }
        }

        return derivatives;
    }

    std::vector<IntervalType>
    Spans() const
    {
        int firstSpan = Degree() - 1;
        int lastSpan = NbKnots() - Degree() - 1;

        int nbSpans = lastSpan - firstSpan + 1;

        std::vector<IntervalType> result(nbSpans);

        for (int i = 0; i < nbSpans; i++) {
            ScalarType t0 = Knot(firstSpan + i);
            ScalarType t1 = Knot(firstSpan + i + 1);

            result[i] = IntervalType(t0, t1);
        }

        return result;
    }

    VectorType
    PointAt(
        const ScalarType t) const
    {
        auto poles = [&](int i) -> VectorType { return Pole(i); };

        return EvaluateAt<VectorType>(poles, t);
    }

    std::vector<VectorType>
    DerivativesAt(
        const ScalarType t,
        const int order) const
    {
        auto poles = [&](int i) -> VectorType { return Pole(i); };

        return EvaluateAt<VectorType>(poles, t, order);
    }
};

using CurveGeometryBase1D = CurveGeometryBase<Point1D>;
using CurveGeometryBase2D = CurveGeometryBase<Point2D>;
using CurveGeometryBase3D = CurveGeometryBase<Point3D>;

} // namespace ANurbs
