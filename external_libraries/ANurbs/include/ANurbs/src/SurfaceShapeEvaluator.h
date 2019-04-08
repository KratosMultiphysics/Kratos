#pragma once

#include "CurveShapeEvaluator.h"
#include "Math.h"
#include "VectorTraits.h"

#include <vector>

namespace ANurbs {

    template <typename TScalar>
    class SurfaceShapeEvaluator
    {
    public:
        using ScalarType = TScalar;

    private:
        int m_order;
        CurveShapeEvaluator<ScalarType> m_shapeU;
        CurveShapeEvaluator<ScalarType> m_shapeV;
        std::vector<ScalarType> m_weightedSums;
        std::vector<ScalarType> m_values;
        int m_firstNonzeroPoleU;
        int m_firstNonzeroPoleV;

    public:
        static constexpr inline int
            NbShapes(
                const int order) noexcept
        {
            return (1 + order) * (2 + order) / 2;
        }

        static constexpr inline int
            ShapeIndex(
                const int derivativeU,
                const int derivativeV) noexcept
        {
            return derivativeV + (derivativeU + derivativeV) * (1 + derivativeU +
                derivativeV) / 2;
        }

    private:
        ScalarType &
            WeightedSum(
                const int index)
        {
            return m_weightedSums[index];
        }

        ScalarType&
            WeightedSum(
                const int derivativeU,
                const int derivativeV)
        {
            const int index = ShapeIndex(derivativeU, derivativeV);

            return WeightedSum(index);
        }

        inline int
            Index(
                const int derivative,
                const int poleU,
                const int poleV) const
        {
            const int pole = Math::MatrixIndex(NbNonzeroPolesU(), NbNonzeroPolesV(),
                poleU, poleV);
            const int index = Math::MatrixIndex(NbShapes(), NbNonzeroPoles(),
                derivative, pole);

            return index;
        }

        ScalarType&
            Value(
                const int derivative,
                const int pole)
        {
            const int index = Math::MatrixIndex(NbShapes(), NbNonzeroPoles(),
                derivative, pole);

            return m_values[index];
        }

        ScalarType&
            Value(
                const int derivative,
                const int poleU,
                const int poleV)
        {
            const int index = Index(derivative, poleU, poleV);

            return m_values[index];
        }

    public:
        SurfaceShapeEvaluator()
        {
        }

        SurfaceShapeEvaluator(
            const int degreeU,
            const int degreeV,
            const int order)
        {
            Resize(degreeU, degreeV, order);
        }

        void
            Resize(
                const int degreeU,
                const int degreeV,
                const int order)
        {
            const int nbShapes = NbShapes(order);
            const int nbNonzeroPoles = (degreeU + 1) * (degreeV + 1);

            m_shapeU.Resize(degreeU, order);
            m_shapeV.Resize(degreeV, order);
            m_values.resize(nbShapes * nbNonzeroPoles);
            m_weightedSums.resize(nbShapes);

            m_order = order;
        }

        int
            DegreeU() const
        {
            return m_shapeU.Degree();
        }

        int
            DegreeV() const
        {
            return m_shapeV.Degree();
        }

        int
            Order() const
        {
            return m_order;
        }

        int
            NbShapes() const
        {
            return NbShapes(Order());
        }

        int
            NbNonzeroPolesU() const
        {
            return m_shapeU.NbNonzeroPoles();
        }

        int
            NbNonzeroPolesV() const
        {
            return m_shapeV.NbNonzeroPoles();
        }

        int
            NbNonzeroPoles() const
        {
            return NbNonzeroPolesU() * NbNonzeroPolesV();
        }

        std::vector<std::pair<int, int>>
            NonzeroPoleIndices() const
        {
            std::vector<std::pair<int, int>> indices(NbNonzeroPoles());

            for (int i = 0; i < NbNonzeroPolesU(); i++) {
                for (int j = 0; j < NbNonzeroPolesV(); j++) {
                    int poleIndex = Math::MatrixIndex(NbNonzeroPolesU(),
                        NbNonzeroPolesV(), i, j);

                    int poleU = FirstNonzeroPoleU() + i;
                    int poleV = FirstNonzeroPoleV() + j;

                    indices[poleIndex] = { poleU, poleV };
                }
            }

            return indices;
        }

        const ScalarType
            Value(
                const int derivative,
                const int poleU,
                const int poleV) const
        {
            const int index = Index(derivative, poleU, poleV);

            return m_values[index];
        }

        const ScalarType
            Value(
                const int derivative,
                const int pole) const
        {
            const int index = Math::MatrixIndex(NbShapes(), NbNonzeroPoles(),
                derivative, pole);

            return m_values[index];
        }

        ScalarType
            operator()(
                const int derivative,
                const int pole) const
        {
            return Value(derivative, pole);
        }

        ScalarType
            operator()(
                const int derivative,
                const int poleU,
                const int poleV) const
        {
            return Value(derivative, poleU, poleV);
        }

        int
            FirstNonzeroPoleU() const
        {
            return m_firstNonzeroPoleU;
        }

        int
            LastNonzeroPoleU() const
        {
            return FirstNonzeroPoleU() + DegreeU();
        }

        int
            FirstNonzeroPoleV() const
        {
            return m_firstNonzeroPoleV;
        }

        int
            LastNonzeroPoleV() const
        {
            return FirstNonzeroPoleV() + DegreeV();
        }

        template <typename TKnots>
        void
            ComputeAtSpan(
                const TKnots& knotsU,
                const TKnots& knotsV,
                const int spanU,
                const int spanV,
                const ScalarType u,
                const ScalarType v)
        {
            const int nbValues = NbShapes() * NbNonzeroPoles();

            std::fill(m_values.begin(), m_values.begin() + nbValues, 0);

            m_firstNonzeroPoleU = spanU - DegreeU() + 1;
            m_firstNonzeroPoleV = spanV - DegreeV() + 1;

            // compute 1D shape functions

            m_shapeU.ComputeAtSpan(knotsU, spanU, u);
            m_shapeV.ComputeAtSpan(knotsV, spanV, v);

            // compute 2D shape functions

            for (int i = 0; i <= Order(); i++) {
                for (int j = 0; j <= Order() - i; j++) {
                    for (int a = 0; a < NbNonzeroPolesU(); a++) {
                        for (int b = 0; b < NbNonzeroPolesV(); b++) {
                            const int index = ShapeIndex(i, j);

                            Value(index, a, b) = m_shapeU(i, a) * m_shapeV(j, b);
                        }
                    }
                }
            }
        }

        template <typename TKnots>
        void
            Compute(
                const TKnots& knotsU,
                const TKnots& knotsV,
                const ScalarType u,
                const ScalarType v)
        {
            const int spanU = Knots::LowerSpan(DegreeU(), knotsU, u);
            const int spanV = Knots::LowerSpan(DegreeV(), knotsV, v);

            ComputeAtSpan(knotsU, knotsV, spanU, spanV, u, v);
        }

        template <typename TKnots, typename TWeights>
        void
            ComputeAtSpan(
                const TKnots& knotsU,
                const TKnots& knotsV,
                const int spanU,
                const int spanV,
                const TWeights& weights,
                const ScalarType u,
                const ScalarType v)
        {
            using Math::Binom;

            // compute B-Spline shape

            ComputeAtSpan(knotsU, knotsV, spanU, spanV, u, v);

            // apply weights

            for (int shape = 0; shape < NbShapes(); shape++) {
                WeightedSum(shape) = ScalarType(0);

                for (int i = 0; i < NbNonzeroPolesU(); i++) {
                    for (int j = 0; j < NbNonzeroPolesV(); j++) {
                        int poleU = FirstNonzeroPoleU() + i;
                        int poleV = FirstNonzeroPoleV() + j;

                        ScalarType weight = Util::SurfaceWeights<TWeights>::Get(
                            weights, poleU, poleV);
                        Value(shape, i, j) *= weight;
                        WeightedSum(shape) += Value(shape, i, j);
                    }
                }
            }

            for (int k = 0; k <= Order(); k++) {
                for (int l = 0; l <= Order() - k; l++) {
                    const int shape = ShapeIndex(k, l);

                    for (int j = 1; j <= l; j++) {
                        const int index = ShapeIndex(k, l - j);

                        ScalarType a = Binom(l, j) * WeightedSum(0, j);

                        for (int p = 0; p < NbNonzeroPoles(); p++) {
                            Value(shape, p) -= a * Value(index, p);
                        }
                    }

                    for (int i = 1; i <= k; i++) {
                        const int index = ShapeIndex(k - i, l);

                        ScalarType a = Binom(k, i) * WeightedSum(i, 0);

                        for (int p = 0; p < NbNonzeroPoles(); p++) {
                            Value(shape, p) -= a * Value(index, p);
                        }
                    }

                    for (int i = 1; i <= k; i++) {
                        const ScalarType a = Binom(k, i);

                        for (int j = 1; j <= l; j++) {
                            const int index = ShapeIndex(k - i, l - j);

                            const ScalarType b = a * Binom(l, j) *
                                WeightedSum(i, j);

                            for (int p = 0; p < NbNonzeroPoles(); p++) {
                                Value(shape, p) -= b * Value(index, p);
                            }
                        }
                    }

                    for (int p = 0; p < NbNonzeroPoles(); p++) {
                        Value(shape, p) /= WeightedSum(0);
                    }
                }
            }
        }

        template <typename TKnots, typename TWeights>
        void
            Compute(
                const TKnots& knotsU,
                const TKnots& knotsV,
                const TWeights& weights,
                const ScalarType u,
                const ScalarType v)
        {
            const int spanU = Knots::LowerSpan(DegreeU(), knotsU, u);
            const int spanV = Knots::LowerSpan(DegreeV(), knotsV, v);

            ComputeAtSpan(knotsU, knotsV, spanU, spanV, weights, u, v);
        }
    };

} // namespace ANurbs