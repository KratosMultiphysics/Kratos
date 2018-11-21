#pragma once

#include "Knots.h"
#include "Math.h"
#include "Util.h"
#include "VectorTraits.h"

#include <functional>
#include <vector>

namespace ANurbs {

template <typename TScalar>
class CurveShapeEvaluator
{
public:
    using ScalarType = TScalar;

private:
    int m_degree;
    int m_order;
    std::vector<ScalarType> m_values;
    std::vector<ScalarType> m_left;
    std::vector<ScalarType> m_right;
    std::vector<ScalarType> m_ndu;
    std::vector<ScalarType> m_a;
    std::vector<ScalarType> m_b;
    int m_firstNonzeroPole;

private:
    ScalarType&
    Values(
        const int order,
        const int pole)
    {
        int index = Math::MatrixIndex(NbShapes(), NbNonzeroPoles(),
            order, pole);

        return m_values[index];
    }

    ScalarType&
    Ndu(
        const int i,
        const int j)
    {
        int index = Math::MatrixIndex(NbShapes(), NbNonzeroPoles(),
            i, j);

        return m_ndu[index];
    }

    void
    ClearValues()
    {
        int nbValues = NbNonzeroPoles() * NbShapes();

        std::fill(m_values.begin(), m_values.begin() + nbValues, 0);
    }

public:
    CurveShapeEvaluator()
    {
    }

    CurveShapeEvaluator(
        const int degree,
        const int order)
    {
        Resize(degree, order);
    }

    void
    Resize(
        const int degree,
        const int order)
    {
        m_values.resize((order + 1) * (degree + 1));
        m_left.resize(degree);
        m_right.resize(degree);
        m_ndu.resize((degree + 1) * (degree + 1));
        m_a.resize(degree + 1);
        m_b.resize(degree + 1);

        m_degree = degree;
        m_order = order;
    }

    int
    Degree() const
    {
        return m_degree;
    }

    int
    Order() const
    {
        return m_order;
    }

    int
    NbNonzeroPoles() const
    {
        return Degree() + 1;
    }

    int
    NbShapes() const
    {
        return Order() + 1;
    }

    ScalarType
    operator()(
        const int order,
        const int pole) const
    {
        return Value(order, pole);
    }

    ScalarType
    Value(
        const int order,
        const int pole) const
    {
        int index = Math::MatrixIndex(NbShapes(), NbNonzeroPoles(),
            order, pole);

        return m_values[index];
    }

    int
    FirstNonzeroPole() const
    {
        return m_firstNonzeroPole;
    }

    int
    LastNonzeroPole() const
    {
        return FirstNonzeroPole() + Degree();
    }

    std::vector<int>
    NonZeroPoleIndices() const
    {
        std::vector<int> indices(NbNonzeroPoles());

        for (int i = 0; i < NbNonzeroPoles(); i++) {
            indices[i] = FirstNonzeroPole() + i;
        }

        return indices;
    }

    template <typename Knots>
    void
    ComputeAtSpan(
        const Knots& knots,
        const int span,
        const ScalarType t)
    {
        ClearValues();

        m_firstNonzeroPole = span - Degree() + 1;

        // compute B-Spline shape

        Ndu(0, 0) = 1.0;

        for (int j = 0; j < Degree(); j++) {
            m_left[j] = t - knots[span - j];
            m_right[j] = knots[span + j + 1] - t;

            ScalarType saved = 0.0;

            for (int r = 0; r <= j; r++) {
                Ndu(j + 1, r) = m_right[r] + m_left[j - r];

                ScalarType temp = Ndu(r, j) / Ndu(j + 1, r);

                Ndu(r, j + 1) = saved + m_right[r] * temp;

                saved = m_left[j - r] * temp;
            }

            Ndu(j + 1, j + 1) = saved;
        }

        for (int j = 0; j < NbNonzeroPoles(); j++) {
            Values(0, j) = Ndu(j, Degree());
        }

        auto& a = m_a;
        auto& b = m_b;

        for (int r = 0; r < NbNonzeroPoles(); r++) {
            a[0] = 1.0;

            for (int k = 1; k < NbShapes(); k++) {
                ScalarType& value = Values(k, r);

                int rk = r - k;
                int pk = Degree() - k;

                if (r >= k) {
                    b[0] = a[0] / Ndu(pk + 1, rk);
                    value = b[0] * Ndu(rk, pk);
                }

                int j1 = r >= k - 1 ? 1 : k - r;
                int j2 = r <= pk + 1 ? k : NbNonzeroPoles() - r;

                for (int j = j1; j < j2; j++) {
                    b[j] = (a[j] - a[j - 1]) / Ndu(pk + 1, rk + j);
                    value += b[j] * Ndu(rk + j, pk);
                }

                if (r <= pk) {
                    b[k] = -a[k - 1] / Ndu(pk + 1, r);
                    value += b[k] * Ndu(r, pk);
                }

                std::swap(a, b);
            }
        }

        int s = Degree();

        for (int k = 1; k < NbShapes(); k++) {
            for (int j = 0; j < NbNonzeroPoles(); j++) {
                Values(k, j) *= s;
            }
            s *= Degree() - k;
        }
    }

    template <typename TKnots, typename TWeights>
    void
    ComputeAtSpan(
        const TKnots& knots,
        const int span,
        const TWeights& weights,
        const ScalarType t)
    {
        // compute B-Spline shape

        ComputeAtSpan(knots, span, t);

        // compute weighted sum

        ScalarType weightedSum{ 0 };

        for (int i = 0; i < NbNonzeroPoles(); i++) {
            m_values[i] *= Util::CurveWeights<TWeights>::Get(weights, i);
            weightedSum += m_values[i];
        }

        // apply weights

        for (int i = 0; i < NbNonzeroPoles(); i++) {
            m_values[i] /= weightedSum;
        }
    }

    template <typename TKnots>
    void
    Compute(
        const TKnots& knots,
        const ScalarType t)
    {
        int span = Knots::UpperSpan(Degree(), knots, t);

        ComputeAtSpan(knots, span, t);
    }

    template <typename TKnots, typename TWeights>
    void
    Compute(
        const TKnots& knots,
        const TWeights& weights,
        const ScalarType t)
    {
        int span = Knots::UpperSpan(Degree(), knots, t);

        ComputeAtSpan(knots, span, weights, t);
    }
};

} // namespace ANurbs
