//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_CURVE_SHAPE_FUNCTIONS_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_SHAPE_FUNCTIONS_H_INCLUDED

#include "nurbs_utility.h"

#include <vector>

namespace Kratos {

class NurbsCurveShapeFunction
{
private:    // variables
    int mDegree;
    int mOrder;
    std::vector<double> mValues;
    std::vector<double> mLeft;
    std::vector<double> mRight;
    std::vector<double> mNdu;
    std::vector<double> mA;
    std::vector<double> mB;
    int mFirstNonzeroPole;

private:    // methods
    double& Value(const int Order, const int Pole)
    {
        const int index = NurbsUtility::GetSingleIndex(GetNbShapeFunctions(),
            GetNbNonzeroPoles(), Order, Pole);

        return mValues[index];
    }

    double& Ndu(const int IndexI, const int IndexJ)
    {
        const int index = NurbsUtility::GetSingleIndex(GetNbShapeFunctions(),
            GetNbNonzeroPoles(), IndexI, IndexJ);

        return mNdu[index];
    }

    void ClearValues()
    {
        const int nb_values = GetNbNonzeroPoles() * GetNbShapeFunctions();

        std::fill(mValues.begin(), mValues.begin() + nb_values, 0);
    }

public:     // constructors
    NurbsCurveShapeFunction()
    {
    }

    NurbsCurveShapeFunction(const int Degree, const int Order)
    {
        Resize(Degree, Order);
    }

public:     // methods
    void Resize(const int Degree, const int Order)
    {
        mValues.resize((Order + 1) * (Degree + 1));
        mLeft.resize(Degree);
        mRight.resize(Degree);
        mNdu.resize((Degree + 1) * (Degree + 1));
        mA.resize(Degree + 1);
        mB.resize(Degree + 1);

        mDegree = Degree;
        mOrder = Order;
    }

    int GetDegree() const
    {
        return mDegree;
    }

    int GetOrder() const
    {
        return mOrder;
    }

    int GetNbNonzeroPoles() const
    {
        return GetDegree() + 1;
    }

    int GetNbShapeFunctions() const
    {
        return GetOrder() + 1;
    }

    double operator()(const int Order, const int Pole) const
    {
        return Value(Order, Pole);
    }

    double Value(const int Order, const int Pole) const
    {
        int index = NurbsUtility::GetSingleIndex(GetNbShapeFunctions(),
            GetNbNonzeroPoles(), Order, Pole);

        return mValues[index];
    }

    int GetFirstNonzeroPole() const
    {
        return mFirstNonzeroPole;
    }

    int GetLastNonzeroPole() const
    {
        return GetFirstNonzeroPole() + GetDegree();
    }

    std::vector<int> GetNonzeroPoleIndices() const
    {
        std::vector<int> indices(GetNbNonzeroPoles());

        for (int i = 0; i < GetNbNonzeroPoles(); i++) {
            indices[i] = GetFirstNonzeroPole() + i;
        }

        return indices;
    }

    template <typename TKnots>
    void ComputeAtSpan(const TKnots& rKnots, const int Span,
        const double ParameterT)
    {
        ClearValues();

        mFirstNonzeroPole = Span - GetDegree() + 1;

        // compute B-Spline shape

        Ndu(0, 0) = 1.0;

        for (int j = 0; j < GetDegree(); j++) {
            mLeft[j] = ParameterT - rKnots[Span - j];
            mRight[j] = rKnots[Span + j + 1] - ParameterT;

            double saved = 0.0;

            for (int r = 0; r <= j; r++) {
                Ndu(j + 1, r) = mRight[r] + mLeft[j - r];

                double temp = Ndu(r, j) / Ndu(j + 1, r);

                Ndu(r, j + 1) = saved + mRight[r] * temp;

                saved = mLeft[j - r] * temp;
            }

            Ndu(j + 1, j + 1) = saved;
        }

        for (int j = 0; j < GetNbNonzeroPoles(); j++) {
            Value(0, j) = Ndu(j, GetDegree());
        }

        auto& a = mA;
        auto& b = mB;

        for (int r = 0; r < GetNbNonzeroPoles(); r++) {
            a[0] = 1.0;

            for (int k = 1; k < GetNbShapeFunctions(); k++) {
                double& value = Value(k, r);

                int rk = r - k;
                int pk = GetDegree() - k;

                if (r >= k) {
                    b[0] = a[0] / Ndu(pk + 1, rk);
                    value = b[0] * Ndu(rk, pk);
                }

                int j1 = r >= k - 1 ? 1 : k - r;
                int j2 = r <= pk + 1 ? k : GetNbNonzeroPoles() - r;

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

        int s = GetDegree();

        for (int k = 1; k < GetNbShapeFunctions(); k++) {
            for (int j = 0; j < GetNbNonzeroPoles(); j++) {
                Value(k, j) *= s;
            }
            s *= GetDegree() - k;
        }
    }

    template <typename TKnots, typename TWeights>
    void ComputeAtSpan(const TKnots& rKnots, const int Span,
        const TWeights& rWeights, const double ParameterT)
    {
        // compute B-Spline shape

        ComputeAtSpan(rKnots, Span, ParameterT);

        // compute weighted sum

        double weightedSum {0};

        for (int i = 0; i < GetNbNonzeroPoles(); i++) {
            mValues[i] *= rWeights(i);
            weightedSum += mValues[i];
        }

        // apply weights

        for (int i = 0; i < GetNbNonzeroPoles(); i++) {
            mValues[i] /= weightedSum;
        }
    }

    void Compute(const std::vector<double>& rKnots, const double ParameterT)
    {
        const int span = NurbsUtility::GetUpperSpan(GetDegree(), rKnots,
            ParameterT);

        ComputeAtSpan(rKnots, span, ParameterT);
    }

    template <typename TWeights>
    void Compute(const std::vector<double>& rKnots, const TWeights& rWeights,
        const double ParameterT)
    {
        const int span = NurbsUtility::GetUpperSpan(GetDegree(), rKnots,
            ParameterT);

        ComputeAtSpan(rKnots, span, rWeights, ParameterT);
    }
};

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_SHAPE_FUNCTIONS_H_INCLUDED defined