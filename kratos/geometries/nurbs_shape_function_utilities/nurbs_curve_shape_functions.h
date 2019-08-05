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
//                   Tobias Teschemacher
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_CURVE_SHAPE_FUNCTIONS_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_SHAPE_FUNCTIONS_H_INCLUDED

#include "nurbs_utilities.h"

#include <vector>

namespace Kratos {

class NurbsCurveShapeFunction
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of NurbsSurfaceShapeFunction
    //KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NurbsCurveShapeFunction);
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NurbsCurveShapeFunction()
    {
    };

    /* Constructor using the degree u, degree v and the order of shape functions.
       This is required to make an optimized memory management. */
    NurbsCurveShapeFunction(
        const int PolynomialDegree,
        const int DerivativeOrder)
    {
        ResizeDataContainers(PolynomialDegree, DerivativeOrder);
    }

    ///@}
    ///@name Operators
    ///@{

    double operator()(const int DerivativeOrder, const int ControlPoint) const
    {
        return Value(DerivativeOrder, ControlPoint);
    }

    ///@}
    ///@name Operations
    ///@{

    /* Resizes the data containers which are needed to compute the values of the
    shape functions. */
    void ResizeDataContainers(const int PolynomialDegree, const int DerivativeOrder)
    {
        mValues.resize((DerivativeOrder + 1) * (PolynomialDegree + 1));
        mLeft.resize(PolynomialDegree);
        mRight.resize(PolynomialDegree);
        mNdu.resize((PolynomialDegree + 1) * (PolynomialDegree + 1));
        mA.resize(PolynomialDegree + 1);
        mB.resize(PolynomialDegree + 1);

        mPolynomialDegree = PolynomialDegree;
        mDerivativeOrder = DerivativeOrder;
    }

    int PolynomialDegree() const
    {
        return mPolynomialDegree;
    }

    int NumberOfNonzeroControlPoints() const
    {
        return PolynomialDegree() + 1;
    }

    int NumberOfShapeFunctionRows() const
    {
        return mDerivativeOrder + 1;
    }

    double Value(const int DerivativeOrder, const int ControlPoint) const
    {
        int index = NurbsUtilities::GetVectorIndexFromMatrixIndices(NumberOfShapeFunctionRows(),
            NumberOfNonzeroControlPoints(), DerivativeOrder, ControlPoint);

        return mValues[index];
    }

    int GetFirstNonzeroControlPoint() const
    {
        return mFirstNonzeroControlPoint;
    }

    int GetLastNonzeroControlPoint() const
    {
        return GetFirstNonzeroControlPoint() + PolynomialDegree();
    }

    std::vector<int> GetNonzeroControlPointIndices() const
    {
        std::vector<int> indices(NumberOfNonzeroControlPoints());

        for (int i = 0; i < NumberOfNonzeroControlPoints(); i++) {
            indices[i] = GetFirstNonzeroControlPoint() + i;
        }

        return indices;
    }

    ///@}
    ///@name Shape Function Computation
    ///@{

    void ComputeBSplineShapeFunctionValuesAtSpan(
        const Vector& rKnots,
        const int Span,
        const double ParameterT)
    {
        ClearValues();

        mFirstNonzeroControlPoint = Span - PolynomialDegree() + 1;

        // compute B-Spline shape

        Ndu(0, 0) = 1.0;

        for (int j = 0; j < PolynomialDegree(); j++) {
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

        for (int j = 0; j < NumberOfNonzeroControlPoints(); j++) {
            Value(0, j) = Ndu(j, PolynomialDegree());
        }

        auto& a = mA;
        auto& b = mB;

        for (int r = 0; r < NumberOfNonzeroControlPoints(); r++) {
            a[0] = 1.0;

            for (int k = 1; k < NumberOfShapeFunctionRows(); k++) {
                double& value = Value(k, r);

                int rk = r - k;
                int pk = PolynomialDegree() - k;

                if (r >= k) {
                    b[0] = a[0] / Ndu(pk + 1, rk);
                    value = b[0] * Ndu(rk, pk);
                }

                int j1 = r >= k - 1 ? 1 : k - r;
                int j2 = r <= pk + 1 ? k : NumberOfNonzeroControlPoints() - r;

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

        int s = PolynomialDegree();

        for (int k = 1; k < NumberOfShapeFunctionRows(); k++) {
            for (int j = 0; j < NumberOfNonzeroControlPoints(); j++) {
                Value(k, j) *= s;
            }
            s *= PolynomialDegree() - k;
        }
    }

    void ComputeNurbsShapeFunctionValuesAtSpan(
        const Vector& rKnots,
        const int Span,
        const Vector& rWeights,
        const double ParameterT)
    {
        // compute B-Spline shape
        ComputeBSplineShapeFunctionValuesAtSpan(
            rKnots, Span, ParameterT);

        // compute weighted sum
        double weightedSum{ 0 };

        for (int i = 0; i < NumberOfNonzeroControlPoints(); i++) {
            mValues[i] *= rWeights(i);
            weightedSum += mValues[i];
        }

        // apply weights

        for (int i = 0; i < NumberOfNonzeroControlPoints(); i++) {
            mValues[i] /= weightedSum;
        }
    }

    void ComputeBSplineShapeFunctionValues(const Vector& rKnots, const double ParameterT)
    {
        const int span = NurbsUtilities::GetUpperSpan(PolynomialDegree(), rKnots,
            ParameterT);

        ComputeBSplineShapeFunctionValuesAtSpan(rKnots, span, ParameterT);
    }

    void ComputeNurbsShapeFunctionValues(
        const Vector& rKnots,
        const Vector& rWeights,
        const double ParameterT)
    {
        const int span = NurbsUtilities::GetUpperSpan(
            PolynomialDegree(), rKnots, ParameterT);

        ComputeNurbsShapeFunctionValuesAtSpan(rKnots, span, rWeights, ParameterT);
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    double& Value(const int DerivativeOrder, const int ControlPointIndex)
    {
        const int index = NurbsUtilities::GetVectorIndexFromMatrixIndices(NumberOfShapeFunctionRows(),
            NumberOfNonzeroControlPoints(), DerivativeOrder, ControlPointIndex);

        return mValues[index];
    }

    double& Ndu(const int IndexI, const int IndexJ)
    {
        const int index = NurbsUtilities::GetVectorIndexFromMatrixIndices(NumberOfShapeFunctionRows(),
            NumberOfNonzeroControlPoints(), IndexI, IndexJ);

        return mNdu[index];
    }

    void ClearValues()
    {
        const int nb_values = NumberOfNonzeroControlPoints() * NumberOfShapeFunctionRows();

        std::fill(mValues.begin(), mValues.begin() + nb_values, 0);
    }

    ///@}
    ///@name Member Variables
    ///@{

    int mPolynomialDegree;
    int mDerivativeOrder;
    Vector mValues;
    Vector mLeft;
    Vector mRight;
    Vector mNdu;
    Vector mA;
    Vector mB;
    int mFirstNonzeroControlPoint;

    ///@}
    ///@name Private Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("PolynomialDegree", mPolynomialDegree);
        rSerializer.save("DerivativeOrder", mDerivativeOrder);
        rSerializer.save("Values", mValues);
        rSerializer.save("Left", mLeft);
        rSerializer.save("Right", mRight);
        rSerializer.save("Ndu", mNdu);
        rSerializer.save("A", mA);
        rSerializer.save("B", mB);
        rSerializer.save("FirstNonzeroControlPoint", mFirstNonzeroControlPoint);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("PolynomialDegree", mPolynomialDegree);
        rSerializer.load("DerivativeOrder", mDerivativeOrder);
        rSerializer.load("Values", mValues);
        rSerializer.load("Left", mLeft);
        rSerializer.load("Right", mRight);
        rSerializer.load("Ndu", mNdu);
        rSerializer.load("A", mA);
        rSerializer.load("B", mB);
        rSerializer.load("FirstNonzeroControlPoint", mFirstNonzeroControlPoint);
    }

    ///@}
};

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_SHAPE_FUNCTIONS_H_INCLUDED defined