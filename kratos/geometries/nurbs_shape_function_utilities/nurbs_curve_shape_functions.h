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

    typedef typename std::size_t IndexType;
    typedef typename std::size_t SizeType;

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
        const SizeType PolynomialDegree,
        const SizeType DerivativeOrder)
    {
        ResizeDataContainers(PolynomialDegree, DerivativeOrder);
    }

    ///@}
    ///@name Operators
    ///@{

    double operator()(const SizeType DerivativeOrder, const IndexType ControlPoint) const
    {
        return Value(DerivativeOrder, ControlPoint);
    }

    ///@}
    ///@name Operations
    ///@{

    /* Resizes the data containers which are needed to compute the values of the
    shape functions. */
    void ResizeDataContainers(const SizeType PolynomialDegree, const SizeType DerivativeOrder)
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

    /* @return provides the polynomial degree of this shape function generator. */
    SizeType PolynomialDegree() const
    {
        return mPolynomialDegree;
    }

    /* @return the number of nonzero control points of this shape function generator. */
    SizeType NumberOfNonzeroControlPoints() const
    {
        return PolynomialDegree() + 1;
    }

    /* @return the number of shape function rows. This is the derivative order + 1.
       rows are defined as: N | dN/de | dN^2/de^2 | ...*/
    SizeType NumberOfShapeFunctionRows() const
    {
        return mDerivativeOrder + 1;
    }

    /* Provides the shape function value depending to the DerivativeOrder and
    the index of the control point.*/
    double Value(const SizeType DerivativeOrder, const IndexType ControlPoint) const
    {
        IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(NumberOfShapeFunctionRows(),
            NumberOfNonzeroControlPoints(), DerivativeOrder, ControlPoint);

        return mValues[index];
    }

    /* @return the index of the first nonzero control point.*/
    IndexType GetFirstNonzeroControlPoint() const
    {
        return mFirstNonzeroControlPoint;
    }

    /* @return the index of the last nonzero control point.*/
    IndexType GetLastNonzeroControlPoint() const
    {
        return GetFirstNonzeroControlPoint() + PolynomialDegree();
    }

    /* @return the indices of all nonzero control points.*/
    std::vector<IndexType> GetNonzeroControlPointIndices() const
    {
        std::vector<IndexType> indices(NumberOfNonzeroControlPoints());

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
        const IndexType Span,
        const double ParameterT)
    {
        ClearValues();

        mFirstNonzeroControlPoint = Span - PolynomialDegree() + 1;

        // compute B-Spline shape

        Ndu(0, 0) = 1.0;

        for (IndexType j = 0; j < PolynomialDegree(); j++) {
            mLeft[j] = ParameterT - rKnots[Span - j];
            mRight[j] = rKnots[Span + j + 1] - ParameterT;

            double saved = 0.0;

            for (IndexType r = 0; r <= j; r++) {
                Ndu(j + 1, r) = mRight[r] + mLeft[j - r];

                double temp = Ndu(r, j) / Ndu(j + 1, r);

                Ndu(r, j + 1) = saved + mRight[r] * temp;

                saved = mLeft[j - r] * temp;
            }

            Ndu(j + 1, j + 1) = saved;
        }

        for (IndexType j = 0; j < NumberOfNonzeroControlPoints(); j++) {
            Value(0, j) = Ndu(j, PolynomialDegree());
        }

        auto& a = mA;
        auto& b = mB;

        for (IndexType r = 0; r < NumberOfNonzeroControlPoints(); r++) {
            a[0] = 1.0;

            for (IndexType k = 1; k < NumberOfShapeFunctionRows(); k++) {
                double& value = Value(k, r);

                int rk = r - k;
                int pk = PolynomialDegree() - k;

                if (r >= k) {
                    b[0] = a[0] / Ndu(pk + 1, rk);
                    value = b[0] * Ndu(rk, pk);
                }

                int j1 = r >= k - 1 ? 1 : k - r;
                int j2 = r <= pk + 1 ? k : NumberOfNonzeroControlPoints() - r;

                for (IndexType j = j1; j < j2; j++) {
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

        for (IndexType k = 1; k < NumberOfShapeFunctionRows(); k++) {
            for (IndexType j = 0; j < NumberOfNonzeroControlPoints(); j++) {
                Value(k, j) *= s;
            }
            s *= PolynomialDegree() - k;
        }
    }

    void ComputeNurbsShapeFunctionValuesAtSpan(
        const Vector& rKnots,
        const IndexType Span,
        const Vector& rWeights,
        const double ParameterT)
    {
        // compute B-Spline shape
        ComputeBSplineShapeFunctionValuesAtSpan(
            rKnots, Span, ParameterT);

        // compute weighted sum
        double weightedSum{ 0 };

        for (IndexType i = 0; i < NumberOfNonzeroControlPoints(); i++) {
            mValues[i] *= rWeights(i);
            weightedSum += mValues[i];
        }

        // apply weights

        for (IndexType i = 0; i < NumberOfNonzeroControlPoints(); i++) {
            mValues[i] /= weightedSum;
        }
    }

    void ComputeBSplineShapeFunctionValues(const Vector& rKnots, const double ParameterT)
    {
        const IndexType span = NurbsUtilities::GetUpperSpan(PolynomialDegree(), rKnots,
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

    /* Provides the shape function value depending to the DerivativeOrder and
    the index of the control point.*/
    double& Value(const SizeType DerivativeOrder, const IndexType ControlPointIndex)
    {
        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(NumberOfShapeFunctionRows(),
            NumberOfNonzeroControlPoints(), DerivativeOrder, ControlPointIndex);

        return mValues[index];
    }

    double& Ndu(const IndexType IndexI, const IndexType IndexJ)
    {
        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(NumberOfShapeFunctionRows(),
            NumberOfNonzeroControlPoints(), IndexI, IndexJ);

        return mNdu[index];
    }

    void ClearValues()
    {
        const SizeType nb_values = NumberOfNonzeroControlPoints() * NumberOfShapeFunctionRows();

        std::fill(mValues.begin(), mValues.begin() + nb_values, 0);
    }

    ///@}
    ///@name Member Variables
    ///@{

    SizeType mPolynomialDegree;
    SizeType mDerivativeOrder;
    Vector mValues;
    Vector mLeft;
    Vector mRight;
    Vector mNdu;
    Vector mA;
    Vector mB;
    IndexType mFirstNonzeroControlPoint;

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