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

#if !defined(KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED )
#define  KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED

// Project includes
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"

namespace Kratos 
{
///@name Kratos Classes
///@{
/// Computes the shape function for 2D geometrical objects as surfaces.
/** This class is to be initialized to optimize. It creates the containers for the 
shape functions and the derivatives for an optimized data treatment.
*/
class NurbsSurfaceShapeFunction
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
    NurbsSurfaceShapeFunction()
    {
    };

    // Constructor using the degree u, degree v and the order of shape functions.
    // This is required to make an optimized memory management.
    NurbsSurfaceShapeFunction(
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const SizeType DerivativeOrder)
    {
        ResizeDataContainers(PolynomialDegreeU, PolynomialDegreeV, DerivativeOrder);
    }

    ///@}
    ///@name Static Operations
    ///@{
     /* @brief Returns the number of shape function rows for a given order.
        @param DerivativesOrder 0 the shape functions N only, 1 N and first derivatives, ...
               Thus, for DerivativesOrder 0 the return value is 1, for 1 the return value is 3, ...
        @return NumberOfShapeFunctionRows, shape functions are provided as:
                N | dN/du, dN/dv | dN^2/du^2, dN^2/dn*dv, dN^2/dv^2 | ... */
    static constexpr inline SizeType NumberOfShapeFunctionRows(const SizeType DerivativeOrder) noexcept
    {
        return (1 + DerivativeOrder) * (2 + DerivativeOrder) / 2;
    }

    /* @brief Returns the index of the shape function row for a given derivative index.
       @param DerivativesOrder 0 the shape functions N only, 1 N and first derivatives, ...
              Thus, for DerivativesOrder 0 the return value is 1, for 1 the return value is 3, ...
       @return NumberOfShapeFunctionRows, shape functions are provided as:
               N | dN/du, dN/dv | dN^2/du^2, dN^2/dn*dv, dN^2/dv^2 | ... */
    static constexpr inline IndexType IndexOfShapeFunctionRow(
        const SizeType DerivativeOrderU,
        const SizeType DerivativeOrderV) noexcept
    {
        return DerivativeOrderV + (DerivativeOrderU + DerivativeOrderV) * 
            (1 + DerivativeOrderU + DerivativeOrderV) / 2;
    }
    ///@}
    ///@name Operators
    ///@{

    double operator()(
        const IndexType DerivativeRow,
        const IndexType ControlPointIndex) const
    {
        return ShapeFunctionValue(DerivativeRow, ControlPointIndex);
    }

    double operator()(
        const IndexType DerivativeRow,
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV) const
    {
        return ShapeFunctionValue(DerivativeRow, ControlPointIndexU, ControlPointIndexV);
    }

    ///@}
    ///@name Operations
    ///@{
    void ResizeDataContainers(
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const SizeType DerivativeOrder)
    {
        const SizeType number_of_shape_function_rows = this->NumberOfShapeFunctionRows(DerivativeOrder);
        const SizeType number_of_nonzero_control_points = (PolynomialDegreeU + 1) * (PolynomialDegreeV + 1);

        mShapeFunctionsU.ResizeDataContainers(PolynomialDegreeU, std::min(DerivativeOrder, PolynomialDegreeU));
        mShapeFunctionsV.ResizeDataContainers(PolynomialDegreeV, std::min(DerivativeOrder, PolynomialDegreeV));
        mShapeFunctionValues.resize(number_of_shape_function_rows * number_of_nonzero_control_points);
        mWeightedSums.resize(number_of_shape_function_rows);

        mDerivativeOrder = DerivativeOrder;
    }

    SizeType PolynomialDegreeU() const
    {
        return mShapeFunctionsU.PolynomialDegree();
    }

    SizeType PolynomialDegreeV() const
    {
        return mShapeFunctionsV.PolynomialDegree();
    }

    SizeType DerivativeOrder() const
    {
        return mDerivativeOrder;
    }

    SizeType NumberOfShapeFunctionRows() const
    {
        return NumberOfShapeFunctionRows(DerivativeOrder());
    }

    SizeType NumberOfNonzeroControlPointsU() const
    {
        return mShapeFunctionsU.NumberOfNonzeroControlPoints();
    }

    SizeType NumberOfNonzeroControlPointsV() const
    {
        return mShapeFunctionsV.NumberOfNonzeroControlPoints();
    }

    SizeType NumberOfNonzeroControlPoints() const
    {
        return NumberOfNonzeroControlPointsU() * NumberOfNonzeroControlPointsV();
    }

    std::vector<std::pair<int, int>> NumberOfNonzeroControlPointIndices() const
    {
        std::vector<std::pair<int, int>> indices(NumberOfNonzeroControlPoints());

        for (int i = 0; i < NumberOfNonzeroControlPointsU(); i++) {
            for (int j = 0; j < NumberOfNonzeroControlPointsV(); j++) {
                int poleIndex = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                    NumberOfNonzeroControlPointsU(), NumberOfNonzeroControlPointsV(), i, j);

                int poleU = GetFirstNonzeroControlPointU() + i;
                int poleV = GetFirstNonzeroControlPointV() + j;

                indices[poleIndex] = { poleU, poleV };
            }
        }

        return indices;
    }

    const double ShapeFunctionValue(
        const SizeType DerivativeRow,
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV) const
    {
        const IndexType index = this->GetIndex(
            DerivativeRow,
            ControlPointIndexU,
            ControlPointIndexV);

        return mShapeFunctionValues[index];
    }

    const double ShapeFunctionValue(
        const SizeType DerivativeOrder,
        const IndexType ControlPointIndex) const
    {
        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NumberOfShapeFunctionRows(), NumberOfNonzeroControlPoints(),
            DerivativeOrder, ControlPointIndex);

        return mShapeFunctionValues[index];
    }

    IndexType GetFirstNonzeroControlPointU() const
    {
        return mFirstNonzeroControlPointU;
    }

    IndexType GetLastNonzeroControlPointU() const
    {
        return GetFirstNonzeroControlPointU() + PolynomialDegreeU();
    }

    IndexType GetFirstNonzeroControlPointV() const
    {
        return mFirstNonzeroControlPointV;
    }

    IndexType GetLastNonzeroControlPointV() const
    {
        return GetFirstNonzeroControlPointV() + PolynomialDegreeV();
    }

    void ComputeBSplineShapeFunctionValuesAtSpan(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const int SpanU,
        const int SpanV,
        const double ParameterU,
        const double ParameterV)
    {
        const int number_of_values = NumberOfShapeFunctionRows() * NumberOfNonzeroControlPoints();

        mShapeFunctionValues = ZeroVector(mShapeFunctionValues.size());

        mFirstNonzeroControlPointU = SpanU - PolynomialDegreeU() + 1;
        mFirstNonzeroControlPointV = SpanV - PolynomialDegreeV() + 1;

        // compute 1D shape functions
        mShapeFunctionsU.ComputeBSplineShapeFunctionValuesAtSpan(rKnotsU, SpanU, ParameterU);
        mShapeFunctionsV.ComputeBSplineShapeFunctionValuesAtSpan(rKnotsV, SpanV, ParameterV);

        // compute 2D shape functions
        for (int i = 0; i <= DerivativeOrder(); i++) {
            for (int j = 0; j <= DerivativeOrder() - i; j++) {
                for (int a = 0; a < NumberOfNonzeroControlPointsU(); a++) {
                    for (int b = 0; b < NumberOfNonzeroControlPointsV(); b++) {
                        const int index = IndexOfShapeFunctionRow(i, j);

                        ShapeFunctionValue(index, a, b) = mShapeFunctionsU(i, a) * mShapeFunctionsV(j, b);
                    }
                }
            }
        }
    }

    void ComputeBSplineShapeFunctionValues(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const double ParameterU,
        const double ParameterV)
    {
        const IndexType SpanU = NurbsUtilities::GetLowerSpan(PolynomialDegreeU(), rKnotsU, ParameterU);
        const IndexType SpanV = NurbsUtilities::GetLowerSpan(PolynomialDegreeV(), rKnotsV, ParameterV);

        ComputeBSplineShapeFunctionValuesAtSpan(
            rKnotsU,
            rKnotsV,
            SpanU,
            SpanV,
            ParameterU,
            ParameterV);
    }

    void ComputeNurbsShapeFunctionValuesAtSpan(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const IndexType SpanU,
        const IndexType SpanV,
        const Vector& Weights,
        const double ParameterU,
        const double ParameterV)
    {
        // compute B-Spline shape functions
        ComputeBSplineShapeFunctionValuesAtSpan(
            rKnotsU, rKnotsV, SpanU, SpanV, ParameterU, ParameterV);

        // apply weights
        for (int shape_row_index = 0; shape_row_index < NumberOfShapeFunctionRows(); shape_row_index++) {
            GetWeightedSum(shape_row_index) = double(0);

            for (int u = 0; u < NumberOfNonzeroControlPointsU(); u++) {
                for (int v = 0; v < NumberOfNonzeroControlPointsV(); v++) {
                    const IndexType ControlPointIndexU = GetFirstNonzeroControlPointU() + u;
                    const IndexType ControlPointIndexV = GetFirstNonzeroControlPointV() + v;

                    const double weight = Weights(GetControlPointIndex(rKnotsU.size(), rKnotsV.size(), ControlPointIndexU, ControlPointIndexV));
                    ShapeFunctionValue(shape_row_index, u, v) *= weight;
                    GetWeightedSum(shape_row_index) += ShapeFunctionValue(shape_row_index, u, v);
                }
            }
        }

        for (int k = 0; k <= DerivativeOrder(); k++) {
            for (int l = 0; l <= DerivativeOrder() - k; l++) {
                const int shape = IndexOfShapeFunctionRow(k, l);

                for (int j = 1; j <= l; j++) {
                    const int index = IndexOfShapeFunctionRow(k, l - j);

                    double a = NurbsUtilities::GetBinomCoefficient(l, j) * GetWeightedSum(0, j);

                    for (int p = 0; p < NumberOfNonzeroControlPoints(); p++) {
                        ShapeFunctionValue(shape, p) -= a * ShapeFunctionValue(index, p);
                    }
                }

                for (int i = 1; i <= k; i++) {
                    const int index = IndexOfShapeFunctionRow(k - i, l);

                    double a = NurbsUtilities::GetBinomCoefficient(k, i) * GetWeightedSum(i, 0);

                    for (int p = 0; p < NumberOfNonzeroControlPoints(); p++) {
                        ShapeFunctionValue(shape, p) -= a * ShapeFunctionValue(index, p);
                    }
                }

                for (int i = 1; i <= k; i++) {
                    const double a = NurbsUtilities::GetBinomCoefficient(k, i);

                    for (int j = 1; j <= l; j++) {
                        const int index = IndexOfShapeFunctionRow(k - i, l - j);

                        const double b = a * NurbsUtilities::GetBinomCoefficient(l, j) *
                            GetWeightedSum(i, j);

                        for (int p = 0; p < NumberOfNonzeroControlPoints(); p++) {
                            ShapeFunctionValue(shape, p) -= b * ShapeFunctionValue(index, p);
                        }
                    }
                }

                for (int p = 0; p < NumberOfNonzeroControlPoints(); p++) {
                    ShapeFunctionValue(shape, p) /= GetWeightedSum(0);
                }
            }
        }
    }

    void ComputeNurbsShapeFunctionValues(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& Weights,
        const double ParameterU,
        const double ParameterV)
    {
        const IndexType SpanU = NurbsUtilities::GetLowerSpan(PolynomialDegreeU(), rKnotsU, ParameterU);
        const IndexType SpanV = NurbsUtilities::GetLowerSpan(PolynomialDegreeV(), rKnotsV, ParameterV);

        ComputeNurbsShapeFunctionValuesAtSpan(
            rKnotsU, rKnotsV, SpanU, SpanV, Weights, ParameterU, ParameterV);
    }
    ///@}
private:
    ///@name Operations
    ///@{
    double& GetWeightedSum(const IndexType Index)
    {
        return mWeightedSums[Index];
    }

    double& GetWeightedSum(const SizeType DerivativeOrderU, const SizeType DerivativeOrderV)
    {
        const int index = IndexOfShapeFunctionRow(DerivativeOrderU, DerivativeOrderV);

        return GetWeightedSum(index);
    }

    inline int GetControlPointIndex(
        const SizeType NumberOfKnotsU,
        const SizeType NumberOfKnotsV,
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV)
        const
    {
        return NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NurbsUtilities::GetNumberOfControlPoints(PolynomialDegreeU(), NumberOfKnotsU),
            NurbsUtilities::GetNumberOfControlPoints(PolynomialDegreeV(), NumberOfKnotsV),
            ControlPointIndexU, ControlPointIndexV);
    }

    inline int GetNonzeroControlPointIndex(
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV)
        const
    {
        return NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NumberOfNonzeroControlPointsU(), NumberOfNonzeroControlPointsV(),
            ControlPointIndexU, ControlPointIndexV);
    }

    inline int GetIndex(
        const SizeType DerivativeRow,
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV)
        const
    {
        const IndexType control_point_index = GetNonzeroControlPointIndex(
            ControlPointIndexU, ControlPointIndexV);
        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NumberOfShapeFunctionRows(), NumberOfNonzeroControlPoints(),
            DerivativeRow, control_point_index);

        return index;
    }

    double& ShapeFunctionValue(
        const IndexType DerivativeRow,
        const IndexType ControlPointIndex)
    {
        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NumberOfShapeFunctionRows(), NumberOfNonzeroControlPoints(),
            DerivativeRow, ControlPointIndex);

        return mShapeFunctionValues[index];
    }

    double& ShapeFunctionValue(
        const SizeType DerivativeRow,
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV)
    {
        const IndexType index = this->GetIndex(DerivativeRow, ControlPointIndexU, ControlPointIndexV);

        return mShapeFunctionValues[index];
    }

    ///@}
    ///@name Member Variables
    ///@{
    int mDerivativeOrder;
    NurbsCurveShapeFunction mShapeFunctionsU;
    NurbsCurveShapeFunction mShapeFunctionsV;
    Vector mWeightedSums;
    Vector mShapeFunctionValues;
    IndexType mFirstNonzeroControlPointU;
    IndexType mFirstNonzeroControlPointV;
    ///@}
    }; // class NurbsSurfaceShapeFunction

} // namespace Kratos

#endif // KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED defined 