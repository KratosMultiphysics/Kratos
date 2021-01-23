//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_NURBS_VOLUME_SHAPE_FUNCTIONS_H_INCLUDED )
#define  KRATOS_NURBS_VOLUME_SHAPE_FUNCTIONS_H_INCLUDED

// Project includes
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"
#include "input_output/logger.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/**
 * @class NurbsVolumeShapeFunctions
 * @ingroup KratosCore
 * @brief Provides the shape functions and evaluation methods for trivariant BSpline volumes.
 * @details The implementation corresponds to a 3D-BSpline, as weights are not yet implemented.
 **/
class NurbsVolumeShapeFunction
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NurbsVolumeShapeFunction()
    {
    };

    // Constructor using the degree u, degree v, degree w and the order of shape functions.
    // This is required to make an optimized memory management.
    NurbsVolumeShapeFunction(
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const SizeType PolynomialDegreeW,
        const SizeType DerivativeOrder)
    {
        ResizeDataContainers(PolynomialDegreeU, PolynomialDegreeV, PolynomialDegreeW, DerivativeOrder);
    }

    ///@}
    ///@name Static Operations
    ///@{
    /**
     * @brief   Returns the number of shape function rows for a given order.
     * @details The shape functions are provided in the following order:
     *          N | dN/du, dN/dv, dN/dw | dN^2/du^2, dN^2/du*dv, dN^2/du*dw, dN^2/v^2, dN^2/dv*dw, dN^2/dw^2 | ...
     *          0 | (1,0,0), (0,1,0), (0,0,1) | (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2) | ...
     * @param   DerivativesOrder 0: the shape functions N only, 1: N and first derivatives, ...
     *          Thus, for DerivativesOrder 0 the return value is 1, for 1 the return value is 4, ...
     * @return  NumberOfShapeFunctionRows.
     */
    static inline SizeType NumberOfShapeFunctionRows(const SizeType DerivativeOrder)
    {
        IndexType number_of_shape_function_rows = 0.0;
        for( IndexType i = 0; i < DerivativeOrder + 1; ++i){
            number_of_shape_function_rows += (1 + i) * (2 + i) / 2;
        }

        return number_of_shape_function_rows;
    }
    /**
     * @brief   Returns the index of the shape function row for the given derivative indices.
     * @details The shape functions are provided in the following order:
     *          N | dN/du, dN/dv, dN/dw | dN^2/du^2, dN^2/du*dv, dN^2/du*dw, dN^2/v^2, dN^2/dv*dw, dN^2/dw^2 | ...
     *          0 | (1,0,0), (0,1,0), (0,0,1) | (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2) | ...
     * @param   DerivativesOrderU Derivative order along u-direction.
     * @param   DerivativesOrderV Derivative order along v-direction.
     * @param   DerivativesOrderW Derivative order along w-direction.
     * @return  IndexOfShapeFunctionRow
     **/
    static inline IndexType IndexOfShapeFunctionRow(
        const SizeType DerivativeOrderU,
        const SizeType DerivativeOrderV,
        const SizeType DerivativeOrderW)
    {
        // This seraializes the pascals pyramid.
        const SizeType current_level = DerivativeOrderU + DerivativeOrderV + DerivativeOrderW;

        SizeType first_index_of_current_level = 0;
        if( current_level > 0){
            first_index_of_current_level = NumberOfShapeFunctionRows(current_level - 1);
            const SizeType offset = current_level - DerivativeOrderU;
            SizeType index_in_current_level = 0;
            for( IndexType i = 0; i < offset; ++i){
                index_in_current_level += i + 1;
            }
            index_in_current_level += DerivativeOrderW;

            return first_index_of_current_level + index_in_current_level;
        }
        else {
            return 0;
        }
    }

    ///@}
    ///@name Operators
    ///@{

    double operator()(
        const IndexType ControlPointIndex,
        const IndexType DerivativeRow) const
    {
        return ShapeFunctionValue(ControlPointIndex, DerivativeRow);
    }

    double operator()(
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV,
        const IndexType ControlPointIndexW,
        const IndexType DerivativeRow) const
    {
        return ShapeFunctionValue(ControlPointIndexU, ControlPointIndexV, ControlPointIndexW, DerivativeRow);
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Prepares the required shape function data containers.
     **/
    void ResizeDataContainers(
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const SizeType PolynomialDegreeW,
        const SizeType DerivativeOrder)
    {
        const SizeType number_of_shape_function_rows = this->NumberOfShapeFunctionRows(DerivativeOrder);
        const SizeType number_of_nonzero_control_points = (PolynomialDegreeU + 1) * (PolynomialDegreeV + 1) * (PolynomialDegreeW + 1);

        mShapeFunctionsU.ResizeDataContainers(PolynomialDegreeU, DerivativeOrder);
        mShapeFunctionsV.ResizeDataContainers(PolynomialDegreeV, DerivativeOrder);
        mShapeFunctionsW.ResizeDataContainers(PolynomialDegreeW, DerivativeOrder);
        mShapeFunctionValues.resize(number_of_nonzero_control_points * number_of_shape_function_rows);
        // mWeightedSums.resize(number_of_shape_function_rows);

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

    SizeType PolynomialDegreeW() const
    {
        return mShapeFunctionsW.PolynomialDegree();
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

    SizeType NumberOfNonzeroControlPointsW() const
    {
        return mShapeFunctionsW.NumberOfNonzeroControlPoints();
    }

    SizeType NumberOfNonzeroControlPoints() const
    {
        return NumberOfNonzeroControlPointsU() * NumberOfNonzeroControlPointsV() * NumberOfNonzeroControlPointsW();
    }

    /**
     * @brief Computes indices of the nonzero control points.
     * @return std::vector<array_1d<int, 3>> which holds the indices in u,v, and w-direction.
     **/
    std::vector<array_1d<int, 3>> NonzeroControlPointIndices() const
    {
        std::vector<array_1d<int, 3>> indices(NumberOfNonzeroControlPoints());

        for (IndexType i = 0; i < NumberOfNonzeroControlPointsU(); ++i) {
            for (IndexType j = 0; j < NumberOfNonzeroControlPointsV(); ++j) {
                for (IndexType k = 0; k < NumberOfNonzeroControlPointsW(); ++k) {
                    IndexType poleIndex = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        NumberOfNonzeroControlPointsU(), NumberOfNonzeroControlPointsV(), NumberOfNonzeroControlPointsW(), i, j, k);

                    IndexType poleU = GetFirstNonzeroControlPointU() + i;
                    IndexType poleV = GetFirstNonzeroControlPointV() + j;
                    IndexType poleW = GetFirstNonzeroControlPointW() + k;

                    indices[poleIndex][0] = poleU;
                    indices[poleIndex][1] = poleV;
                    indices[poleIndex][2] = poleW;
                }
            }
        }
        return indices;
    }

    /**
     * @return Indices of control points for given discretization.
     * @param Number control points in u-direction.
     * @param Number control points in v-direction.
     * @param Number control points in w-direction.
     * @return Vector with indices.
     **/
    std::vector<int> ControlPointIndices(
        SizeType NumberOfControlPointsU, SizeType NumberOfControlPointsV, SizeType NumberOfControlPointsW) const
    {
        std::vector<int> indices(NumberOfNonzeroControlPoints());

        for (IndexType i = 0; i < NumberOfNonzeroControlPointsU(); ++i) {
            for (IndexType j = 0; j < NumberOfNonzeroControlPointsV(); ++j) {
                for(IndexType k = 0; k < NumberOfNonzeroControlPointsW(); ++k) {
                    IndexType poleIndex = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        NumberOfNonzeroControlPointsU(), NumberOfNonzeroControlPointsV(), NumberOfNonzeroControlPointsW(), i, j, k);

                    IndexType cp_index_u = GetFirstNonzeroControlPointU() + i;
                    IndexType cp_index_v = GetFirstNonzeroControlPointV() + j;
                    IndexType cp_index_w = GetFirstNonzeroControlPointW() + k;

                    indices[poleIndex] = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        NumberOfControlPointsU, NumberOfControlPointsV, NumberOfControlPointsW, cp_index_u, cp_index_v, cp_index_w);
                }
            }
        }

        return indices;
    }

    /**
     * @brief Get shape function value for given CP-indices and derivative row.
     **/
    double ShapeFunctionValue(
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV,
        const IndexType ControlPointIndexW,
        const SizeType DerivativeRow) const
    {
        const IndexType index = this->GetIndex(
            ControlPointIndexU,
            ControlPointIndexV,
            ControlPointIndexW,
            DerivativeRow);

        KRATOS_DEBUG_ERROR_IF(index >= mShapeFunctionValues.size()) << "Index exceeds size of shape function values. Index: "
            << index << ", mShapeFunctionValues.size(): " << mShapeFunctionValues.size()
            << ", ControlPointIndexU: " << ControlPointIndexU << ", ControlPointIndexV: " << ControlPointIndexV
            << ", ControlPointIndexW: " << ControlPointIndexW << ", DerivativeRow: " << DerivativeRow << std::endl;

        return mShapeFunctionValues[index];
    }

    /**
     * @brief Get shape function value for given CP-index and derivative row.
     **/
    double ShapeFunctionValue(
        const IndexType ControlPointIndex,
        const SizeType DerivativeOrder) const
    {
        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NumberOfShapeFunctionRows(), NumberOfNonzeroControlPoints(),
            DerivativeOrder, ControlPointIndex);

        KRATOS_DEBUG_ERROR_IF(index >= mShapeFunctionValues.size()) << "Index exceeds size of shape function values. Index: "
            << index << ", mShapeFunctionValues.size(): " << mShapeFunctionValues.size()
            << ", ControlPointIndex: " << ControlPointIndex << ", DerivativeOrder: " << DerivativeOrder << std::endl;

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

    IndexType GetFirstNonzeroControlPointW() const
    {
        return mFirstNonzeroControlPointW;
    }

    IndexType GetLastNonzeroControlPointW() const
    {
        return GetFirstNonzeroControlPointW() + PolynomialDegreeW();
    }

    /**
     * @brief Computes the shape function values at the given parameter and span.
     **/
    void ComputeBSplineShapeFunctionValuesAtSpan(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rKnotsW,
        const int SpanU,
        const int SpanV,
        const int SpanW,
        const double ParameterU,
        const double ParameterV,
        const double ParameterW)
    {
        mShapeFunctionValues = ZeroVector(mShapeFunctionValues.size());

        mFirstNonzeroControlPointU = SpanU - PolynomialDegreeU() + 1;
        mFirstNonzeroControlPointV = SpanV - PolynomialDegreeV() + 1;
        mFirstNonzeroControlPointW = SpanW - PolynomialDegreeW() + 1;

        // Compute 1D shape functions
        mShapeFunctionsU.ComputeBSplineShapeFunctionValuesAtSpan(rKnotsU, SpanU, ParameterU);
        mShapeFunctionsV.ComputeBSplineShapeFunctionValuesAtSpan(rKnotsV, SpanV, ParameterV);
        mShapeFunctionsW.ComputeBSplineShapeFunctionValuesAtSpan(rKnotsW, SpanW, ParameterW);
        // Compute 3D shape functions
        for( IndexType currrent_order = 0; currrent_order <= DerivativeOrder(); ++currrent_order){
            for (int OrderU = currrent_order; OrderU >= 0; OrderU = OrderU - 1) {
                const IndexType difference = currrent_order - OrderU;
                for (IndexType OrderW = 0; OrderW <= difference; ++OrderW) {
                    const IndexType OrderV = difference - OrderW;
                    for (IndexType a = 0; a < NumberOfNonzeroControlPointsU(); ++a) {
                        for (IndexType b = 0; b < NumberOfNonzeroControlPointsV(); ++b) {
                            for (IndexType c = 0; c < NumberOfNonzeroControlPointsW(); ++c) {
                                const IndexType index = IndexOfShapeFunctionRow(OrderU, OrderV, OrderW);
                                ShapeFunctionValue(a, b, c, index) = mShapeFunctionsU(a, OrderU) * mShapeFunctionsV(b, OrderV) * mShapeFunctionsW(c, OrderW);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Computes the shape function values at the given parameter.
     * @details Finds the span, where the given parameters are located and
     *          computes then the shape function values.
     **/
    void ComputeBSplineShapeFunctionValues(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rKnotsW,
        const double ParameterU,
        const double ParameterV,
        const double ParameterW)
    {
        const IndexType SpanU = NurbsUtilities::GetLowerSpan(PolynomialDegreeU(), rKnotsU, ParameterU);
        const IndexType SpanV = NurbsUtilities::GetLowerSpan(PolynomialDegreeV(), rKnotsV, ParameterV);
        const IndexType SpanW = NurbsUtilities::GetLowerSpan(PolynomialDegreeW(), rKnotsW, ParameterW);

        ComputeBSplineShapeFunctionValuesAtSpan(
            rKnotsU,
            rKnotsV,
            rKnotsW,
            SpanU,
            SpanV,
            SpanW,
            ParameterU,
            ParameterV,
            ParameterW);
    }

    void ComputeNurbsShapeFunctionValuesAtSpan(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rKnotsW,
        const IndexType SpanU,
        const IndexType SpanV,
        const IndexType SpanW,
        const Vector& Weights,
        const double ParameterU,
        const double ParameterV,
        const double ParameterW)
    {
        // Check input
        KRATOS_WARNING("NURBSVolumeShapeFunctions: ComputeNurbsShapeFunctionValuesAtSpan")
            << "The NURBS weights are not yet imlemented. "
            << "Therefore the 'ComputeBSplineShapeFunctionValuesAtSpan' is invoked." << std::endl;

        // Compute B-Spline shape functions
        ComputeBSplineShapeFunctionValuesAtSpan(
            rKnotsU, rKnotsV, rKnotsW, SpanU, SpanV, SpanW, ParameterU, ParameterV, ParameterW);
    }

    void ComputeNurbsShapeFunctionValues(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rKnotsW,
        const Vector& Weights,
        const double ParameterU,
        const double ParameterV,
        const double ParameterW)
    {
        const IndexType SpanU = NurbsUtilities::GetLowerSpan(PolynomialDegreeU(), rKnotsU, ParameterU);
        const IndexType SpanV = NurbsUtilities::GetLowerSpan(PolynomialDegreeV(), rKnotsV, ParameterV);
        const IndexType SpanW = NurbsUtilities::GetLowerSpan(PolynomialDegreeW(), rKnotsW, ParameterW);

        ComputeNurbsShapeFunctionValuesAtSpan(
            rKnotsU, rKnotsV, rKnotsW, SpanU, SpanV, SpanW, Weights, ParameterU, ParameterV, ParameterW);
    }
    ///@}
private:
    ///@name Operations
    ///@{

    // Attention: Weight ar not yet implemented.    
    // double& GetWeightedSum(const IndexType Index)
    // {
    //     return mWeightedSums[Index];
    // }

    // double& GetWeightedSum(const SizeType DerivativeOrderU, const SizeType DerivativeOrderV, const SizeType DerivativeOrderW)
    // {
    //     const int index = IndexOfShapeFunctionRow(DerivativeOrderU, DerivativeOrderV, DerivativeOrderW);

    //     return GetWeightedSum(index);
    // }

    inline int GetControlPointIndex(
        const SizeType NumberOfKnotsU,
        const SizeType NumberOfKnotsV,
        const SizeType NumberOfKnotsW,
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV,
        const IndexType ControlPointIndexW)
        const
    {
        return NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NurbsUtilities::GetNumberOfControlPoints(PolynomialDegreeU(), NumberOfKnotsU),
            NurbsUtilities::GetNumberOfControlPoints(PolynomialDegreeV(), NumberOfKnotsV),
            NurbsUtilities::GetNumberOfControlPoints(PolynomialDegreeW(), NumberOfKnotsW),
            ControlPointIndexU, ControlPointIndexV, ControlPointIndexW);
    }

    inline int GetNonzeroControlPointIndex(
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV,
        const IndexType ControlPointIndexW)
        const
    {
        return NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NumberOfNonzeroControlPointsU(), NumberOfNonzeroControlPointsV(), NumberOfNonzeroControlPointsW(),
            ControlPointIndexU, ControlPointIndexV, ControlPointIndexW);
    }

    inline int GetIndex(
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV,
        const IndexType ControlPointIndexW,
        const SizeType DerivativeRow)
        const
    {
        const IndexType control_point_index = GetNonzeroControlPointIndex(
            ControlPointIndexU, ControlPointIndexV, ControlPointIndexW);

        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NumberOfShapeFunctionRows(), NumberOfNonzeroControlPoints(),
            DerivativeRow, control_point_index);

        return index;
    }

    double& ShapeFunctionValue(
        const IndexType ControlPointIndex,
        const IndexType DerivativeRow)
    {
        const IndexType index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
            NumberOfShapeFunctionRows(), NumberOfNonzeroControlPoints(),
            DerivativeRow, ControlPointIndex);

        return mShapeFunctionValues[index];
    }

    double& ShapeFunctionValue(
        const IndexType ControlPointIndexU,
        const IndexType ControlPointIndexV,
        const IndexType ControlPointIndexW,
        const SizeType DerivativeRow)
    {
        const IndexType index = this->GetIndex(ControlPointIndexU, ControlPointIndexV, ControlPointIndexW, DerivativeRow);

        return mShapeFunctionValues[index];
    }

    ///@}
    ///@name Member Variables
    ///@{
        
    int mDerivativeOrder;
    NurbsCurveShapeFunction mShapeFunctionsU;
    NurbsCurveShapeFunction mShapeFunctionsV;
    NurbsCurveShapeFunction mShapeFunctionsW;
    Vector mShapeFunctionValues;
    IndexType mFirstNonzeroControlPointU;
    IndexType mFirstNonzeroControlPointV;
    IndexType mFirstNonzeroControlPointW;
    ///@}
    }; // class NurbsVolumeShapeFunction

} // namespace Kratos

#endif // KRATOS_NURBS_VOLUME_SHAPE_FUNCTIONS_H_INCLUDED defined