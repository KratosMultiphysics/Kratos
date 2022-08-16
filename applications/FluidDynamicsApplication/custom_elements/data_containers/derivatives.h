//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_DERIVATIVES_H)
#define KRATOS_DERIVATIVES_H

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

namespace DerivativeHelperUtilities
{

template<unsigned int TNumNodes, unsigned int TColumnOffset, unsigned int TSubVectorBlockSize, unsigned int TAssemblyBlockSize>
static void inline AssembleSubVectorToMatrix(
    Matrix& rOutput,
    const IndexType RowIndex,
    const BoundedVector<double, TSubVectorBlockSize * TNumNodes>& rSubVector)
{
    KRATOS_TRY

    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TSubVectorBlockSize; ++j) {
            rOutput(RowIndex, i * TAssemblyBlockSize + j + TColumnOffset) += rSubVector[i * TSubVectorBlockSize + j];
        }
    }

    KRATOS_CATCH("");
}

}

class ZeroDerivatives
{
public:
    ///@name Local type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Public static operations
    ///@{

    template<IndexType TSize, class... TArgs>
    static void inline CalculateGaussPointResidualsDerivativeContributions(BoundedVector<double, TSize>& rResidualDerivative, TArgs&... rArgs)
    {
        rResidualDerivative.clear();
    }

    ///@}
};

template<unsigned int TRowStartingIndex, unsigned int TColumnStartingIndex, class TSubAssemblyType>
class SubAssembly
{
public:
    ///@name Type definitions
    ///@{

    static constexpr unsigned int RowStartingIndex = TRowStartingIndex;

    static constexpr unsigned int ColumnStartingIndex = TColumnStartingIndex;

    using SubAssemblyType = TSubAssemblyType;

    ///@}
};

template <
    class TResidualData,
    class... TResidualDerivatives
>
class FirstDerivatives
{
public:
    ///@name Local type definitions
    ///@{

    using IndexType = std::size_t;

    static constexpr IndexType TBlockSize = TResidualData::TBlockSize;

    static constexpr IndexType TNumNodes = TResidualData::TLNumNodes;

    static constexpr IndexType TResidualsSize = TNumNodes * TBlockSize;

    using ResidualVector = BoundedVector<double, TResidualsSize>;

    ///@}
    ///@name Derivative type definitions
    ///@{

    using ElementDataType = TResidualData;

    ///@}
    ///@name Public static operations
    ///@{

    template<unsigned int TAssemblyBlockSize>
    static void inline CalculateAndAddGaussPointNodalVariableDerivativeContributions(
        Matrix& rOutput,
        TResidualData& rData,
        ResidualVector& rResidualDerivative,
        const int NodeIndex,
        const double W,
        const Vector& rN,
        const Matrix& rdNdX,
        const double WDerivative,
        const double DetJDerivative,
        const Matrix& rdNdXDerivative,
        const double MassTermsDerivativesWeight)
    {
        (CalculateAndAddGaussPointNodalVariableDerivativeContributions<TAssemblyBlockSize, TResidualDerivatives>(
            rOutput,
            rData,
            rResidualDerivative,
            NodeIndex,
            W,
            rN,
            rdNdX,
            WDerivative,
            DetJDerivative,
            rdNdXDerivative,
            MassTermsDerivativesWeight
        ), ...);
    }

    ///@}

private:
    ///@name Private static operations
    ///@{

    template<unsigned int TAssemblyBlockSize, class TDerivative>
    static void inline CalculateAndAddGaussPointNodalVariableDerivativeContributions(
        Matrix& rOutput,
        TResidualData& rData,
        ResidualVector& rResidualDerivative,
        const int NodeIndex,
        const double W,
        const Vector& rN,
        const Matrix& rdNdX,
        const double WDerivative,
        const double DetJDerivative,
        const Matrix& rdNdXDerivative,
        const double MassTermsDerivativesWeight)
    {
        TDerivative::SubAssemblyType::CalculateGaussPointResidualsDerivativeContributions(
            rResidualDerivative, rData, NodeIndex, W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative, MassTermsDerivativesWeight);

        DerivativeHelperUtilities::AssembleSubVectorToMatrix<TNumNodes, TDerivative::ColumnStartingIndex, TBlockSize,  TAssemblyBlockSize>(
            rOutput, TAssemblyBlockSize * NodeIndex + TDerivative::RowStartingIndex, rResidualDerivative);
    }

    ///@}
};

template <
    class TResidualData,
    class... TResidualDerivatives
>
class SecondDerivatives
{
public:
    ///@name Local type definitions
    ///@{

    using IndexType = std::size_t;

    static constexpr IndexType TBlockSize = TResidualData::TBlockSize;

    static constexpr IndexType TNumNodes = TResidualData::TLNumNodes;

    static constexpr IndexType TResidualsSize = TNumNodes * TBlockSize;

    using ResidualVector = BoundedVector<double, TResidualsSize>;

    ///@}
    ///@name Derivative type definitions
    ///@{

    using ElementDataType = TResidualData;

    ///@}
    ///@name Public static operations
    ///@{

    template<unsigned int TAssemblyBlockSize>
    static void inline CalculateAndAddGaussPointNodalVariableDerivativeContributions(
        Matrix& rOutput,
        TResidualData& rData,
        ResidualVector& rResidualDerivative,
        const int NodeIndex,
        const double W,
        const Vector& rN,
        const Matrix& rdNdX)
    {
        (CalculateAndAddGaussPointNodalVariableDerivativeContributions<TAssemblyBlockSize, TResidualDerivatives>(
            rOutput,
            rData,
            rResidualDerivative,
            NodeIndex,
            W,
            rN,
            rdNdX
        ), ...);
    }

    ///@}

private:
    ///@name Private static operations
    ///@{

    template<unsigned int TAssemblyBlockSize, class TDerivative>
    static void inline CalculateAndAddGaussPointNodalVariableDerivativeContributions(
        Matrix& rOutput,
        TResidualData& rData,
        ResidualVector& rResidualDerivative,
        const int NodeIndex,
        const double W,
        const Vector& rN,
        const Matrix& rdNdX)
    {
        TDerivative::SubAssemblyType::CalculateGaussPointResidualsDerivativeContributions(
            rResidualDerivative, rData, NodeIndex, W, rN, rdNdX);

        DerivativeHelperUtilities::AssembleSubVectorToMatrix<TNumNodes, TDerivative::ColumnStartingIndex, TBlockSize,  TAssemblyBlockSize>(
            rOutput, TAssemblyBlockSize * NodeIndex + TDerivative::RowStartingIndex, rResidualDerivative);
    }

    ///@}
};

}

#endif // KRATOS_DERIVATIVES_H