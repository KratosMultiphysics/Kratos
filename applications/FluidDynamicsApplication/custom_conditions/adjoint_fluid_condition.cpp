//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes

// Application includes
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "adjoint_fluid_condition.h"

namespace Kratos
{
template <unsigned int TDim>
class AdjointFluidConditionUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    using IndexType = std::size_t;

    using EquationIdVectorType = std::vector<std::size_t>;

    ///@}
    ///@name Static Operations
    ///@{

    static void SetEquationIds(EquationIdVectorType& rResult,
                               IndexType& rLocalIndex,
                               const NodeType& rNode);

    static void SetAdjointValues(Vector& rResult,
                                 IndexType& rLocalIndex,
                                 const NodeType& rNode,
                                 const int Step);

    ///@}
};

template <>
void AdjointFluidConditionUtilities<2>::SetEquationIds(EquationIdVectorType& rResult,
                                                       IndexType& rLocalIndex,
                                                       const NodeType& rNode)
{
    rResult[rLocalIndex++] = rNode.GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId();
    rResult[rLocalIndex++] = rNode.GetDof(ADJOINT_FLUID_VECTOR_1_Y).EquationId();
    rResult[rLocalIndex++] = rNode.GetDof(ADJOINT_FLUID_SCALAR_1).EquationId();
}

template <>
void AdjointFluidConditionUtilities<3>::SetEquationIds(EquationIdVectorType& rResult,
                                                       IndexType& rLocalIndex,
                                                       const NodeType& rNode)
{
    rResult[rLocalIndex++] = rNode.GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId();
    rResult[rLocalIndex++] = rNode.GetDof(ADJOINT_FLUID_VECTOR_1_Y).EquationId();
    rResult[rLocalIndex++] = rNode.GetDof(ADJOINT_FLUID_VECTOR_1_Z).EquationId();
    rResult[rLocalIndex++] = rNode.GetDof(ADJOINT_FLUID_SCALAR_1).EquationId();
}

template <>
void AdjointFluidConditionUtilities<2>::SetAdjointValues(Vector& rResult,
                                                         IndexType& rLocalIndex,
                                                         const NodeType& rNode,
                                                         const int Step)
{
    rResult[rLocalIndex++] = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X, Step);
    rResult[rLocalIndex++] = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_Y, Step);
    rResult[rLocalIndex++] = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
}

template <>
void AdjointFluidConditionUtilities<3>::SetAdjointValues(Vector& rResult,
                                                         IndexType& rLocalIndex,
                                                         const NodeType& rNode,
                                                         const int Step)
{
    rResult[rLocalIndex++] = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X, Step);
    rResult[rLocalIndex++] = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_Y, Step);
    rResult[rLocalIndex++] = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_Z, Step);
    rResult[rLocalIndex++] = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointFluidCondition<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                              const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        AdjointFluidConditionUtilities<TDim>::SetEquationIds(
            rResult, local_index, this->GetGeometry()[i]);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointFluidCondition<TDim, TNumNodes>::CalculateFirstDerivativesLHS(
    Matrix& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != LocalSize ||
        rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    rLeftHandSideMatrix.clear();
}


template <unsigned int TDim, unsigned int TNumNodes>
void AdjointFluidCondition<TDim, TNumNodes>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rOutput.size1() != 0)
        rOutput.resize((TDim) * TNumNodes, LocalSize, false);

    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointFluidCondition<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rResult,
                                                                       int Step) const
{
    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        AdjointFluidConditionUtilities<TDim>::SetAdjointValues(
            rResult, local_index, this->GetGeometry()[i], Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointFluidCondition<TDim, TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    this->GetFirstDerivativesVector(rValues, Step);
}

// template instantiations
template class AdjointFluidCondition<2, 2>;
template class AdjointFluidCondition<3, 3>;
template class AdjointFluidCondition<3, 4>;

} // namespace Kratos.
