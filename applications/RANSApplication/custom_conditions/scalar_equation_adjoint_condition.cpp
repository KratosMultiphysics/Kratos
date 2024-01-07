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

// Project includes
#include "includes/cfd_variables.h"
#include "includes/variables.h"
#include "utilities/indirect_scalar.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "scalar_equation_adjoint_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(1);
    rVector[0] = MakeIndirectScalar(r_node, RANS_SCALAR_1_ADJOINT_1.GetTimeDerivative(), Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(1);
    rVector[0] = MakeIndirectScalar(r_node, RANS_SCALAR_1_ADJOINT_1.GetTimeDerivative().GetTimeDerivative(), Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(1);
    rVector[0] = MakeIndirectScalar(r_node, RANS_SCALAR_1_ADJOINT_1.GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(), Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &RANS_SCALAR_1_ADJOINT_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &RANS_SCALAR_1_ADJOINT_3;
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &RANS_AUX_ADJOINT_SCALAR_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes);
    }

    for (IndexType i = 0; i < TNumNodes; ++i) {
        rResult[i] = this->GetGeometry()[i].GetDof(RANS_SCALAR_1_ADJOINT_1).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes);
    }

    for (IndexType i = 0; i < TNumNodes; ++i) {
        rResult[i] = this->GetGeometry()[i].pGetDof(RANS_SCALAR_1_ADJOINT_1);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::GetValuesVector(
    Vector& Values,
    int Step) const
{
    if (Values.size() != TNumNodes) {
        Values.resize(TNumNodes);
    }

    for (IndexType i = 0; i < TNumNodes; ++i) {
        Values[i] = this->GetGeometry()[i].FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_1);
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::GetFirstDerivativesVector(
    Vector& Values,
    int Step) const
{
    if (Values.size() != TNumNodes) {
        Values.resize(TNumNodes, false);
    }

    Values.clear();
}

template<unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::GetSecondDerivativesVector(
    Vector& Values,
    int Step) const
{
    if (Values.size() != TNumNodes) {
        Values.resize(TNumNodes, false);
    }

    for (IndexType i = 0; i < TNumNodes; ++i) {
        Values[i] = this->GetGeometry()[i].FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int ScalarEquationAdjointCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize and initialize output
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    rRightHandSideVector.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    Matrix& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    rLeftHandSideMatrix.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::CalculateFirstDerivativesLHS(
    Matrix& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    rLeftHandSideMatrix.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::CalculateSecondDerivativesLHS(
    Matrix& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    rLeftHandSideMatrix.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ScalarEquationAdjointCondition<TDim, TNumNodes>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rDesignVariable == SHAPE_SENSITIVITY) {
        if (rOutput.size1() != TCoordsLocalSize || rOutput.size2() != TNumNodes) {
            rOutput.resize(TCoordsLocalSize, TNumNodes, false);
        }

        rOutput.clear();
    } else {
        KRATOS_ERROR << "Partial sensitivity w.r.t. " << rDesignVariable.Name()
                     << " not supported.";
    }

    KRATOS_CATCH("");
}

template class ScalarEquationAdjointCondition<2, 2>;

} // namespace Kratos
