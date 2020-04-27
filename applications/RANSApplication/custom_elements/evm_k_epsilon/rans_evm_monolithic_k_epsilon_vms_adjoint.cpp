//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <algorithm>
#include <iterator>

// External includes

// Project includes

// Application includes
#include "rans_evm_monolithic_k_epsilon_vms_adjoint.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::Initialize()
{
    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmMonolithicKEpsilonVMSAdjoint>(
        NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmMonolithicKEpsilonVMSAdjoint>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmMonolithicKEpsilonVMSAdjoint>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
int RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Check(rCurrentProcessInfo);

    AdjointFluidElement adjoint_fluid_element(this->Id(), this->pGetGeometry());
    AdjointKElement adjoint_k_element(this->Id(), this->pGetGeometry());
    AdjointEpsilonElement adjoint_epsilon_element(this->Id(), this->pGetGeometry());

    adjoint_fluid_element.Check(rCurrentProcessInfo);
    adjoint_k_element.Check(rCurrentProcessInfo);
    adjoint_epsilon_element.Check(rCurrentProcessInfo);

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        NodeType& r_node = this->GetGeometry()[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_3, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_SCALAR_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(AUX_ADJOINT_FLUID_VECTOR_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_3, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUX_ADJOINT_SCALAR_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_3, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUX_ADJOINT_SCALAR_2, r_node);
    }

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TElementLocalSize)
        rResult.resize(TElementLocalSize);

    AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
    AdjointKElement k_element(this->Id(), this->pGetGeometry());
    AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

    fluid_element.SetData(this->Data());
    k_element.SetData(this->Data());
    epsilon_element.SetData(this->Data());

    std::array<std::size_t, TFluidLocalSize> fluid_ids;
    fluid_element.EquationIdArray(fluid_ids, rCurrentProcessInfo);
    AssignSubArray(fluid_ids, rResult, VelPresBlock());
    std::array<std::size_t, TKLocalSize> k_ids;
    k_element.EquationIdArray(k_ids, rCurrentProcessInfo);
    AssignSubArray(k_ids, rResult, KBlock());
    std::array<std::size_t, TEpsilonLocalSize> epsilon_ids;
    epsilon_element.EquationIdArray(epsilon_ids, rCurrentProcessInfo);
    AssignSubArray(epsilon_ids, rResult, EpsilonBlock());
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != TElementLocalSize)
        rElementalDofList.resize(TElementLocalSize);

    AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
    AdjointKElement k_element(this->Id(), this->pGetGeometry());
    AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

    fluid_element.SetData(this->Data());
    k_element.SetData(this->Data());
    epsilon_element.SetData(this->Data());

    std::array<Dof<double>::Pointer, TFluidLocalSize> fluid_dofs;
    fluid_element.GetDofArray(fluid_dofs, rCurrentProcessInfo);
    AssignSubArray(fluid_dofs, rElementalDofList, VelPresBlock());
    std::array<Dof<double>::Pointer, TKLocalSize> k_dofs;
    k_element.GetDofArray(k_dofs, rCurrentProcessInfo);
    AssignSubArray(k_dofs, rElementalDofList, KBlock());
    std::array<Dof<double>::Pointer, TEpsilonLocalSize> epsilon_dofs;
    epsilon_element.GetDofArray(epsilon_dofs, rCurrentProcessInfo);
    AssignSubArray(epsilon_dofs, rElementalDofList, EpsilonBlock());
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::GetValuesVector(VectorType& rValues,
                                                                           int Step)
{
    if (rValues.size() != TElementLocalSize)
        rValues.resize(TElementLocalSize, false);

    AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
    AdjointKElement k_element(this->Id(), this->pGetGeometry());
    AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

    fluid_element.SetData(this->Data());
    k_element.SetData(this->Data());
    epsilon_element.SetData(this->Data());

    std::array<double, TFluidLocalSize> fluid_values;
    fluid_element.GetValuesArray(fluid_values, Step);
    AssignSubArray(fluid_values, rValues, VelPresBlock());
    std::array<double, TKLocalSize> k_values;
    k_element.GetValuesArray(k_values, Step);
    AssignSubArray(k_values, rValues, KBlock());
    std::array<double, TEpsilonLocalSize> epsilon_values;
    epsilon_element.GetValuesArray(epsilon_values, Step);
    AssignSubArray(epsilon_values, rValues, EpsilonBlock());
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::GetFirstDerivativesVector(
    VectorType& rValues, int Step)
{
    if (rValues.size() != TElementLocalSize)
        rValues.resize(TElementLocalSize, false);

    AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
    AdjointKElement k_element(this->Id(), this->pGetGeometry());
    AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

    fluid_element.SetData(this->Data());
    k_element.SetData(this->Data());
    epsilon_element.SetData(this->Data());

    std::array<double, TFluidLocalSize> fluid_values;
    fluid_element.GetFirstDerivativesArray(fluid_values, Step);
    AssignSubArray(fluid_values, rValues, VelPresBlock());
    std::array<double, TKLocalSize> k_values;
    k_element.GetFirstDerivativesArray(k_values, Step);
    AssignSubArray(k_values, rValues, KBlock());
    std::array<double, TEpsilonLocalSize> epsilon_values;
    epsilon_element.GetFirstDerivativesArray(epsilon_values, Step);
    AssignSubArray(epsilon_values, rValues, EpsilonBlock());
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::GetSecondDerivativesVector(
    VectorType& rValues, int Step)
{
    if (rValues.size() != TElementLocalSize)
        rValues.resize(TElementLocalSize, false);

    AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
    AdjointKElement k_element(this->Id(), this->pGetGeometry());
    AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

    fluid_element.SetData(this->Data());
    k_element.SetData(this->Data());
    epsilon_element.SetData(this->Data());

    std::array<double, TFluidLocalSize> fluid_values;
    fluid_element.GetSecondDerivativesArray(fluid_values, Step);
    AssignSubArray(fluid_values, rValues, VelPresBlock());
    std::array<double, TKLocalSize> k_values;
    k_element.GetSecondDerivativesArray(k_values, Step);
    AssignSubArray(k_values, rValues, KBlock());
    std::array<double, TEpsilonLocalSize> epsilon_values;
    epsilon_element.GetSecondDerivativesArray(epsilon_values, Step);
    AssignSubArray(epsilon_values, rValues, EpsilonBlock());
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmMonolithicKEpsilonVMSAdjoint::"
                    "CalculateLocalSystem method is not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& /*rCurrentProcessInfo*/)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize)
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& /*rCurrentProcessInfo*/)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmMonolithicKEpsilonVMSAdjoint::"
                    "CalculateRightHandSide method is not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
    AdjointKElement k_element(this->Id(), this->pGetGeometry());
    AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

    fluid_element.SetData(this->Data());
    k_element.SetData(this->Data());
    epsilon_element.SetData(this->Data());

    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize)
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

    rLeftHandSideMatrix.clear();

    BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> vms_vms;
    fluid_element.CalculateFirstDerivativesLHS(vms_vms, rCurrentProcessInfo);
    AssignSubMatrix(vms_vms, rLeftHandSideMatrix, VelPresBlock(), VelPresBlock());

    BoundedMatrix<double, TNumNodes, TFluidLocalSize> vms_k;
    fluid_element.CalculateResidualScalarDerivatives(
        TURBULENT_KINETIC_ENERGY, vms_k, rCurrentProcessInfo);
    AssignSubMatrix(vms_k, rLeftHandSideMatrix, KBlock(), VelPresBlock());

    BoundedMatrix<double, TNumNodes, TFluidLocalSize> vms_epsilon;
    fluid_element.CalculateResidualScalarDerivatives(
        TURBULENT_ENERGY_DISSIPATION_RATE, vms_epsilon, rCurrentProcessInfo);
    AssignSubMatrix(vms_epsilon, rLeftHandSideMatrix, EpsilonBlock(), VelPresBlock());

    BoundedMatrix<double, TCoordLocalSize, TNumNodes> k_vms;
    k_element.CalculateElementTotalResidualVelocityDerivatives(k_vms, rCurrentProcessInfo);
    AssignSubMatrix(k_vms, rLeftHandSideMatrix, VelBlock(), KBlock());

    BoundedMatrix<double, TNumNodes, TNumNodes> k_k;
    k_element.CalculateFirstDerivativesLHS(k_k, rCurrentProcessInfo);
    AssignSubMatrix(k_k, rLeftHandSideMatrix, KBlock(), KBlock());

    BoundedMatrix<double, TNumNodes, TNumNodes> k_epsilon;
    k_element.CalculateResidualScalarDerivatives(
        TURBULENT_ENERGY_DISSIPATION_RATE, k_epsilon, rCurrentProcessInfo);
    AssignSubMatrix(k_epsilon, rLeftHandSideMatrix, EpsilonBlock(), KBlock());

    BoundedMatrix<double, TCoordLocalSize, TNumNodes> epsilon_vms;
    epsilon_element.CalculateElementTotalResidualVelocityDerivatives(
        epsilon_vms, rCurrentProcessInfo);
    AssignSubMatrix(epsilon_vms, rLeftHandSideMatrix, VelBlock(), EpsilonBlock());

    BoundedMatrix<double, TNumNodes, TNumNodes> epsilon_k;
    epsilon_element.CalculateResidualScalarDerivatives(
        TURBULENT_KINETIC_ENERGY, epsilon_k, rCurrentProcessInfo);
    AssignSubMatrix(epsilon_k, rLeftHandSideMatrix, KBlock(), EpsilonBlock());

    BoundedMatrix<double, TNumNodes, TNumNodes> epsilon_epsilon;
    epsilon_element.CalculateFirstDerivativesLHS(epsilon_epsilon, rCurrentProcessInfo);
    AssignSubMatrix(epsilon_epsilon, rLeftHandSideMatrix, EpsilonBlock(),
                    EpsilonBlock());

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
    AdjointKElement k_element(this->Id(), this->pGetGeometry());
    AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

    fluid_element.SetData(this->Data());
    k_element.SetData(this->Data());
    epsilon_element.SetData(this->Data());

    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize)
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

    rLeftHandSideMatrix.clear();

    BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> vms_vms;
    fluid_element.CalculateSecondDerivativesLHS(vms_vms, rCurrentProcessInfo);
    AssignSubMatrix(vms_vms, rLeftHandSideMatrix, VelPresBlock(), VelPresBlock());

    BoundedMatrix<double, TNumNodes, TNumNodes> k_k;
    k_element.CalculateSecondDerivativesLHS(k_k, rCurrentProcessInfo);
    AssignSubMatrix(k_k, rLeftHandSideMatrix, KBlock(), KBlock());

    BoundedMatrix<double, TNumNodes, TNumNodes> epsilon_epsilon;
    epsilon_element.CalculateSecondDerivativesLHS(epsilon_epsilon, rCurrentProcessInfo);
    AssignSubMatrix(epsilon_epsilon, rLeftHandSideMatrix, EpsilonBlock(), EpsilonBlock());

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::CalculateMassMatrix(
    MatrixType& rMassMatrix, ProcessInfo& /*rCurrentProcessInfo*/)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmMonolithicKEpsilonVMSAdjoint::"
                    "CalculateMassMatrix method is not implemented.";

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& /*rCurrentProcessInfo*/)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmMonolithicKEpsilonVMSAdjoint::"
                    "CalculateDampingMatrix method is not implemented.";

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmMonolithicKEpsilonVMSAdjoint<TDim, TNumNodes>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY)
    {
        if (rOutput.size1() != TCoordLocalSize || rOutput.size2() != TElementLocalSize)
            rOutput.resize(TCoordLocalSize, TElementLocalSize, false);

        rOutput.clear();

        AdjointFluidElement fluid_element(this->Id(), this->pGetGeometry());
        AdjointKElement k_element(this->Id(), this->pGetGeometry());
        AdjointEpsilonElement epsilon_element(this->Id(), this->pGetGeometry());

        fluid_element.SetData(this->Data());
        k_element.SetData(this->Data());
        epsilon_element.SetData(this->Data());

        BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize> vms_residuals;
        fluid_element.CalculateSensitivityMatrix(
            rSensitivityVariable, vms_residuals, rCurrentProcessInfo);
        AssignSubMatrix(vms_residuals, rOutput, CoordBlock(), VelPresBlock());

        BoundedMatrix<double, TCoordLocalSize, TNumNodes> k_residuals;
        k_element.CalculateSensitivityMatrix(rSensitivityVariable, k_residuals,
                                             rCurrentProcessInfo);
        AssignSubMatrix(k_residuals, rOutput, CoordBlock(), KBlock());

        BoundedMatrix<double, TCoordLocalSize, TNumNodes> epsilon_residuals;
        epsilon_element.CalculateSensitivityMatrix(
            rSensitivityVariable, epsilon_residuals, rCurrentProcessInfo);
        AssignSubMatrix(epsilon_residuals, rOutput, CoordBlock(), EpsilonBlock());
    }
    else
    {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

// template instantiations
template class RansEvmMonolithicKEpsilonVMSAdjoint<2>;
template class RansEvmMonolithicKEpsilonVMSAdjoint<3>;

///@}
///@}

} // namespace Kratos
