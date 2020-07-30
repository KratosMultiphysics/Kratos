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

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/serializer.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/line_sensitivity_utility.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "includes/cfd_variables.h"
#include "rans_evm_k_adjoint_wall_condition.h"
#include "rans_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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


template <unsigned int TNumNodes, unsigned int TDim>
RansEvmKAdjointWallCondition<TNumNodes, TDim>& RansEvmKAdjointWallCondition<TNumNodes, TDim>::operator=(
    RansEvmKAdjointWallCondition<TNumNodes, TDim> const& rOther)
{
    Condition::operator=(rOther);

    return *this;
}

template <unsigned int TNumNodes, unsigned int TDim>
Condition::Pointer RansEvmKAdjointWallCondition<TNumNodes, TDim>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmKAdjointWallCondition>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TNumNodes, unsigned int TDim>
Condition::Pointer RansEvmKAdjointWallCondition<TNumNodes, TDim>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmKAdjointWallCondition>(NewId, pGeom, pProperties);
}

template <unsigned int TNumNodes, unsigned int TDim>
Condition::Pointer RansEvmKAdjointWallCondition<TNumNodes, TDim>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmKAdjointWallCondition::"
                    "Clone method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmKAdjointWallCondition::"
                    "CalculateLocalSystem method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

    rLeftHandSideMatrix.clear();
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateLeftHandSide(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.clear();
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmKAdjointWallCondition::"
                    "CalculateRightHandSide method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmKAdjointWallCondition::"
                    "CalculateDampingMatrix method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
int RansEvmKAdjointWallCondition<TNumNodes, TDim>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

    return Check;

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    std::array<std::size_t, TNumNodes> ids;
    this->EquationIdArray(ids, rCurrentProcessInfo);

    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);
    std::copy(ids.begin(), ids.end(), rResult.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::EquationIdArray(
    std::array<std::size_t, TNumNodes>& rResult, ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        rResult[i] = this->GetGeometry()[i].GetDof(RANS_SCALAR_1_ADJOINT_1).EquationId();
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::GetDofList(DofsVectorType& rConditionDofList,
                                                               ProcessInfo& rCurrentProcessInfo)
{
    std::array<Dof<double>::Pointer, TNumNodes> dofs;
    this->GetDofArray(dofs, rCurrentProcessInfo);

    if (rConditionDofList.size() != TNumNodes)
        rConditionDofList.resize(TNumNodes);

    std::copy(dofs.begin(), dofs.end(), rConditionDofList.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::GetDofArray(
    std::array<Dof<double>::Pointer, TNumNodes>& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        rConditionDofList[i] = this->GetGeometry()[i].pGetDof(RANS_SCALAR_1_ADJOINT_1);
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::GetValuesVector(VectorType& rValues, int Step)
{
    std::array<double, TNumNodes> values;
    this->GetValuesArray(values, Step);

    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::GetValuesArray(
    std::array<double, TNumNodes>& rValues, int Step)
{
    const GeometryType& r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rValues[i] = r_geometry[i].FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_1, Step);
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    std::array<double, TNumNodes> values;
    this->GetFirstDerivativesArray(values, Step);
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::GetFirstDerivativesArray(
    std::array<double, TNumNodes>& rValues, int Step)
{
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rValues[i] = 0.0;
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    std::array<double, TNumNodes> values;
    this->GetSecondDerivativesArray(values, Step);
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::GetSecondDerivativesArray(
    std::array<double, TNumNodes>& rValues, int Step)
{
    const GeometryType& r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rValues[i] = r_geometry[i].FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3, Step);
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes ||
        rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

    rLeftHandSideMatrix.clear();
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateFirstDerivativesLHS(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.clear();
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::Calculate(const Variable<Matrix>& rVariable,
                                                              Matrix& rOutput,
                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
    {
        if (rOutput.size1() != TNumNodes || rOutput.size2() != TNumNodes)
            rOutput.resize(TNumNodes, TNumNodes, false);

        rOutput.clear();
    }
    else if (rVariable == RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
    {
        if (rOutput.size1() != TNumNodes || rOutput.size2() != TNumNodes)
            rOutput.resize(TNumNodes, TNumNodes, false);

        rOutput.clear();
    }
    else if (rVariable == RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE)
    {
        constexpr unsigned int local_velocity_size = TNumNodes * TDim;
        if (rOutput.size1() != local_velocity_size || rOutput.size2() != TNumNodes)
            rOutput.resize(local_velocity_size, TNumNodes);
        rOutput.clear();
    }
    else
    {
        KRATOS_ERROR << "Unsupported variable "
                     << rVariable.Name() << " requested at RansEvmKAdjointWallCondition::Calculate.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    rLeftHandSideMatrix.clear();
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY)
    {
        constexpr unsigned int local_coords_size = TNumNodes * TDim;
        if (rOutput.size1() != local_coords_size || rOutput.size2() != TNumNodes)
            rOutput.resize(local_coords_size, TNumNodes, false);

        rOutput.clear();
    }
    else
    {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmKAdjointWallCondition<TNumNodes, TDim>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY)
    {
        rOutput.clear();
    }
    else
    {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

///@}

///@} addtogroup block

// template instantiations
template class RansEvmKAdjointWallCondition<2>;
template class RansEvmKAdjointWallCondition<3>;
} // namespace Kratos.
