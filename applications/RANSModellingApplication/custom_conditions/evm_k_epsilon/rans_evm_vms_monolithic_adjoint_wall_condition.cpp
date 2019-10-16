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
#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "rans_evm_vms_monolithic_adjoint_wall_condition.h"
#include "rans_modelling_application_variables.h"

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

template <unsigned int TDim>
RansEvmVmsMonolithicAdjointWallCondition<TDim>::RansEvmVmsMonolithicAdjointWallCondition(IndexType NewId)
    : Condition(NewId)
{
}

template <unsigned int TDim>
RansEvmVmsMonolithicAdjointWallCondition<TDim>::RansEvmVmsMonolithicAdjointWallCondition(
    IndexType NewId, const NodesArrayType& ThisNodes)
    : Condition(NewId, ThisNodes)
{
}

template <unsigned int TDim>
RansEvmVmsMonolithicAdjointWallCondition<TDim>::RansEvmVmsMonolithicAdjointWallCondition(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

template <unsigned int TDim>
RansEvmVmsMonolithicAdjointWallCondition<TDim>::RansEvmVmsMonolithicAdjointWallCondition(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim>
RansEvmVmsMonolithicAdjointWallCondition<TDim>::RansEvmVmsMonolithicAdjointWallCondition(
    RansEvmVmsMonolithicAdjointWallCondition const& rOther)
    : Condition(rOther)
{
}

template <unsigned int TDim>
RansEvmVmsMonolithicAdjointWallCondition<TDim>::~RansEvmVmsMonolithicAdjointWallCondition()
{
}

template <unsigned int TDim>
RansEvmVmsMonolithicAdjointWallCondition<TDim>& RansEvmVmsMonolithicAdjointWallCondition<TDim>::operator=(
    RansEvmVmsMonolithicAdjointWallCondition<TDim> const& rOther)
{
    Condition::operator=(rOther);
    return *this;
}

template <unsigned int TDim>
Condition::Pointer RansEvmVmsMonolithicAdjointWallCondition<TDim>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmVmsMonolithicAdjointWallCondition>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TDim>
Condition::Pointer RansEvmVmsMonolithicAdjointWallCondition<TDim>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmVmsMonolithicAdjointWallCondition>(
        NewId, pGeom, pProperties);
}

template <unsigned int TDim>
Condition::Pointer RansEvmVmsMonolithicAdjointWallCondition<TDim>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition =
        Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmVmsMonolithicAdjointWallCondition::"
                    "CalculateLocalSystem method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateLeftHandSide(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmVmsMonolithicAdjointWallCondition::"
                    "CalculateRightHandSide method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmVmsMonolithicAdjointWallCondition::"
                    "CalculateDampingMatrix method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
int RansEvmVmsMonolithicAdjointWallCondition<TDim>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

    return Check;

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    EquationIdArrayType ids;
    this->EquationIdArray(ids, rCurrentProcessInfo);

    if (rResult.size() != ids.size())
        rResult.resize(ids.size(), false);
    std::copy(ids.begin(), ids.end(), rResult.begin());
}

template <>
void RansEvmVmsMonolithicAdjointWallCondition<2>::EquationIdArray(EquationIdArrayType& rResult,
                                                                  ProcessInfo& rCurrentProcessInfo)
{
    IndexType LocalIndex = 0;
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        rResult[LocalIndex++] =
            this->GetGeometry()[i].GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[i].GetDof(ADJOINT_FLUID_VECTOR_1_Y).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[i].GetDof(ADJOINT_FLUID_SCALAR_1).EquationId();
    }
}

template <>
void RansEvmVmsMonolithicAdjointWallCondition<3>::EquationIdArray(EquationIdArrayType& rResult,
                                                                  ProcessInfo& rCurrentProcessInfo)
{
    IndexType LocalIndex = 0;
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        rResult[LocalIndex++] =
            this->GetGeometry()[i].GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[i].GetDof(ADJOINT_FLUID_VECTOR_1_Y).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[i].GetDof(ADJOINT_FLUID_VECTOR_1_Z).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[i].GetDof(ADJOINT_FLUID_SCALAR_1).EquationId();
    }
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::GetDofList(DofsVectorType& rConditionDofList,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    DofsArrayType dofs;
    this->GetDofArray(dofs, rCurrentProcessInfo);

    if (rConditionDofList.size() != dofs.size())
        rConditionDofList.resize(dofs.size());

    std::copy(dofs.begin(), dofs.end(), rConditionDofList.begin());
}

template <>
void RansEvmVmsMonolithicAdjointWallCondition<2>::GetDofArray(DofsArrayType& rConditionDofList,
                                                              ProcessInfo& rCurrentProcessInfo)
{
    IndexType LocalIndex = 0;
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        rConditionDofList[LocalIndex++] =
            this->GetGeometry()[i].pGetDof(ADJOINT_FLUID_VECTOR_1_X);
        rConditionDofList[LocalIndex++] =
            this->GetGeometry()[i].pGetDof(ADJOINT_FLUID_VECTOR_1_Y);
        rConditionDofList[LocalIndex++] =
            this->GetGeometry()[i].pGetDof(ADJOINT_FLUID_SCALAR_1);
    }
}

template <>
void RansEvmVmsMonolithicAdjointWallCondition<3>::GetDofArray(DofsArrayType& rConditionDofList,
                                                              ProcessInfo& rCurrentProcessInfo)
{
    IndexType LocalIndex = 0;
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        rConditionDofList[LocalIndex++] =
            this->GetGeometry()[i].pGetDof(ADJOINT_FLUID_VECTOR_1_X);
        rConditionDofList[LocalIndex++] =
            this->GetGeometry()[i].pGetDof(ADJOINT_FLUID_VECTOR_1_Y);
        rConditionDofList[LocalIndex++] =
            this->GetGeometry()[i].pGetDof(ADJOINT_FLUID_VECTOR_1_Z);
        rConditionDofList[LocalIndex++] =
            this->GetGeometry()[i].pGetDof(ADJOINT_FLUID_SCALAR_1);
    }
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::GetValuesVector(VectorType& rValues, int Step)
{
    ArrayType values;
    this->GetValuesArray(values, Step);

    if (rValues.size() != values.size())
        rValues.resize(values.size(), false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::GetValuesArray(ArrayType& rValues, int Step)
{
    const GeometryType& r_geometry = this->GetGeometry();
    unsigned int LocalIndex = 0;
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const array_1d<double, 3>& r_adjoint_vector =
            r_geometry[i].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1);
        for (unsigned int dim = 0; dim < TDim; ++dim)
        {
            rValues[LocalIndex++] = r_adjoint_vector[dim];
        }
        rValues[LocalIndex++] =
            r_geometry[i].FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
    }
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    ArrayType values;
    this->GetFirstDerivativesArray(values, Step);
    if (rValues.size() != values.size())
        rValues.resize(values.size(), false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::GetFirstDerivativesArray(ArrayType& rValues,
                                                                              int Step)
{
    for (IndexType i = 0; i < rValues.size(); ++i)
    {
        rValues[i] = 0.0;
    }
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    ArrayType values;
    this->GetSecondDerivativesArray(values, Step);
    if (rValues.size() != values.size())
        rValues.resize(values.size(), false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::GetSecondDerivativesArray(ArrayType& rValues,
                                                                               int Step)
{
    const GeometryType& r_geometry = this->GetGeometry();
    unsigned int LocalIndex = 0;
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const array_1d<double, 3>& r_adjoint_vector =
            r_geometry[i].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
        for (unsigned int dim = 0; dim < TDim; ++dim)
        {
            rValues[LocalIndex++] = r_adjoint_vector[dim];
        }
        rValues[LocalIndex++] = 0.0;
    }
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> local_matrix;
    this->CalculateFirstDerivativesLHS(local_matrix, rCurrentProcessInfo);
    if (rLeftHandSideMatrix.size1() != local_matrix.size1() ||
        rLeftHandSideMatrix.size2() != local_matrix.size2())
        rLeftHandSideMatrix.resize(local_matrix.size1(), local_matrix.size2(), false);

    noalias(rLeftHandSideMatrix) = local_matrix;
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateFirstDerivativesLHS(
    BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    // TODO:
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::Calculate(
    const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
    {
        // TODO:
    }
    else if (rVariable == RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
    {
        // TODO:
    }
    else if (rVariable == RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE)
    {
        // TODO:
    }
    else
    {
        KRATOS_ERROR << "Unsupported variable "
                     << rVariable.Name() << " requested at RansEvmVmsMonolithicAdjointWallCondition::Calculate.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateConditionResidualTurbulentKineticEnergyDerivatives(
    BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // TODO:

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateConditionResidualTurbulentEnergyDissipationRateDerivatives(
    BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // TODO:

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY)
    {
        BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize> local_matrix;
        this->CalculateResidualShapeSensitivity(local_matrix, rCurrentProcessInfo);
        if (rOutput.size1() != local_matrix.size1() ||
            rOutput.size2() != local_matrix.size2())
            rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rOutput) = local_matrix;
    }
    else
    {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY)
    {
        this->CalculateResidualShapeSensitivity(rOutput, rCurrentProcessInfo);
    }
    else
    {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateResidualShapeSensitivity(
    BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // TODO:

    KRATOS_CATCH("");
}

template <unsigned int TDim>
double RansEvmVmsMonolithicAdjointWallCondition<TDim>::EvaluateInPoint(
    const Variable<double>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return RansCalculationUtilities().EvaluateInPoint(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TDim>
typename RansEvmVmsMonolithicAdjointWallCondition<TDim>::MatrixType RansEvmVmsMonolithicAdjointWallCondition<TDim>::GetJacobian(
    GeometryData::IntegrationMethod QuadratureOrder, unsigned int IntegrationPointIndex) const
{
    const auto& r_geometry = this->GetGeometry();
    const auto& rDN_De =
        r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex, QuadratureOrder);
    MatrixType jacobian(r_geometry.WorkingSpaceDimension(),
                        r_geometry.LocalSpaceDimension());
    MatrixType coordinates(r_geometry.WorkingSpaceDimension(), r_geometry.PointsNumber());

    for (unsigned int i = 0; i < r_geometry.PointsNumber(); i++)
    {
        const auto& r_coordinates = r_geometry[i].Coordinates();
        for (unsigned int d = 0; d < r_geometry.WorkingSpaceDimension(); d++)
        {
            coordinates(d, i) = r_coordinates[d];
        }
    }

    noalias(jacobian) = prod(coordinates, rDN_De);
    return jacobian;
}

///@}

///@} addtogroup block

// template instantiations
template class RansEvmVmsMonolithicAdjointWallCondition<2>;
template class RansEvmVmsMonolithicAdjointWallCondition<3>;
} // namespace Kratos.
