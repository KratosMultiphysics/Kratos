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
    const ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.clear();
    if (this->Is(SLIP) && rCurrentProcessInfo[IS_CO_SOLVING_PROCESS_ACTIVE])
    {
        this->ApplyRansBasedWallLawFirstDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
    }
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::Calculate(
    const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
    {
        BoundedMatrix<double, TNumNodes, TFluidLocalSize> local_matrix;
        this->CalculateConditionResidualTurbulentKineticEnergyDerivatives(
            local_matrix, rCurrentProcessInfo);
        if (rOutput.size1() != local_matrix.size1() ||
            rOutput.size2() != local_matrix.size2())
            rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rOutput) = local_matrix;
    }
    else if (rVariable == RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
    {
        BoundedMatrix<double, TNumNodes, TFluidLocalSize> local_matrix;
        this->CalculateConditionResidualTurbulentEnergyDissipationRateDerivatives(
            local_matrix, rCurrentProcessInfo);
        if (rOutput.size1() != local_matrix.size1() ||
            rOutput.size2() != local_matrix.size2())
            rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rOutput) = local_matrix;
    }
    else if (rVariable == RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE)
    {
        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> local_matrix;
        this->CalculateFirstDerivativesLHS(local_matrix, rCurrentProcessInfo);
        if (rOutput.size1() != local_matrix.size1() ||
            rOutput.size2() != local_matrix.size2())
            rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rOutput) = local_matrix;
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

    rOutput.clear();
    if (this->Is(SLIP) && rCurrentProcessInfo[IS_CO_SOLVING_PROCESS_ACTIVE])
    {
        this->ApplyRansBasedWallLawTurbulentKineticEnergyDerivatives(rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateConditionResidualTurbulentEnergyDissipationRateDerivatives(
    BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rOutput.clear();
    if (this->Is(SLIP) && rCurrentProcessInfo[IS_CO_SOLVING_PROCESS_ACTIVE])
    {
        this->ApplyRansBasedWallLawTurbulentEnergyDissipationRateDerivatives(
            rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::ApplyRansBasedWallLawFirstDerivatives(
    BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rLocalMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geometry = this->GetGeometry();

    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const std::size_t number_of_gauss_points = integration_points.size();
    MatrixType shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    array_1d<double, 3> normal;
    this->CalculateNormal(normal); // this already contains the area
    double A = norm_2(normal);

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    double J = (TDim == 2) ? 0.5 * A : 2.0 * A;

    const size_t block_size = TDim + 1;

    const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    const double inv_von_karman = 1.0 / rCurrentProcessInfo[WALL_VON_KARMAN];
    const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
    const double eps = std::numeric_limits<double>::epsilon();

    for (size_t g = 0; g < number_of_gauss_points; ++g)
    {
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const double weight = J * integration_points[g].Weight();

        const double tke =
            this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const array_1d<double, 3>& r_velocity =
            this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
        const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
        const double velocity_magnitude = norm_2(r_velocity);
        const double velocity_magnitude_2 = std::pow(velocity_magnitude, 2);
        const double rho = this->EvaluateInPoint(DENSITY, gauss_shape_functions);

        if (velocity_magnitude > eps)
        {
            BoundedMatrix<double, TNumNodes, TDim> velocity_magnitude_velocity_derivatives;
            this->CalculateVelocityMagnitudeVelocityDerivative(
                velocity_magnitude_velocity_derivatives, velocity_magnitude,
                r_velocity, gauss_shape_functions);

            const double u_plus = inv_von_karman * std::log(y_plus) + beta;
            const double u_tau_tke = c_mu_25 * std::sqrt(std::max(tke, 0.0));
            const double u_tau_vel = velocity_magnitude / u_plus;
            const double u_tau = std::max(u_tau_tke, u_tau_vel);
            const double u_tau_2 = std::pow(u_tau, 2);

            BoundedMatrix<double, TNumNodes, TDim> u_tau_derivatives;
            u_tau_derivatives.clear();
            if (u_tau_tke < u_tau_vel)
            {
                noalias(u_tau_derivatives) = velocity_magnitude_velocity_derivatives / u_plus;
            }

            for (IndexType c = 0; c < TNumNodes; ++c)
            {
                const int c_block = c * block_size;
                for (IndexType a = 0; a < TNumNodes; ++a)
                {
                    const int a_block = a * block_size;
                    for (IndexType k = 0; k < TDim; ++k)
                    {
                        double value = 0.0;

                        value += gauss_shape_functions[a] * gauss_shape_functions[c] *
                                 u_tau_2 / velocity_magnitude;

                        for (IndexType i = 0; i < TDim; ++i)
                        {
                            double value_ik = 0.0;
                            value_ik += gauss_shape_functions[a] * 2 * u_tau *
                                        u_tau_derivatives(c, k) *
                                        r_velocity[i] / velocity_magnitude;
                            value_ik -= gauss_shape_functions[a] * u_tau_2 *
                                        velocity_magnitude_velocity_derivatives(c, k) *
                                        r_velocity[i] / velocity_magnitude_2;
                            rLocalMatrix(c_block + k, a_block + i) -=
                                value_ik * weight * rho;
                        }
                        rLocalMatrix(c_block + k, a_block + k) -= value * weight * rho;
                    }
                }
            }
        }
    }
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::ApplyRansBasedWallLawTurbulentKineticEnergyDerivatives(
    BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rLocalMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geometry = this->GetGeometry();

    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const std::size_t number_of_gauss_points = integration_points.size();
    MatrixType shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    array_1d<double, 3> normal;
    this->CalculateNormal(normal); // this already contains the area
    double A = norm_2(normal);

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    double J = (TDim == 2) ? 0.5 * A : 2.0 * A;

    const size_t block_size = TDim + 1;

    const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    const double inv_von_karman = 1.0 / rCurrentProcessInfo[WALL_VON_KARMAN];
    const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
    const double eps = std::numeric_limits<double>::epsilon();

    for (size_t g = 0; g < number_of_gauss_points; ++g)
    {
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const double weight = J * integration_points[g].Weight();

        const double tke = std::max(
            this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions), 0.0);
        const array_1d<double, 3>& r_velocity =
            this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
        const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
        const double velocity_magnitude = norm_2(r_velocity);
        const double rho = this->EvaluateInPoint(DENSITY, gauss_shape_functions);

        if (velocity_magnitude > eps)
        {
            const double u_plus = inv_von_karman * std::log(y_plus) + beta;
            const double u_tau_tke = c_mu_25 * std::sqrt(std::max(tke, 0.0));
            const double u_tau_vel = velocity_magnitude / u_plus;
            const double u_tau = std::max(u_tau_tke, u_tau_vel);

            BoundedVector<double, TNumNodes> u_tau_derivatives;
            u_tau_derivatives.clear();
            if (u_tau_tke >= u_tau_vel)
            {
                noalias(u_tau_derivatives) =
                    gauss_shape_functions * (c_mu_25 * 0.5 / std::sqrt(tke));
            }

            for (IndexType c = 0; c < TNumNodes; ++c)
            {
                for (IndexType a = 0; a < TNumNodes; ++a)
                {
                    const int a_block = a * block_size;
                    for (IndexType i = 0; i < TDim; ++i)
                    {
                        double value_ik = 0.0;

                        value_ik += gauss_shape_functions[a] * 2 * u_tau *
                                    u_tau_derivatives[c] * r_velocity[i] / velocity_magnitude;

                        rLocalMatrix(c, a_block + i) -= value_ik * weight * rho;
                    }
                }
            }
        }
    }
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::ApplyRansBasedWallLawTurbulentEnergyDissipationRateDerivatives(
    BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rLocalMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TFluidLocalSize ||
        rLeftHandSideMatrix.size2() != TFluidLocalSize)
        rLeftHandSideMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);
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

    rOutput.clear();

    GeometryType& r_geometry = this->GetGeometry();
    const unsigned int local_dimension = r_geometry.LocalSpaceDimension();
    const unsigned int num_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Get Shape function data
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const IndexType num_gauss_points = integration_points.size();
    MatrixType shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    Matrix shape_function_local_gradients(num_nodes, local_dimension);
    Matrix jacobian(dimension, local_dimension);

    array_1d<double, 3> normal;
    this->CalculateNormal(normal); // this already contains the area
    double A = norm_2(normal);
    normal /= A;
    double J = (TDim == 2) ? 0.5 * A : 2.0 * A;

    const size_t block_size = TDim + 1;

    const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    const double inv_von_karman = 1.0 / rCurrentProcessInfo[WALL_VON_KARMAN];
    const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
    const double eps = std::numeric_limits<double>::epsilon();

    BoundedMatrix<double, TCoordLocalSize, TDim> unit_normal_derivatives;
    this->CalculateUnitNormalShapeSensitivities(unit_normal_derivatives);

    for (size_t g = 0; g < num_gauss_points; ++g)
    {
        const Vector& gauss_shape_functions = row(shape_functions, g);
        double Weight = J * integration_points[g].Weight();
        noalias(shape_function_local_gradients) =
            r_geometry.ShapeFunctionLocalGradient(g, GeometryData::GI_GAUSS_2);
        noalias(jacobian) = this->GetJacobian(GeometryData::GI_GAUSS_2, g);

        LineSensitivityUtility sensitivity_utility(jacobian, shape_function_local_gradients);

        if (!this->Is(SLIP))
        {
            const double external_pressure =
                this->EvaluateInPoint(EXTERNAL_PRESSURE, gauss_shape_functions);

            for (auto s = ShapeParameter::Sequence(num_nodes, dimension); s; ++s)
            {
                const auto& deriv = s.CurrentValue();
                double J_deriv;
                sensitivity_utility.CalculateSensitivity(deriv, J_deriv);

                for (unsigned int i = 0; i < num_nodes; ++i)
                {
                    for (unsigned int d = 0; d < dimension; ++d)
                    {
                        rOutput(deriv.NodeIndex * dimension + deriv.Direction,
                                i * block_size + d) -=
                            J_deriv * gauss_shape_functions[i] *
                            external_pressure * normal[d];
                        rOutput(deriv.NodeIndex * dimension + deriv.Direction,
                                i * block_size + d) -=
                            Weight * gauss_shape_functions[i] * external_pressure *
                            unit_normal_derivatives(
                                deriv.NodeIndex * dimension + deriv.Direction, d);
                    }
                }
            }
        }
        else if (this->Is(SLIP) && rCurrentProcessInfo[IS_CO_SOLVING_PROCESS_ACTIVE])
        {
            const array_1d<double, 3>& r_velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            const double velocity_magnitude = norm_2(r_velocity);

            const double tke =
                this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
            const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
            const double rho = this->EvaluateInPoint(DENSITY, gauss_shape_functions);

            if (velocity_magnitude > eps)
            {
                const double u_tau = std::max(
                    c_mu_25 * std::sqrt(std::max(tke, 0.0)),
                    velocity_magnitude / (inv_von_karman * std::log(y_plus) + beta));
                const double value = rho * std::pow(u_tau, 2) / velocity_magnitude;

                for (auto s = ShapeParameter::Sequence(num_nodes, dimension); s; ++s)
                {
                    const auto& deriv = s.CurrentValue();
                    double J_deriv;
                    sensitivity_utility.CalculateSensitivity(deriv, J_deriv);

                    for (unsigned int i = 0; i < num_nodes; ++i)
                    {
                        for (unsigned int d = 0; d < dimension; ++d)
                        {
                            rOutput(deriv.NodeIndex * dimension + deriv.Direction,
                                    i * block_size + d) -=
                                J_deriv * value * gauss_shape_functions[i] * r_velocity[d];
                        }
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
double RansEvmVmsMonolithicAdjointWallCondition<TDim>::EvaluateInPoint(
    const Variable<double>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TDim>
array_1d<double, 3> RansEvmVmsMonolithicAdjointWallCondition<TDim>::EvaluateInPoint(
    const Variable<array_1d<double, 3>>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TDim>
void RansEvmVmsMonolithicAdjointWallCondition<TDim>::CalculateVelocityMagnitudeVelocityDerivative(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double VelocityMagnitude,
    const array_1d<double, 3>& rVelocity,
    const Vector& rGaussShapeFunctions) const
{
    if (VelocityMagnitude <= std::numeric_limits<double>::epsilon())
    {
        rOutput.clear();
    }
    else
    {
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) =
                    rVelocity[i_dim] * rGaussShapeFunctions[i_node] / VelocityMagnitude;
    }
}

template <>
void RansEvmVmsMonolithicAdjointWallCondition<2>::CalculateNormal(array_1d<double, 3>& An)
{
    Geometry<Node<3>>& pGeometry = this->GetGeometry();

    An[0] = pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = -(pGeometry[1].X() - pGeometry[0].X());
    An[2] = 0.00;
}

template <>
void RansEvmVmsMonolithicAdjointWallCondition<3>::CalculateNormal(array_1d<double, 3>& An)
{
    Geometry<Node<3>>& pGeometry = this->GetGeometry();

    array_1d<double, 3> v1, v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An, v1, v2);
    An *= 0.5;
}

template <>
void RansEvmVmsMonolithicAdjointWallCondition<2>::CalculateUnitNormalShapeSensitivities(
    BoundedMatrix<double, 4, 2>& rOutput)
{
    const int number_of_nodes = TNumNodes;
    const int dimension = TNumNodes;

    array_1d<double, 3> normal;
    this->CalculateNormal(normal);
    const double A = norm_2(normal);

    BoundedMatrix<double, number_of_nodes * dimension, dimension> normal_shape_sensitivity;
    // normal direction - x
    normal_shape_sensitivity(0, 0) = 0.0; // derivative w.r.t. node 0, direction x
    normal_shape_sensitivity(1, 0) = -1.0; // derivative w.r.t. node 0, direction y
    normal_shape_sensitivity(2, 0) = 0.0; // derivative w.r.t. node 1, direction x
    normal_shape_sensitivity(3, 0) = 1.0; // derivative w.r.t. node 1, direction y

    // normal direction - y
    normal_shape_sensitivity(0, 1) = 1.0; // derivative w.r.t. node 0, direction x
    normal_shape_sensitivity(1, 1) = 0.0; // derivative w.r.t. node 0, direction y
    normal_shape_sensitivity(2, 1) = -1.0; // derivative w.r.t. node 1, direction x
    normal_shape_sensitivity(3, 1) = 0.0; // derivative w.r.t. node 1, direction y

    noalias(rOutput) = normal_shape_sensitivity * (1 / A);

    BoundedVector<double, number_of_nodes * dimension> n_i_n_i_deriv;
    n_i_n_i_deriv.clear();
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const int i_block = i_node * dimension;
        for (int i_dim = 0; i_dim < dimension; ++i_dim)
        {
            for (int normal_dim = 0; normal_dim < dimension; ++normal_dim)
                n_i_n_i_deriv[i_block + i_dim] +=
                    normal[normal_dim] * normal_shape_sensitivity(i_block + i_dim, normal_dim);

            for (int normal_dim = 0; normal_dim < dimension; ++normal_dim)
            {
                rOutput(i_block + i_dim, normal_dim) -=
                    normal[normal_dim] * n_i_n_i_deriv[i_block + i_dim] / std::pow(A, 3);
            }
        }
    }
}

template <>
void RansEvmVmsMonolithicAdjointWallCondition<3>::CalculateUnitNormalShapeSensitivities(
    BoundedMatrix<double, 9, 3>& rOutput)
{
    const int number_of_nodes = TNumNodes;
    const int dimension = TNumNodes;
    const GeometryType& r_geometry = this->GetGeometry();

    array_1d<double, 3> normal;
    this->CalculateNormal(normal);
    const double A = norm_2(normal);

    const double a1 = r_geometry[1].X() - r_geometry[0].X();
    const double a2 = r_geometry[1].Y() - r_geometry[0].Y();
    const double a3 = r_geometry[1].Z() - r_geometry[0].Z();

    const double b1 = r_geometry[2].X() - r_geometry[0].X();
    const double b2 = r_geometry[2].Y() - r_geometry[0].Y();
    const double b3 = r_geometry[2].Z() - r_geometry[0].Z();

    BoundedMatrix<double, number_of_nodes * dimension, dimension> normal_shape_sensitivity;
    // normal direction - x
    normal_shape_sensitivity(0, 0) = 0.0; // derivative w.r.t. node 0, direction x
    normal_shape_sensitivity(1, 0) = -b3 + a3; // derivative w.r.t. node 0, direction y
    normal_shape_sensitivity(2, 0) = -a2 + b2; // derivative w.r.t. node 0, direction z
    normal_shape_sensitivity(3, 0) = 0.0; // derivative w.r.t. node 1, direction x
    normal_shape_sensitivity(4, 0) = b3; // derivative w.r.t. node 1, direction y
    normal_shape_sensitivity(5, 0) = -b2; // derivative w.r.t. node 1, direction z
    normal_shape_sensitivity(6, 0) = 0.0; // derivative w.r.t. node 2, direction x
    normal_shape_sensitivity(7, 0) = -a3; // derivative w.r.t. node 2, direction y
    normal_shape_sensitivity(8, 0) = a2; // derivative w.r.t. node 2, direction z

    // normal direction - y
    normal_shape_sensitivity(0, 1) = -a3 + b3; // derivative w.r.t. node 0, direction x
    normal_shape_sensitivity(1, 1) = 0.0; // derivative w.r.t. node 0, direction y
    normal_shape_sensitivity(2, 1) = -b1 + a1; // derivative w.r.t. node 0, direction z
    normal_shape_sensitivity(3, 1) = -b3; // derivative w.r.t. node 1, direction x
    normal_shape_sensitivity(4, 1) = 0.0; // derivative w.r.t. node 1, direction y
    normal_shape_sensitivity(5, 1) = b1; // derivative w.r.t. node 1, direction z
    normal_shape_sensitivity(6, 1) = a3; // derivative w.r.t. node 2, direction x
    normal_shape_sensitivity(7, 1) = 0.0; // derivative w.r.t. node 2, direction y
    normal_shape_sensitivity(8, 1) = -a1; // derivative w.r.t. node 2, direction z

    // normal direction - z
    normal_shape_sensitivity(0, 2) = -b2 + a2; // derivative w.r.t. node 0, direction x
    normal_shape_sensitivity(1, 2) = -a1 + b1; // derivative w.r.t. node 0, direction y
    normal_shape_sensitivity(2, 2) = 0.0; // derivative w.r.t. node 0, direction z
    normal_shape_sensitivity(3, 2) = b2; // derivative w.r.t. node 1, direction x
    normal_shape_sensitivity(4, 2) = -b1; // derivative w.r.t. node 1, direction y
    normal_shape_sensitivity(5, 2) = 0.0; // derivative w.r.t. node 1, direction z
    normal_shape_sensitivity(6, 2) = -a2; // derivative w.r.t. node 2, direction x
    normal_shape_sensitivity(7, 2) = a1; // derivative w.r.t. node 2, direction y
    normal_shape_sensitivity(8, 2) = 0.0; // derivative w.r.t. node 2, direction z

    noalias(normal_shape_sensitivity) = normal_shape_sensitivity * 0.5;

    noalias(rOutput) = normal_shape_sensitivity * (1 / A);

    BoundedVector<double, number_of_nodes * dimension> n_i_n_i_deriv;
    n_i_n_i_deriv.clear();
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const int i_block = i_node * dimension;
        for (int i_dim = 0; i_dim < dimension; ++i_dim)
        {
            for (int normal_dim = 0; normal_dim < dimension; ++normal_dim)
                n_i_n_i_deriv[i_block + i_dim] +=
                    normal[normal_dim] * normal_shape_sensitivity(i_block + i_dim, normal_dim);

            for (int normal_dim = 0; normal_dim < dimension; ++normal_dim)
            {
                rOutput(i_block + i_dim, normal_dim) -=
                    normal[normal_dim] * n_i_n_i_deriv[i_block + i_dim] / std::pow(A, 3);
            }
        }
    }
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

    for (unsigned int i = 0; i < r_geometry.PointsNumber(); ++i)
    {
        const auto& r_coordinates = r_geometry[i].Coordinates();
        for (unsigned int d = 0; d < r_geometry.WorkingSpaceDimension(); ++d)
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
