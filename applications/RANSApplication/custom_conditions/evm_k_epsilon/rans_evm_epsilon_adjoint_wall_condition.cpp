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
#include "rans_evm_epsilon_adjoint_wall_condition.h"
#include "rans_application_variables.h"
#include "custom_elements/stabilized_convection_diffusion_reaction_adjoint_utilities.h"

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
RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>& RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::operator=(
    RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim> const& rOther)
{
    Condition::operator=(rOther);

    return *this;
}

template <unsigned int TNumNodes, unsigned int TDim>
Condition::Pointer RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmEpsilonAdjointWallCondition>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TNumNodes, unsigned int TDim>
Condition::Pointer RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmEpsilonAdjointWallCondition>(NewId, pGeom, pProperties);
}

template <unsigned int TNumNodes, unsigned int TDim>
Condition::Pointer RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition =
        Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmEpsilonAdjointWallCondition::"
                    "CalculateLocalSystem method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

    rLeftHandSideMatrix.clear();
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateLeftHandSide(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.clear();
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmEpsilonAdjointWallCondition::"
                    "CalculateRightHandSide method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "RansEvmEpsilonAdjointWallCondition::"
                    "CalculateDampingMatrix method not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
int RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

    return Check;

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    std::array<std::size_t, TNumNodes> ids;
    this->EquationIdArray(ids, rCurrentProcessInfo);

    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);
    std::copy(ids.begin(), ids.end(), rResult.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::EquationIdArray(
    std::array<std::size_t, TNumNodes>& rResult, ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        rResult[i] = this->GetGeometry()[i].GetDof(RANS_SCALAR_2_ADJOINT_1).EquationId();
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetDofList(DofsVectorType& rConditionDofList,
                                                                     ProcessInfo& rCurrentProcessInfo)
{
    std::array<Dof<double>::Pointer, TNumNodes> dofs;
    this->GetDofArray(dofs, rCurrentProcessInfo);

    if (rConditionDofList.size() != TNumNodes)
        rConditionDofList.resize(TNumNodes);

    std::copy(dofs.begin(), dofs.end(), rConditionDofList.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetDofArray(
    std::array<Dof<double>::Pointer, TNumNodes>& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        rConditionDofList[i] = this->GetGeometry()[i].pGetDof(RANS_SCALAR_2_ADJOINT_1);
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetValuesVector(VectorType& rValues,
                                                                          int Step)
{
    std::array<double, TNumNodes> values;
    this->GetValuesArray(values, Step);

    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetValuesArray(
    std::array<double, TNumNodes>& rValues, int Step)
{
    const GeometryType& r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rValues[i] = r_geometry[i].FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_1, Step);
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    std::array<double, TNumNodes> values;
    this->GetFirstDerivativesArray(values, Step);
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetFirstDerivativesArray(
    std::array<double, TNumNodes>& rValues, int Step)
{
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rValues[i] = 0.0;
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    std::array<double, TNumNodes> values;
    this->GetSecondDerivativesArray(values, Step);
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);
    std::copy(values.begin(), values.end(), rValues.begin());
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetSecondDerivativesArray(
    std::array<double, TNumNodes>& rValues, int Step)
{
    const GeometryType& r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rValues[i] = r_geometry[i].FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_3, Step);
    }
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
    this->CalculateFirstDerivativesLHS(local_matrix, rCurrentProcessInfo);
    if (rLeftHandSideMatrix.size1() != local_matrix.size1() ||
        rLeftHandSideMatrix.size2() != local_matrix.size2())
        rLeftHandSideMatrix.resize(local_matrix.size1(), local_matrix.size2(), false);

    noalias(rLeftHandSideMatrix) = local_matrix;
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateFirstDerivativesLHS(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    CalculateConditionResidualTurbulentEnergyDissipationRateDerivatives(
        rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::Calculate(
    const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculateConditionResidualTurbulentKineticEnergyDerivatives(
            local_matrix, rCurrentProcessInfo);
        if (rOutput.size1() != local_matrix.size1() ||
            rOutput.size2() != local_matrix.size2())
            rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rOutput) = local_matrix;
    }
    else if (rVariable == RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculateConditionResidualTurbulentEnergyDissipationRateDerivatives(
            local_matrix, rCurrentProcessInfo);
        if (rOutput.size1() != local_matrix.size1() ||
            rOutput.size2() != local_matrix.size2())
            rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rOutput) = local_matrix;
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
                     << rVariable.Name() << " requested at RansEvmEpsilonAdjointWallCondition::Calculate.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateConditionResidualTurbulentKineticEnergyDerivatives(
    BoundedMatrixNN& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rOutput.clear();

    const GeometryType& r_geometry = this->GetGeometry();

    if (!this->Is(SLIP))
        return;

    // Get Shape function data
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const IndexType num_gauss_points = integration_points.size();
    MatrixType shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    const double area = r_geometry.DomainSize();

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    double J = (TNumNodes == 2) ? 0.5 * area : 2.0 * area;

    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    const double eps = std::numeric_limits<double>::epsilon();

    BoundedVector<double, TNumNodes> turbulent_kinematic_viscosity_tke_nodal_sensitivities;

    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = this->GetGeometry()[i_node];

        const Vector& turbulent_kinematic_viscosity_sensitivities =
            r_node.GetValue(RANS_NUT_SCALAR_PARTIAL_DERIVATIVES);

        KRATOS_ERROR_IF(turbulent_kinematic_viscosity_sensitivities.size() != 2) << "RANS_NUT_SCALAR_PARTIAL_DERIVATIVES variable is not specified for node "
                                                                                 << r_node
                                                                                        .Info()
                                                                                 << "\n Please use available NutKEpsilonHighReSensitivitiesProcess to calculate RANS_NUT_SCALAR_PARTIAL_DERIVATIVES.\n";

        turbulent_kinematic_viscosity_tke_nodal_sensitivities[i_node] =
            turbulent_kinematic_viscosity_sensitivities[0] / epsilon_sigma;
    }

    for (unsigned int g = 0; g < num_gauss_points; ++g)
    {
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const double weight = J * integration_points[g].Weight();

        BoundedVector<double, TNumNodes> turbulent_kinematic_viscosity_tke_gauss_sensitivities;
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateGaussSensitivities(
            turbulent_kinematic_viscosity_tke_gauss_sensitivities,
            turbulent_kinematic_viscosity_tke_nodal_sensitivities, gauss_shape_functions);

        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t =
            this->EvaluateInPoint(TURBULENT_VISCOSITY, gauss_shape_functions);
        const double tke = std::max(
            this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions), 0.0);
        const double epsilon = this->EvaluateInPoint(
            TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
        const double effective_kinematic_viscosity = nu + nu_t / epsilon_sigma;

        if (y_plus > eps)
        {
            const double coeff_1 = weight * c_mu_25 * epsilon / (y_plus * nu);
            const double sqrt_tke = std::sqrt(tke);

            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double value = 0.0;

                    value += gauss_shape_functions[a] *
                             turbulent_kinematic_viscosity_tke_gauss_sensitivities[c] * sqrt_tke;

                    value += gauss_shape_functions[a] * gauss_shape_functions[c] *
                             effective_kinematic_viscosity / (2.0 * sqrt_tke);

                    rOutput(c, a) += coeff_1 * value;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateConditionResidualTurbulentEnergyDissipationRateDerivatives(
    BoundedMatrixNN& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rOutput.clear();

    const GeometryType& r_geometry = this->GetGeometry();

    if (!this->Is(SLIP))
        return;

    // Get Shape function data
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const IndexType num_gauss_points = integration_points.size();
    MatrixType shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    const double area = r_geometry.DomainSize();

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    double J = (TNumNodes == 2) ? 0.5 * area : 2.0 * area;

    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    const double eps = std::numeric_limits<double>::epsilon();

    BoundedVector<double, TNumNodes> turbulent_kinematic_viscosity_epsilon_nodal_sensitivities;

    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = this->GetGeometry()[i_node];
        const Vector& turbulent_kinematic_viscosity_sensitivities =
            r_node.GetValue(RANS_NUT_SCALAR_PARTIAL_DERIVATIVES);

        KRATOS_ERROR_IF(turbulent_kinematic_viscosity_sensitivities.size() != 2) << "RANS_NUT_SCALAR_PARTIAL_DERIVATIVES variable is not specified for node "
                                                                                 << r_node
                                                                                        .Info()
                                                                                 << "\n Please use available NutKEpsilonHighReSensitivitiesProcess to calculate RANS_NUT_SCALAR_PARTIAL_DERIVATIVES.\n";

        turbulent_kinematic_viscosity_epsilon_nodal_sensitivities[i_node] =
            turbulent_kinematic_viscosity_sensitivities[1] / epsilon_sigma;
    }

    for (unsigned int g = 0; g < num_gauss_points; ++g)
    {
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const double weight = J * integration_points[g].Weight();

        BoundedVector<double, TNumNodes> turbulent_kinematic_viscosity_epsilon_gauss_sensitivities;
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateGaussSensitivities(
            turbulent_kinematic_viscosity_epsilon_gauss_sensitivities,
            turbulent_kinematic_viscosity_epsilon_nodal_sensitivities, gauss_shape_functions);

        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t =
            this->EvaluateInPoint(TURBULENT_VISCOSITY, gauss_shape_functions);
        const double tke = std::max(
            this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions), 0.0);
        const double epsilon = this->EvaluateInPoint(
            TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
        const double effective_kinematic_viscosity = nu + nu_t / epsilon_sigma;

        if (y_plus > eps)
        {
            const double coeff_1 = weight * c_mu_25 * std::sqrt(tke) / (y_plus * nu);

            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double value = 0.0;

                    value += gauss_shape_functions[a] * gauss_shape_functions[c] *
                             effective_kinematic_viscosity;

                    value += gauss_shape_functions[a] *
                             turbulent_kinematic_viscosity_epsilon_gauss_sensitivities[c] *
                             epsilon;

                    rOutput(c, a) += coeff_1 * value;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    rLeftHandSideMatrix.clear();
}

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY)
    {
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes> local_matrix;
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

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
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

template <unsigned int TNumNodes, unsigned int TDim>
void RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::CalculateResidualShapeSensitivity(
    BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rOutput.clear();

    const GeometryType& r_geometry = this->GetGeometry();
    const unsigned int local_dimension = r_geometry.LocalSpaceDimension();
    const unsigned int num_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    if (!this->Is(SLIP))
        return;

    // Get Shape function data
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const IndexType num_gauss_points = integration_points.size();
    MatrixType shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    Matrix shape_function_local_gradients(num_nodes, local_dimension);
    Matrix jacobian(dimension, local_dimension);

    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    const double eps = std::numeric_limits<double>::epsilon();

    for (unsigned int g = 0; g < num_gauss_points; ++g)
    {
        const Vector& gauss_shape_functions = row(shape_functions, g);
        noalias(shape_function_local_gradients) =
            r_geometry.ShapeFunctionLocalGradient(g, GeometryData::GI_GAUSS_2);
        noalias(jacobian) = this->GetJacobian(GeometryData::GI_GAUSS_2, g);

        LineSensitivityUtility sensitivity_utility(jacobian, shape_function_local_gradients);

        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t =
            this->EvaluateInPoint(TURBULENT_VISCOSITY, gauss_shape_functions);
        const double tke = std::max(
            this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions), 0.0);
        const double epsilon = this->EvaluateInPoint(
            TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
        const double effective_kinematic_viscosity = nu + nu_t / epsilon_sigma;

        const double weight = integration_points[g].Weight();

        if (y_plus > eps)
        {
            const double coeff_1 = weight * c_mu_25 * std::sqrt(tke) *
                                   effective_kinematic_viscosity * epsilon /
                                   (y_plus * nu);

            for (auto s = ShapeParameter::Sequence(num_nodes, dimension); s; ++s)
            {
                const auto& deriv = s.CurrentValue();
                double jacobian_sensitivity;
                sensitivity_utility.CalculateSensitivity(deriv, jacobian_sensitivity);

                for (unsigned int i = 0; i < num_nodes; ++i)
                {
                    rOutput(deriv.NodeIndex * dimension + deriv.Direction, i) +=
                        coeff_1 * gauss_shape_functions[i] * jacobian_sensitivity;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TNumNodes, unsigned int TDim>
double RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::EvaluateInPoint(
    const Variable<double>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TNumNodes, unsigned int TDim>
typename RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::MatrixType RansEvmEpsilonAdjointWallCondition<TNumNodes, TDim>::GetJacobian(
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
template class RansEvmEpsilonAdjointWallCondition<2>;
template class RansEvmEpsilonAdjointWallCondition<3>;
} // namespace Kratos.
