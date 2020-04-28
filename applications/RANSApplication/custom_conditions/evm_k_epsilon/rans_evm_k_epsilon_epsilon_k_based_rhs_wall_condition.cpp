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
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_evm_k_epsilon_epsilon_k_based_rhs_wall_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>& RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::operator=(
    RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes> const& rOther)
{
    Condition::operator=(rOther);

    return *this;
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmKEpsilonEpsilonKBasedRHSWallCondition>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmKEpsilonEpsilonKBasedRHSWallCondition>(
        NewId, pGeom, pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition =
        Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes);

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    VectorType RHS;
    this->CalculateLocalVelocityContribution(rDampingMatrix, RHS, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes);

    if (rDampingMatrix.size1() != TNumNodes || rDampingMatrix.size2() != TNumNodes)
        rDampingMatrix.resize(TNumNodes, TNumNodes);

    rRightHandSideVector.clear();
    rDampingMatrix.clear();

    if (RansCalculationUtilities::IsWall(*this))
    {
        this->AddLocalVelocityContribution(rDampingMatrix, rRightHandSideVector,
                                           rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
int RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

    const GeometryType& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
    }

    return Check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (IndexType i = 0; i < TNumNodes; ++i)
        rResult[i] = Condition::GetGeometry()[i]
                         .GetDof(TURBULENT_ENERGY_DISSIPATION_RATE)
                         .EquationId();
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo)
{
    if (ConditionDofList.size() != TNumNodes)
        ConditionDofList.resize(TNumNodes);

    for (IndexType i = 0; i < TNumNodes; ++i)
        ConditionDofList[i] =
            Condition::GetGeometry()[i].pGetDof(TURBULENT_ENERGY_DISSIPATION_RATE);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::GetValuesVector(
    VectorType& rValues, int Step)
{
    this->GetFirstDerivativesVector(rValues, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::GetFirstDerivativesVector(
    Vector& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
            TURBULENT_ENERGY_DISSIPATION_RATE, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::GetSecondDerivativesVector(
    Vector& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
            TURBULENT_ENERGY_DISSIPATION_RATE_2, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "RansEvmKEpsilonEpsilonKBasedRHSWallCondition" << TNumNodes << "N";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "RansEvmKEpsilonEpsilonKBasedRHSWallCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::AddLocalVelocityContribution(
    MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geometry = this->GetGeometry();
    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    RansCalculationUtilities::CalculateConditionGeometryData(
        r_geometry, this->GetIntegrationMethod(), gauss_weights, shape_functions);
    const IndexType num_gauss_points = gauss_weights.size();

    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    const double kappa = rCurrentProcessInfo[WALL_VON_KARMAN];
    const double eps = std::numeric_limits<double>::epsilon();

    KRATOS_ERROR_IF(!(this->Has(RANS_Y_PLUS)))
        << "RANS_Y_PLUS value is not set in " << this->Info() << " at "
        << this->GetGeometry() << "\n";

    const double y_plus_limit = rCurrentProcessInfo[RANS_Y_PLUS_LIMIT];
    const double y_plus = std::max(this->GetValue(RANS_Y_PLUS), y_plus_limit);

    if (y_plus > eps)
    {
        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const double nu = RansCalculationUtilities::EvaluateInPoint(
                r_geometry, KINEMATIC_VISCOSITY, gauss_shape_functions);
            const double nu_t = RansCalculationUtilities::EvaluateInPoint(
                r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions);
            const double tke = RansCalculationUtilities::EvaluateInPoint(
                r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);

            const double u_tau = c_mu_25 * std::sqrt(std::max(tke, 0.0));

            const double value = gauss_weights[g] * (nu + nu_t / epsilon_sigma) *
                                 std::pow(u_tau, 5) /
                                 (kappa * std::pow(y_plus * nu, 2));
            noalias(rRightHandSideVector) += gauss_shape_functions * value;
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonKBasedRHSWallCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiation

template class RansEvmKEpsilonEpsilonKBasedRHSWallCondition<2>;
template class RansEvmKEpsilonEpsilonKBasedRHSWallCondition<3>;

} // namespace Kratos.
