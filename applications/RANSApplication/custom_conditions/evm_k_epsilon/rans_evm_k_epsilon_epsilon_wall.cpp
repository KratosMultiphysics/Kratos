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
#include <limits>
#include <cmath>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"

// Application includes
#include "rans_evm_k_epsilon_epsilon_wall.h"

#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>& RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::operator=(
    RansEvmKEpsilonEpsilonWall<TDim, TNumNodes> const& rOther)
{
    Condition::operator=(rOther);

    return *this;
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmKEpsilonEpsilonWall>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmKEpsilonEpsilonWall>(NewId, pGeom, pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition =
        Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::CalculateLocalSystem(
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
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes);

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    VectorType RHS;
    this->CalculateLocalVelocityContribution(rDampingMatrix, RHS, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes);

    if (rDampingMatrix.size1() != TNumNodes || rDampingMatrix.size2() != TNumNodes)
        rDampingMatrix.resize(TNumNodes, TNumNodes);

    rRightHandSideVector.clear();
    rDampingMatrix.clear();

    if (this->Is(SLIP))
    {
        this->AddLocalVelocityContribution(rDampingMatrix, rRightHandSideVector,
                                           rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
int RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
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
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                                   ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (IndexType i = 0; i < TNumNodes; ++i)
        rResult[i] = Condition::GetGeometry()[i]
                         .GetDof(TURBULENT_ENERGY_DISSIPATION_RATE)
                         .EquationId();
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::GetDofList(DofsVectorType& ConditionDofList,
                                                             ProcessInfo& CurrentProcessInfo)
{
    if (ConditionDofList.size() != TNumNodes)
        ConditionDofList.resize(TNumNodes);

    for (IndexType i = 0; i < TNumNodes; ++i)
        ConditionDofList[i] =
            Condition::GetGeometry()[i].pGetDof(TURBULENT_ENERGY_DISSIPATION_RATE);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::GetValuesVector(VectorType& rValues, int Step)
{
    this->GetFirstDerivativesVector(rValues, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step)
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
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step)
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
std::string RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "RansEvmKEpsilonEpsilonWall" << TNumNodes << "N";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "RansEvmKEpsilonEpsilonWall";
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::AddLocalVelocityContribution(
    MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geometry = this->GetGeometry();

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
    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const double weight = J * integration_points[g].Weight();

        const double nu = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions);
        const double tke = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double epsilon = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double y_plus = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, RANS_Y_PLUS, gauss_shape_functions);

        const double u_tau = c_mu_25 * std::sqrt(std::max(tke, 0.0));

        if (y_plus > eps)
        {
            const double value =
                weight * (nu + nu_t / epsilon_sigma) * u_tau / (y_plus * nu);

            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                for (IndexType b = 0; b < TNumNodes; ++b)
                {
                    rDampingMatrix(a, b) -=
                        gauss_shape_functions[a] * gauss_shape_functions[b] * value;
                }
                rRightHandSideVector[a] += value * gauss_shape_functions[a] * epsilon;
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonEpsilonWall<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiation

template class RansEvmKEpsilonEpsilonWall<2>;
template class RansEvmKEpsilonEpsilonWall<3>;

} // namespace Kratos.
