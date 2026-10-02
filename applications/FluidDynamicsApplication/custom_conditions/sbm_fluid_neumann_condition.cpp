//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes
#include <limits>

// Project includes
#include "custom_conditions/sbm_fluid_neumann_condition.h"
#include "custom_utilities/sbm_boundary_integration_utility.h"
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "includes/variables.h"

namespace Kratos {

Condition::Pointer SbmFluidNeumannCondition2D4N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmFluidNeumannCondition2D4N>(
        NewId, pGeometry, pProperties);
}

Condition::Pointer SbmFluidNeumannCondition2D4N::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmFluidNeumannCondition2D4N>(
        NewId, GetGeometry().Create(rThisNodes), pProperties);
}

void SbmFluidNeumannCondition2D4N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType local_size = GetGeometry().size() * 3;
    rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        CalculateNeumannLoad(rRightHandSideVector);
    } else {
        rRightHandSideVector = ZeroVector(local_size);
    }
}

void SbmFluidNeumannCondition2D4N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType local_size = GetGeometry().size() * 3;
    rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
}

void SbmFluidNeumannCondition2D4N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        CalculateNeumannLoad(rRightHandSideVector);
    } else {
        rRightHandSideVector = ZeroVector(GetGeometry().size() * 3);
    }
}

void SbmFluidNeumannCondition2D4N::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType local_size = GetGeometry().size() * 3;
    rDampingMatrix = ZeroMatrix(local_size, local_size);
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        if (rRightHandSideVector.size() != local_size) {
            rRightHandSideVector = ZeroVector(local_size);
        }
    } else {
        CalculateNeumannLoad(rRightHandSideVector);
    }
}

void SbmFluidNeumannCondition2D4N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType local_size = GetGeometry().size() * 3;
    rDampingMatrix = ZeroMatrix(local_size, local_size);
}

void SbmFluidNeumannCondition2D4N::CalculateNeumannLoad(
    VectorType& rRightHandSideVector) const
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const SizeType local_size = r_geometry.size() * 3;
    rRightHandSideVector = ZeroVector(local_size);

    SbmBoundaryIntegrationUtility::CheckConditionData(*this, true);
    for (IndexType g = 0;
         g < SbmBoundaryIntegrationUtility::NumberOfIntegrationPoints;
         ++g) {
        const auto integration_data =
            SbmBoundaryIntegrationUtility::CalculateIntegrationPointData(
                *this, g);
        if (integration_data.Weight <=
            std::numeric_limits<double>::epsilon()) {
            continue;
        }

        const auto traction = GetPrescribedTraction(
            g, integration_data.Normal);

        // TraceShapeFunctions makes the physical trace explicit. For Gap-SBM
        // it is the directly extended Q1 basis evaluated on the reconstructed
        // boundary; no Taylor correction is required here.
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            rRightHandSideVector[3 * i] +=
                integration_data.Weight *
                integration_data.TraceShapeFunctions[i] * traction[0];
            rRightHandSideVector[3 * i + 1] +=
                integration_data.Weight *
                integration_data.TraceShapeFunctions[i] * traction[1];
        }
    }
    KRATOS_CATCH("")
}

array_1d<double, 3> SbmFluidNeumannCondition2D4N::GetPrescribedTraction(
    const SizeType IntegrationPointIndex,
    const array_1d<double, 2>& rNormal) const
{
    if (Has(CAUCHY_STRESS_TENSOR)) {
        const Matrix& r_stresses = GetValue(CAUCHY_STRESS_TENSOR);
        return array_1d<double, 3>{
            r_stresses(IntegrationPointIndex, 0) * rNormal[0] +
                r_stresses(IntegrationPointIndex, 1) * rNormal[1],
            r_stresses(IntegrationPointIndex, 2) * rNormal[0] +
                r_stresses(IntegrationPointIndex, 3) * rNormal[1],
            0.0};
    }
    const auto& r_face_load = GetValue(FACE_LOAD);
    return array_1d<double, 3>{r_face_load[0], r_face_load[1], 0.0};
}

void SbmFluidNeumannCondition2D4N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = GetGeometry();
    rResult.resize(r_geometry.size() * 3, false);
    IndexType index = 0;
    for (const auto& r_node : r_geometry) {
        rResult[index++] = r_node.GetDof(VELOCITY_X).EquationId();
        rResult[index++] = r_node.GetDof(VELOCITY_Y).EquationId();
        rResult[index++] = r_node.GetDof(PRESSURE).EquationId();
    }
}

void SbmFluidNeumannCondition2D4N::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = GetGeometry();
    rConditionDofList.clear();
    rConditionDofList.reserve(r_geometry.size() * 3);
    for (const auto& r_node : r_geometry) {
        rConditionDofList.push_back(r_node.pGetDof(VELOCITY_X));
        rConditionDofList.push_back(r_node.pGetDof(VELOCITY_Y));
        rConditionDofList.push_back(r_node.pGetDof(PRESSURE));
    }
}

void SbmFluidNeumannCondition2D4N::GetFirstDerivativesVector(
    Vector& rValues,
    const int Step) const
{
    const auto& r_geometry = GetGeometry();
    rValues.resize(r_geometry.size() * 3, false);
    IndexType index = 0;
    for (const auto& r_node : r_geometry) {
        const auto& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY, Step);
        rValues[index++] = r_velocity[0];
        rValues[index++] = r_velocity[1];
        rValues[index++] = r_node.FastGetSolutionStepValue(PRESSURE, Step);
    }
}

void SbmFluidNeumannCondition2D4N::GetSecondDerivativesVector(
    Vector& rValues,
    const int Step) const
{
    const auto& r_geometry = GetGeometry();
    rValues.resize(r_geometry.size() * 3, false);
    IndexType index = 0;
    for (const auto& r_node : r_geometry) {
        const auto& r_acceleration =
            r_node.FastGetSolutionStepValue(ACCELERATION, Step);
        rValues[index++] = r_acceleration[0];
        rValues[index++] = r_acceleration[1];
        rValues[index++] = 0.0;
    }
}

int SbmFluidNeumannCondition2D4N::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    const int base_check = Condition::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF(
        GetGeometry().size() != 4 ||
        GetGeometry().LocalSpaceDimension() != 2)
        << Info() << " #" << Id() << " requires its owner Quad4 geometry."
        << std::endl;
    SbmBoundaryIntegrationUtility::CheckConditionData(*this, true);
    KRATOS_ERROR_IF_NOT(Has(CAUCHY_STRESS_TENSOR) || Has(FACE_LOAD))
        << Info() << " #" << Id()
        << " requires CAUCHY_STRESS_TENSOR or FACE_LOAD." << std::endl;
    if (Has(CAUCHY_STRESS_TENSOR)) {
        const Matrix& r_stresses = GetValue(CAUCHY_STRESS_TENSOR);
        KRATOS_ERROR_IF(r_stresses.size1() != 2 || r_stresses.size2() != 4)
            << Info() << " #" << Id()
            << " expects CAUCHY_STRESS_TENSOR as a 2x4 matrix."
            << std::endl;
    }

    for (const auto& r_node : GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node);
    }
    return base_check;
    KRATOS_CATCH("")
}

} // namespace Kratos
