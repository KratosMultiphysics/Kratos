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
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>

// Project includes
#include "custom_conditions/sbm_fluid_gap_interface_condition.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "utilities/element_size_calculator.h"
#include "utilities/math_utils.h"

namespace Kratos {

namespace {

Point ConnectorPoint(const Matrix& rCoordinates, const double Coordinate)
{
    const double first_weight = 0.5 * (1.0 - Coordinate);
    const double second_weight = 0.5 * (1.0 + Coordinate);
    return Point(
        first_weight * rCoordinates(0, 0) +
            second_weight * rCoordinates(1, 0),
        first_weight * rCoordinates(0, 1) +
            second_weight * rCoordinates(1, 1),
        first_weight * rCoordinates(0, 2) +
            second_weight * rCoordinates(1, 2));
}

} // unnamed namespace

Condition::Pointer SbmFluidGapInterfaceCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmFluidGapInterfaceCondition2D>(
        NewId, pGeometry, pProperties);
}

Condition::Pointer SbmFluidGapInterfaceCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmFluidGapInterfaceCondition2D>(
        NewId,
        Kratos::make_shared<Geometry<Node>>(rThisNodes),
        pProperties);
}

void SbmFluidGapInterfaceCondition2D::Initialize(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    Condition::Initialize(rCurrentProcessInfo);
    const auto& r_neighbours = GetValue(NEIGHBOUR_ELEMENTS);
    const Matrix& r_connector =
        GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
    const auto& r_integration_points =
        LineGaussLegendreIntegrationPoints2::IntegrationPoints();
    mConstitutiveLaws.resize(2 * r_integration_points.size());

    for (IndexType g = 0; g < r_integration_points.size(); ++g) {
        const Point global_point = ConnectorPoint(
            r_connector, r_integration_points[g].X());
        for (IndexType side = 0; side < 2; ++side) {
            const auto& r_parent_geometry =
                r_neighbours[side].GetGeometry();
            Vector N;
            Matrix DN_DX;
            CalculateParentShapeFunctions(
                r_parent_geometry, global_point, N, DN_DX);
            const IndexType constitutive_index = 2 * g + side;
            mConstitutiveLaws[constitutive_index] =
                GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLaws[constitutive_index]->InitializeMaterial(
                GetProperties(), r_parent_geometry, N);
        }
    }
    KRATOS_CATCH("")
}

void SbmFluidGapInterfaceCondition2D::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        CalculateAll(
            rLeftHandSideMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            true,
            true);
    } else {
        const SizeType local_size = GetGeometry().size() * 3;
        rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
        rRightHandSideVector = ZeroVector(local_size);
    }
}

void SbmFluidGapInterfaceCondition2D::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        VectorType rhs;
        CalculateAll(
            rLeftHandSideMatrix, rhs, rCurrentProcessInfo, true, false);
    } else {
        rLeftHandSideMatrix = ZeroMatrix(
            GetGeometry().size() * 3, GetGeometry().size() * 3);
    }
}

void SbmFluidGapInterfaceCondition2D::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        MatrixType lhs;
        CalculateAll(
            lhs, rRightHandSideVector, rCurrentProcessInfo, false, true);
    } else {
        rRightHandSideVector = ZeroVector(GetGeometry().size() * 3);
    }
}

void SbmFluidGapInterfaceCondition2D::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        const SizeType local_size = GetGeometry().size() * 3;
        rDampingMatrix = ZeroMatrix(local_size, local_size);
        if (rRightHandSideVector.size() != local_size) {
            rRightHandSideVector = ZeroVector(local_size);
        }
    } else {
        CalculateAll(
            rDampingMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            true,
            true);
    }
}

void SbmFluidGapInterfaceCondition2D::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        const SizeType local_size = GetGeometry().size() * 3;
        rDampingMatrix = ZeroMatrix(local_size, local_size);
    } else {
        VectorType rhs;
        CalculateAll(rDampingMatrix, rhs, rCurrentProcessInfo, true, false);
    }
}

void SbmFluidGapInterfaceCondition2D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    const auto& r_union_geometry = GetGeometry();
    const auto& r_neighbours = GetValue(NEIGHBOUR_ELEMENTS);
    const Matrix& r_connector =
        GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
    const SizeType number_of_union_nodes = r_union_geometry.size();
    const SizeType local_size = number_of_union_nodes * 3;
    MatrixType local_lhs = ZeroMatrix(local_size, local_size);
    VectorType local_rhs = ZeroVector(local_size);
    const auto& r_integration_points =
        LineGaussLegendreIntegrationPoints2::IntegrationPoints();
    KRATOS_ERROR_IF(
        mConstitutiveLaws.size() != 2 * r_integration_points.size())
        << Info() << " #" << Id() << " has not been initialized."
        << std::endl;

    std::unordered_map<std::size_t, std::size_t> union_local_index;
    union_local_index.reserve(number_of_union_nodes);
    for (IndexType i = 0; i < number_of_union_nodes; ++i) {
        union_local_index.emplace(r_union_geometry[i].Id(), i);
    }

    const double tangent_x = r_connector(1, 0) - r_connector(0, 0);
    const double tangent_y = r_connector(1, 1) - r_connector(0, 1);
    const double connector_length = std::hypot(tangent_x, tangent_y);
    KRATOS_ERROR_IF(connector_length <= std::numeric_limits<double>::epsilon())
        << Info() << " #" << Id() << " has a zero-length connector."
        << std::endl;
    const array_1d<double, 2> normal{
        tangent_y / connector_length, -tangent_x / connector_length};
    const double penalty_factor = GetProperties().Has(PENALTY_COEFFICIENT)
        ? GetProperties()[PENALTY_COEFFICIENT]
        : 0.0;
    // The interior-penalty scaling is based on the background-cell size, not
    // on the connector length. The latter can become arbitrarily small when a
    // projected vertex happens to lie close to the surrogate boundary and
    // would otherwise produce an artificial conditioning problem.
    const double first_parent_size =
        ElementSizeCalculator<2, 4>::MinimumElementSize(
            r_neighbours[0].GetGeometry());
    const double second_parent_size =
        ElementSizeCalculator<2, 4>::MinimumElementSize(
            r_neighbours[1].GetGeometry());
    const double characteristic_length =
        std::min(first_parent_size, second_parent_size);
    KRATOS_ERROR_IF(
        characteristic_length <= std::numeric_limits<double>::epsilon())
        << Info() << " #" << Id()
        << " has a parent element with zero characteristic length."
        << std::endl;
    const double penalty = penalty_factor / characteristic_length;

    struct SideData
    {
        const Geometry<Node>* pGeometry = nullptr;
        Vector N;
        Matrix DN_DX;
        Matrix ConstitutiveB;
        Matrix StressTensor;
        array_1d<double, 2> Velocity = ZeroVector(2);
        double Pressure = 0.0;
        std::vector<std::size_t> UnionIndices;
    };

    for (IndexType g = 0; g < r_integration_points.size(); ++g) {
        const Point global_point = ConnectorPoint(
            r_connector, r_integration_points[g].X());
        const double integration_weight =
            0.5 * connector_length * r_integration_points[g].Weight();
        std::array<SideData, 2> sides;

        for (IndexType side = 0; side < 2; ++side) {
            auto& r_side = sides[side];
            r_side.pGeometry = &r_neighbours[side].GetGeometry();
            const SizeType number_of_parent_nodes = r_side.pGeometry->size();
            CalculateParentShapeFunctions(
                *r_side.pGeometry,
                global_point,
                r_side.N,
                r_side.DN_DX);

            Matrix B = ZeroMatrix(3, number_of_parent_nodes * 2);
            CalculateB(B, r_side.DN_DX);
            Vector nodal_velocity(number_of_parent_nodes * 2);
            r_side.UnionIndices.resize(number_of_parent_nodes);
            for (IndexType i = 0; i < number_of_parent_nodes; ++i) {
                const auto& r_node = (*r_side.pGeometry)[i];
                r_side.UnionIndices[i] = union_local_index.at(r_node.Id());
                const auto& r_velocity =
                    r_node.GetSolutionStepValue(VELOCITY);
                nodal_velocity[2 * i] = r_velocity[0];
                nodal_velocity[2 * i + 1] = r_velocity[1];
                r_side.Velocity[0] += r_side.N[i] * r_velocity[0];
                r_side.Velocity[1] += r_side.N[i] * r_velocity[1];
                r_side.Pressure += r_side.N[i] *
                    r_node.GetSolutionStepValue(PRESSURE);
            }

            ConstitutiveVariables variables;
            ConstitutiveLaw::Parameters parameters(
                *r_side.pGeometry,
                GetProperties(),
                rCurrentProcessInfo);
            Flags& r_options = parameters.GetOptions();
            r_options.Set(
                ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
            r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            r_options.Set(
                ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
            noalias(variables.StrainVector) = prod(B, nodal_velocity);
            parameters.SetShapeFunctionsValues(r_side.N);
            parameters.SetStrainVector(variables.StrainVector);
            parameters.SetStressVector(variables.StressVector);
            parameters.SetConstitutiveMatrix(variables.ConstitutiveMatrix);
            mConstitutiveLaws[2 * g + side]->
                CalculateMaterialResponseCauchy(parameters);

            r_side.ConstitutiveB = prod(
                parameters.GetConstitutiveMatrix(), B);
            const Vector& r_stress = parameters.GetStressVector();
            r_side.StressTensor = ZeroMatrix(2, 2);
            r_side.StressTensor(0, 0) = r_stress[0];
            r_side.StressTensor(1, 1) = r_stress[1];
            r_side.StressTensor(0, 1) = r_stress[2];
            r_side.StressTensor(1, 0) = r_stress[2];
        }

        const array_1d<double, 2> jump_velocity =
            sides[0].Velocity - sides[1].Velocity;
        array_1d<double, 2> average_traction = ZeroVector(2);
        for (IndexType side = 0; side < 2; ++side) {
            average_traction += 0.5 * prod(
                sides[side].StressTensor, normal);
            average_traction -=
                0.5 * sides[side].Pressure * normal;
        }

        for (IndexType test_side = 0; test_side < 2; ++test_side) {
            const auto& r_test = sides[test_side];
            const double test_jump_sign = test_side == 0 ? 1.0 : -1.0;
            Matrix test_stress = ZeroMatrix(2, 2);
            for (IndexType i = 0; i < r_test.N.size(); ++i) {
                const IndexType union_i = r_test.UnionIndices[i];
                for (IndexType i_dim = 0; i_dim < 2; ++i_dim) {
                    const IndexType velocity_row = 3 * union_i + i_dim;
                    BuildStressFromVoigtColumn(
                        test_stress,
                        r_test.ConstitutiveB,
                        2 * i + i_dim);
                    const Vector test_traction = prod(test_stress, normal);

                    const double velocity_residual =
                        -test_jump_sign * r_test.N[i] *
                            average_traction[i_dim] +
                        0.5 * inner_prod(test_traction, jump_velocity) +
                        penalty * test_jump_sign * r_test.N[i] *
                            jump_velocity[i_dim];
                    local_rhs[velocity_row] -=
                        velocity_residual * integration_weight;

                    for (IndexType trial_side = 0; trial_side < 2;
                         ++trial_side) {
                        const auto& r_trial = sides[trial_side];
                        const double trial_jump_sign =
                            trial_side == 0 ? 1.0 : -1.0;
                        Matrix trial_stress = ZeroMatrix(2, 2);
                        for (IndexType j = 0; j < r_trial.N.size(); ++j) {
                            const IndexType union_j =
                                r_trial.UnionIndices[j];
                            for (IndexType j_dim = 0; j_dim < 2; ++j_dim) {
                                BuildStressFromVoigtColumn(
                                    trial_stress,
                                    r_trial.ConstitutiveB,
                                    2 * j + j_dim);
                                const Vector trial_traction =
                                    prod(trial_stress, normal);
                                local_lhs(
                                    velocity_row,
                                    3 * union_j + j_dim) +=
                                    (-0.5 * test_jump_sign * r_test.N[i] *
                                         trial_traction[i_dim] +
                                     0.5 * trial_jump_sign * r_trial.N[j] *
                                         test_traction[j_dim]) *
                                    integration_weight;
                            }
                            local_lhs(
                                velocity_row,
                                3 * union_j + i_dim) +=
                                penalty * test_jump_sign * trial_jump_sign *
                                r_test.N[i] * r_trial.N[j] *
                                integration_weight;
                            local_lhs(
                                velocity_row,
                                3 * union_j + 2) +=
                                0.5 * test_jump_sign * r_test.N[i] *
                                r_trial.N[j] * normal[i_dim] *
                                integration_weight;
                        }
                    }
                }

                const IndexType pressure_row = 3 * union_i + 2;
                local_rhs[pressure_row] +=
                    0.5 * r_test.N[i] *
                    inner_prod(normal, jump_velocity) * integration_weight;
                for (IndexType trial_side = 0; trial_side < 2;
                     ++trial_side) {
                    const auto& r_trial = sides[trial_side];
                    const double trial_jump_sign =
                        trial_side == 0 ? 1.0 : -1.0;
                    for (IndexType j = 0; j < r_trial.N.size(); ++j) {
                        const IndexType union_j = r_trial.UnionIndices[j];
                        for (IndexType j_dim = 0; j_dim < 2; ++j_dim) {
                            local_lhs(
                                pressure_row,
                                3 * union_j + j_dim) -=
                                0.5 * r_test.N[i] * trial_jump_sign *
                                r_trial.N[j] * normal[j_dim] *
                                integration_weight;
                        }
                    }
                }
            }
        }
    }

    if (CalculateStiffnessMatrixFlag) {
        rLeftHandSideMatrix.swap(local_lhs);
    }
    if (CalculateResidualVectorFlag) {
        rRightHandSideVector.swap(local_rhs);
    }
    KRATOS_CATCH("")
}

void SbmFluidGapInterfaceCondition2D::EquationIdVector(
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

void SbmFluidGapInterfaceCondition2D::GetDofList(
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

void SbmFluidGapInterfaceCondition2D::GetFirstDerivativesVector(
    Vector& rValues,
    const int Step) const
{
    const auto& r_geometry = GetGeometry();
    rValues.resize(r_geometry.size() * 3, false);
    IndexType index = 0;
    for (const auto& r_node : r_geometry) {
        const auto& r_velocity =
            r_node.FastGetSolutionStepValue(VELOCITY, Step);
        rValues[index++] = r_velocity[0];
        rValues[index++] = r_velocity[1];
        rValues[index++] =
            r_node.FastGetSolutionStepValue(PRESSURE, Step);
    }
}

void SbmFluidGapInterfaceCondition2D::GetSecondDerivativesVector(
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

int SbmFluidGapInterfaceCondition2D::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    KRATOS_ERROR_IF(Id() < 1)
        << "Condition found with Id " << Id() << std::endl;
    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_ELEMENTS))
        << Info() << " #" << Id() << " is missing NEIGHBOUR_ELEMENTS."
        << std::endl;
    const auto& r_neighbours = GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(r_neighbours.size() != 2)
        << Info() << " #" << Id() << " requires exactly two parent elements."
        << std::endl;
    for (IndexType side = 0; side < 2; ++side) {
        KRATOS_ERROR_IF(r_neighbours[side].GetGeometry().size() != 4)
            << Info() << " #" << Id() << " requires two Quad4 parents."
            << std::endl;
    }
    KRATOS_ERROR_IF_NOT(Has(SURROGATE_BOUNDARY_FACE_COORDINATES))
        << Info() << " #" << Id()
        << " is missing its connector coordinates." << std::endl;
    const Matrix& r_connector =
        GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
    KRATOS_ERROR_IF(r_connector.size1() != 2 || r_connector.size2() != 3)
        << Info() << " #" << Id() << " expects a 2x3 connector matrix."
        << std::endl;
    KRATOS_ERROR_IF_NOT(
        GetProperties().Has(CONSTITUTIVE_LAW) &&
        GetProperties()[CONSTITUTIVE_LAW] != nullptr)
        << Info() << " #" << Id() << " requires a fluid CONSTITUTIVE_LAW."
        << std::endl;
    const double penalty_factor = GetProperties().Has(PENALTY_COEFFICIENT)
        ? GetProperties()[PENALTY_COEFFICIENT]
        : 0.0;
    KRATOS_ERROR_IF(penalty_factor < 0.0)
        << Info() << " #" << Id()
        << " requires a non-negative penalty factor." << std::endl;

    for (const auto& r_node : GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node);
    }
    return 0;
    KRATOS_CATCH("")
}

void SbmFluidGapInterfaceCondition2D::CalculateParentShapeFunctions(
    const GeometryType& rParentGeometry,
    const Point& rGlobalPoint,
    Vector& rN,
    Matrix& rDN_DX)
{
    GeometryType::CoordinatesArrayType local_coordinates = ZeroVector(3);
    rParentGeometry.PointLocalCoordinates(local_coordinates, rGlobalPoint);
    rParentGeometry.ShapeFunctionsValues(rN, local_coordinates);
    Matrix DN_De;
    rParentGeometry.ShapeFunctionsLocalGradients(DN_De, local_coordinates);
    Matrix jacobian;
    rParentGeometry.Jacobian(jacobian, local_coordinates);
    Matrix inverse_jacobian;
    double determinant_jacobian = 0.0;
    MathUtils<double>::InvertMatrix(
        jacobian, inverse_jacobian, determinant_jacobian);
    KRATOS_ERROR_IF(
        std::abs(determinant_jacobian) <=
        std::numeric_limits<double>::epsilon())
        << "SbmFluidGapInterfaceCondition2D has a singular parent element."
        << std::endl;
    rDN_DX.resize(DN_De.size1(), inverse_jacobian.size2(), false);
    noalias(rDN_DX) = prod(DN_De, inverse_jacobian);
}

void SbmFluidGapInterfaceCondition2D::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX)
{
    noalias(rB) = ZeroMatrix(rB.size1(), rB.size2());
    for (IndexType i = 0; i < rDN_DX.size1(); ++i) {
        rB(0, 2 * i) = rDN_DX(i, 0);
        rB(1, 2 * i + 1) = rDN_DX(i, 1);
        rB(2, 2 * i) = rDN_DX(i, 1);
        rB(2, 2 * i + 1) = rDN_DX(i, 0);
    }
}

void SbmFluidGapInterfaceCondition2D::BuildStressFromVoigtColumn(
    Matrix& rStressTensor,
    const Matrix& rConstitutiveBMatrix,
    const IndexType Column)
{
    noalias(rStressTensor) = ZeroMatrix(2, 2);
    rStressTensor(0, 0) = rConstitutiveBMatrix(0, Column);
    rStressTensor(1, 1) = rConstitutiveBMatrix(1, Column);
    rStressTensor(0, 1) = rConstitutiveBMatrix(2, Column);
    rStressTensor(1, 0) = rConstitutiveBMatrix(2, Column);
}

} // namespace Kratos
