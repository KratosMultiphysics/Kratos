//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// System includes
#include <array>
#include <cmath>
#include <limits>

// Project includes
#include "custom_conditions/sbm_fluid_dirichlet_condition.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

namespace Kratos {

Condition::Pointer SbmFluidDirichletCondition2D4N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmFluidDirichletCondition2D4N>(NewId, pGeometry, pProperties);
}

Condition::Pointer SbmFluidDirichletCondition2D4N::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SbmFluidDirichletCondition2D4N>(
        NewId, GetGeometry().Create(rThisNodes), pProperties);
}

void SbmFluidDirichletCondition2D4N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    Condition::Initialize(rCurrentProcessInfo);

    const auto& r_parent_geometry = GetParentGeometry();
    const Matrix& r_face_coordinates = GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
    constexpr std::array<double, 2> gauss_coordinates{
        -0.57735026918962576451, 0.57735026918962576451};
    mConstitutiveLaws.resize(gauss_coordinates.size());

    for (IndexType g = 0; g < gauss_coordinates.size(); ++g) {
        const double first_weight = 0.5 * (1.0 - gauss_coordinates[g]);
        const double second_weight = 0.5 * (1.0 + gauss_coordinates[g]);
        Point global_point(0.0, 0.0, 0.0);
        for (IndexType d = 0; d < 3; ++d) {
            global_point[d] = first_weight * r_face_coordinates(0, d) +
                second_weight * r_face_coordinates(1, d);
        }
        Vector N;
        Matrix DN_DX;
        CalculateParentShapeFunctions(global_point, N, DN_DX);
        mConstitutiveLaws[g] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLaws[g]->InitializeMaterial(GetProperties(), r_parent_geometry, N);
    }
    KRATOS_CATCH("")
}

void SbmFluidDirichletCondition2D4N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    } else {
        const SizeType local_size = GetParentGeometry().size() * 3;
        rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
        rRightHandSideVector = ZeroVector(local_size);
    }
}

void SbmFluidDirichletCondition2D4N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        VectorType rhs;
        CalculateAll(rLeftHandSideMatrix, rhs, rCurrentProcessInfo, true, false);
    } else {
        const SizeType local_size = GetParentGeometry().size() * 3;
        rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
    }
}

void SbmFluidDirichletCondition2D4N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        MatrixType lhs;
        CalculateAll(lhs, rRightHandSideVector, rCurrentProcessInfo, false, true);
    } else {
        rRightHandSideVector = ZeroVector(GetParentGeometry().size() * 3);
    }
}

void SbmFluidDirichletCondition2D4N::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        const SizeType local_size = GetParentGeometry().size() * 3;
        rDampingMatrix = ZeroMatrix(local_size, local_size);
        if (rRightHandSideVector.size() != local_size) {
            rRightHandSideVector = ZeroVector(local_size);
        }
    } else {
        CalculateAll(rDampingMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }
}

void SbmFluidDirichletCondition2D4N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        const SizeType local_size = GetParentGeometry().size() * 3;
        rDampingMatrix = ZeroMatrix(local_size, local_size);
    } else {
        VectorType rhs;
        CalculateAll(rDampingMatrix, rhs, rCurrentProcessInfo, true, false);
    }
}

void SbmFluidDirichletCondition2D4N::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    const auto& r_parent_geometry = GetParentGeometry();
    const SizeType number_of_nodes = r_parent_geometry.size();
    const SizeType local_size = number_of_nodes * 3;
    MatrixType local_lhs = ZeroMatrix(local_size, local_size);
    VectorType local_rhs = ZeroVector(local_size);

    constexpr std::array<double, 2> gauss_coordinates{
        -0.57735026918962576451, 0.57735026918962576451};
    const Matrix& r_face_coordinates = GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
    const Matrix& r_projections = GetValue(SURROGATE_BOUNDARY_PROJECTION);
    KRATOS_ERROR_IF(r_projections.size1() != gauss_coordinates.size() || r_projections.size2() != 3)
        << Info() << " #" << Id() << " expects one 3D projection row per face integration point." << std::endl;
    KRATOS_ERROR_IF(mConstitutiveLaws.size() != gauss_coordinates.size())
        << Info() << " #" << Id() << " has not been initialized." << std::endl;

    const double tangent_x = r_face_coordinates(1, 0) - r_face_coordinates(0, 0);
    const double tangent_y = r_face_coordinates(1, 1) - r_face_coordinates(0, 1);
    const double face_length = std::hypot(tangent_x, tangent_y);
    KRATOS_ERROR_IF(face_length <= std::numeric_limits<double>::epsilon())
        << Info() << " #" << Id() << " has a zero-length surrogate face." << std::endl;
    const array_1d<double, 2> normal{tangent_y / face_length, -tangent_x / face_length};
    const double integration_weight = 0.5 * face_length;

    const double penalty_factor = GetProperties().Has(PENALTY_COEFFICIENT)
        ? GetProperties()[PENALTY_COEFFICIENT]
        : 0.0;
    const double penalty = penalty_factor / face_length;

    Vector velocity(number_of_nodes * 2);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const auto& r_velocity = r_parent_geometry[i].GetSolutionStepValue(VELOCITY);
        velocity[2 * i] = r_velocity[0];
        velocity[2 * i + 1] = r_velocity[1];
    }

    for (IndexType g = 0; g < gauss_coordinates.size(); ++g) {
        const double first_weight = 0.5 * (1.0 - gauss_coordinates[g]);
        const double second_weight = 0.5 * (1.0 + gauss_coordinates[g]);
        Point surrogate_point(0.0, 0.0, 0.0);
        for (IndexType d = 0; d < 3; ++d) {
            surrogate_point[d] = first_weight * r_face_coordinates(0, d) +
                second_weight * r_face_coordinates(1, d);
        }

        Vector N;
        Matrix DN_DX;
        CalculateParentShapeFunctions(surrogate_point, N, DN_DX);
        Vector shifted_N = N;
        const array_1d<double, 3> shift{
            r_projections(g, 0) - surrogate_point.X(),
            r_projections(g, 1) - surrogate_point.Y(),
            r_projections(g, 2) - surrogate_point.Z()};
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            shifted_N[i] += DN_DX(i, 0) * shift[0] + DN_DX(i, 1) * shift[1];
        }

        Matrix B = ZeroMatrix(3, number_of_nodes * 2);
        CalculateB(B, DN_DX);
        ConstitutiveVariables constitutive_variables;
        ConstitutiveLaw::Parameters constitutive_parameters(
            r_parent_geometry, GetProperties(), rCurrentProcessInfo);
        Flags& r_options = constitutive_parameters.GetOptions();
        r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        noalias(constitutive_variables.StrainVector) = prod(B, velocity);
        constitutive_parameters.SetShapeFunctionsValues(N);
        constitutive_parameters.SetStrainVector(constitutive_variables.StrainVector);
        constitutive_parameters.SetStressVector(constitutive_variables.StressVector);
        constitutive_parameters.SetConstitutiveMatrix(constitutive_variables.ConstitutiveMatrix);
        mConstitutiveLaws[g]->CalculateMaterialResponseCauchy(constitutive_parameters);

        const Vector& r_stress = constitutive_parameters.GetStressVector();
        const Matrix& r_constitutive = constitutive_parameters.GetConstitutiveMatrix();
        const Matrix constitutive_B = prod(r_constitutive, B);
        Matrix stress_tensor = ZeroMatrix(2, 2);
        stress_tensor(0, 0) = r_stress[0];
        stress_tensor(1, 1) = r_stress[1];
        stress_tensor(0, 1) = r_stress[2];
        stress_tensor(1, 0) = r_stress[2];

        const double penalty_weight = penalty * integration_weight;

        double pressure = 0.0;
        array_1d<double, 2> shifted_velocity = ZeroVector(2);
        for (IndexType j = 0; j < number_of_nodes; ++j) {
            pressure += r_parent_geometry[j].GetSolutionStepValue(PRESSURE) * N[j];
            const auto& r_velocity = r_parent_geometry[j].GetSolutionStepValue(VELOCITY);
            shifted_velocity[0] += r_velocity[0] * shifted_N[j];
            shifted_velocity[1] += r_velocity[1] * shifted_N[j];
        }
        const Vector current_traction = prod(stress_tensor, normal);
        const auto prescribed_velocity = GetPrescribedVelocity(g);

        Matrix test_stress = ZeroMatrix(2, 2);
        Matrix trial_stress = ZeroMatrix(2, 2);
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType i_dim = 0; i_dim < 2; ++i_dim) {
                BuildStressFromVoigtColumn(test_stress, constitutive_B, 2 * i + i_dim);
                const Vector test_traction = prod(test_stress, normal);

                for (IndexType j = 0; j < number_of_nodes; ++j) {
                    local_lhs(3 * i + i_dim, 3 * j + i_dim) +=
                        shifted_N[i] * shifted_N[j] * penalty_weight;
                    for (IndexType j_dim = 0; j_dim < 2; ++j_dim) {
                        BuildStressFromVoigtColumn(trial_stress, constitutive_B, 2 * j + j_dim);
                        const Vector trial_traction = prod(trial_stress, normal);
                        local_lhs(3 * i + i_dim, 3 * j + j_dim) -=
                            N[i] * trial_traction[i_dim] * integration_weight;
                        local_lhs(3 * i + i_dim, 3 * j + j_dim) +=
                            shifted_N[j] * test_traction[j_dim] * integration_weight;
                    }
                    local_lhs(3 * i + i_dim, 3 * j + 2) +=
                        N[j] * N[i] * normal[i_dim] * integration_weight;
                    local_lhs(3 * j + 2, 3 * i + i_dim) -=
                        N[j] * shifted_N[i] * normal[i_dim] * integration_weight;
                }

                local_rhs[3 * i + i_dim] -=
                    shifted_N[i] * shifted_velocity[i_dim] * penalty_weight;
                local_rhs[3 * i + i_dim] +=
                    N[i] * current_traction[i_dim] * integration_weight;
                local_rhs[3 * i + i_dim] -=
                    pressure * N[i] * normal[i_dim] * integration_weight;
                local_rhs[3 * i + i_dim] +=
                    shifted_N[i] * prescribed_velocity[i_dim] * penalty_weight;

                BuildStressFromVoigtColumn(trial_stress, constitutive_B, 2 * i + i_dim);
                const Vector trial_traction = prod(trial_stress, normal);
                for (IndexType j_dim = 0; j_dim < 2; ++j_dim) {
                    local_rhs[3 * i + i_dim] -=
                        shifted_velocity[j_dim] * trial_traction[j_dim] * integration_weight;
                    local_rhs[3 * i + i_dim] +=
                        prescribed_velocity[j_dim] * trial_traction[j_dim] * integration_weight;
                }
                local_rhs[3 * i + 2] +=
                    shifted_velocity[i_dim] * N[i] * normal[i_dim] * integration_weight;
                local_rhs[3 * i + 2] -=
                    prescribed_velocity[i_dim] * N[i] * normal[i_dim] * integration_weight;
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

void SbmFluidDirichletCondition2D4N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = GetParentGeometry();
    rResult.resize(r_geometry.size() * 3, false);
    IndexType index = 0;
    for (const auto& r_node : r_geometry) {
        rResult[index++] = r_node.GetDof(VELOCITY_X).EquationId();
        rResult[index++] = r_node.GetDof(VELOCITY_Y).EquationId();
        rResult[index++] = r_node.GetDof(PRESSURE).EquationId();
    }
}

void SbmFluidDirichletCondition2D4N::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = GetParentGeometry();
    rConditionDofList.clear();
    rConditionDofList.reserve(r_geometry.size() * 3);
    for (const auto& r_node : r_geometry) {
        rConditionDofList.push_back(r_node.pGetDof(VELOCITY_X));
        rConditionDofList.push_back(r_node.pGetDof(VELOCITY_Y));
        rConditionDofList.push_back(r_node.pGetDof(PRESSURE));
    }
}

void SbmFluidDirichletCondition2D4N::GetFirstDerivativesVector(Vector& rValues, const int Step) const
{
    const auto& r_geometry = GetParentGeometry();
    rValues.resize(r_geometry.size() * 3, false);
    IndexType index = 0;
    for (const auto& r_node : r_geometry) {
        const auto& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY, Step);
        rValues[index++] = r_velocity[0];
        rValues[index++] = r_velocity[1];
        rValues[index++] = r_node.FastGetSolutionStepValue(PRESSURE, Step);
    }
}

void SbmFluidDirichletCondition2D4N::GetSecondDerivativesVector(Vector& rValues, const int Step) const
{
    const auto& r_geometry = GetParentGeometry();
    rValues.resize(r_geometry.size() * 3, false);
    IndexType index = 0;
    for (const auto& r_node : r_geometry) {
        const auto& r_acceleration = r_node.FastGetSolutionStepValue(ACCELERATION, Step);
        rValues[index++] = r_acceleration[0];
        rValues[index++] = r_acceleration[1];
        rValues[index++] = 0.0;
    }
}

int SbmFluidDirichletCondition2D4N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    const int base_check = Condition::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF(GetGeometry().size() != 4 || GetGeometry().LocalSpaceDimension() != 2)
        << Info() << " #" << Id() << " requires its owner Quad4 geometry." << std::endl;
    KRATOS_ERROR_IF_NOT(Has(SURROGATE_BOUNDARY_FACE_COORDINATES))
        << Info() << " #" << Id() << " is missing SURROGATE_BOUNDARY_FACE_COORDINATES." << std::endl;
    const Matrix& r_face_coordinates = GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
    KRATOS_ERROR_IF(r_face_coordinates.size1() != 2 || r_face_coordinates.size2() != 3)
        << Info() << " #" << Id() << " expects a 2x3 surrogate-face coordinate matrix." << std::endl;
    KRATOS_ERROR_IF_NOT(Has(SURROGATE_BOUNDARY_PROJECTION))
        << Info() << " #" << Id() << " is missing SURROGATE_BOUNDARY_PROJECTION." << std::endl;
    KRATOS_ERROR_IF_NOT(Has(SBM_BOUNDARY_VELOCITIES))
        << Info() << " #" << Id() << " is missing SBM_BOUNDARY_VELOCITIES." << std::endl;
    const Matrix& r_prescribed_velocities = GetValue(SBM_BOUNDARY_VELOCITIES);
    KRATOS_ERROR_IF(r_prescribed_velocities.size1() != 2 || r_prescribed_velocities.size2() != 3)
        << Info() << " #" << Id() << " expects a 2x3 SBM_BOUNDARY_VELOCITIES matrix." << std::endl;
    const double penalty_factor = GetProperties().Has(PENALTY_COEFFICIENT)
        ? GetProperties()[PENALTY_COEFFICIENT]
        : 0.0;
    KRATOS_ERROR_IF(penalty_factor < 0.0)
        << Info() << " #" << Id() << " requires a non-negative penalty factor." << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(CONSTITUTIVE_LAW) && GetProperties()[CONSTITUTIVE_LAW] != nullptr)
        << Info() << " #" << Id() << " requires a fluid CONSTITUTIVE_LAW." << std::endl;

    for (const auto& r_node : GetParentGeometry()) {
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

const SbmFluidDirichletCondition2D4N::GeometryType&
SbmFluidDirichletCondition2D4N::GetParentGeometry() const
{
    return GetGeometry();
}

void SbmFluidDirichletCondition2D4N::CalculateParentShapeFunctions(
    const Point& rGlobalPoint,
    Vector& rN,
    Matrix& rDN_DX) const
{
    const auto& r_geometry = GetParentGeometry();
    GeometryType::CoordinatesArrayType local_coordinates = ZeroVector(3);
    r_geometry.PointLocalCoordinates(local_coordinates, rGlobalPoint);
    r_geometry.ShapeFunctionsValues(rN, local_coordinates);
    Matrix DN_De;
    r_geometry.ShapeFunctionsLocalGradients(DN_De, local_coordinates);
    Matrix J;
    r_geometry.Jacobian(J, local_coordinates);
    Matrix inverse_J;
    double determinant_J = 0.0;
    MathUtils<double>::InvertMatrix(J, inverse_J, determinant_J);
    KRATOS_ERROR_IF(std::abs(determinant_J) <= std::numeric_limits<double>::epsilon())
        << Info() << " #" << Id() << " has a singular owner element." << std::endl;
    rDN_DX.resize(DN_De.size1(), inverse_J.size2(), false);
    noalias(rDN_DX) = prod(DN_De, inverse_J);
}

array_1d<double, 3> SbmFluidDirichletCondition2D4N::GetPrescribedVelocity(
    const SizeType IntegrationPointIndex) const
{
    const Matrix& r_values = GetValue(SBM_BOUNDARY_VELOCITIES);
    return array_1d<double, 3>{
        r_values(IntegrationPointIndex, 0),
        r_values(IntegrationPointIndex, 1),
        r_values(IntegrationPointIndex, 2)};
}

void SbmFluidDirichletCondition2D4N::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    noalias(rB) = ZeroMatrix(rB.size1(), rB.size2());
    for (IndexType i = 0; i < rDN_DX.size1(); ++i) {
        rB(0, 2 * i) = rDN_DX(i, 0);
        rB(1, 2 * i + 1) = rDN_DX(i, 1);
        rB(2, 2 * i) = rDN_DX(i, 1);
        rB(2, 2 * i + 1) = rDN_DX(i, 0);
    }
}

void SbmFluidDirichletCondition2D4N::BuildStressFromVoigtColumn(
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
