//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti,
//                   Ricky Aristio, Alicia Knauer
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/non_linear_bernoulli_beam_element_3D.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace
{

/*
 * Products with the variation containers of a 3x3 matrix M with respect to the n element dofs:
 *  - first variation:  (3 n) x 3     matrix with M,r(i, j)  at (i * n + r, j)
 *  - second variation: (3 n) x (3 n) matrix with M,rs(i, j) at (i * n + r, j * n + s)
 * A product of a first variation with a constant 3x3 matrix from the right is a plain prod(), the
 * other products below are not, because of this layout.
 */

/// rResult += A * V,r for all dofs r, with V,r a first variation
template<class TMatrixType>
void AddMatrixTimesVariation(
    Matrix& rResult,
    const TMatrixType& rMatrix,
    const Matrix& rVariation,
    const std::size_t NumberOfDofs)
{
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t k = 0; k < 3; ++k) {
                for (std::size_t r = 0; r < NumberOfDofs; ++r) {
                    rResult(i * NumberOfDofs + r, j) += rMatrix(i, k) * rVariation(k * NumberOfDofs + r, j);
                }
            }
        }
    }
}

/// rResult += W,rs * A for all dofs r, s, with W,rs a second variation
template<class TMatrixType>
void AddSecondVariationTimesMatrix(
    Matrix& rResult,
    const Matrix& rSecondVariation,
    const TMatrixType& rMatrix,
    const std::size_t NumberOfDofs)
{
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t k = 0; k < 3; ++k) {
                for (std::size_t r = 0; r < NumberOfDofs; ++r) {
                    for (std::size_t s = 0; s < NumberOfDofs; ++s) {
                        rResult(i * NumberOfDofs + r, j * NumberOfDofs + s) +=
                            rSecondVariation(i * NumberOfDofs + r, k * NumberOfDofs + s) * rMatrix(k, j);
                    }
                }
            }
        }
    }
}

/// rResult += A * W,rs for all dofs r, s, with W,rs a second variation
template<class TMatrixType>
void AddMatrixTimesSecondVariation(
    Matrix& rResult,
    const TMatrixType& rMatrix,
    const Matrix& rSecondVariation,
    const std::size_t NumberOfDofs)
{
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t k = 0; k < 3; ++k) {
                for (std::size_t r = 0; r < NumberOfDofs; ++r) {
                    for (std::size_t s = 0; s < NumberOfDofs; ++s) {
                        rResult(i * NumberOfDofs + r, j * NumberOfDofs + s) +=
                            rMatrix(i, k) * rSecondVariation(k * NumberOfDofs + r, j * NumberOfDofs + s);
                    }
                }
            }
        }
    }
}

/// rResult += V,r * W,s for all dofs r, s: the mixed second variation of the product of two matrices V and W
void AddVariationTimesVariation(
    Matrix& rResult,
    const Matrix& rFirstVariation,
    const Matrix& rSecondVariation,
    const std::size_t NumberOfDofs)
{
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t k = 0; k < 3; ++k) {
                for (std::size_t r = 0; r < NumberOfDofs; ++r) {
                    for (std::size_t s = 0; s < NumberOfDofs; ++s) {
                        rResult(i * NumberOfDofs + r, j * NumberOfDofs + s) +=
                            rFirstVariation(i * NumberOfDofs + r, k) * rSecondVariation(k * NumberOfDofs + s, j);
                    }
                }
            }
        }
    }
}

} // namespace

    void NonLinearBernoulliBeamElement3D::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        InitializeMaterial();


        mDofsPerNode = 4;
        mNumberOfDofs  = this->GetGeometry().size() * mDofsPerNode;
        KRATOS_CATCH("")

    }

    void NonLinearBernoulliBeamElement3D::GetFirstDerivativesVector(
        Vector& rValues,
        const int Step) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        if (rValues.size() != number_of_nodes * 4) {
            rValues.resize(number_of_nodes * 4, false);
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 4;

            const auto& r_velocity =
                r_geometry[i].FastGetSolutionStepValue(
                    VELOCITY,
                    Step);

            rValues[index] = r_velocity[0];
            rValues[index + 1] = r_velocity[1];
            rValues[index + 2] = r_velocity[2];
            rValues[index + 3] =
                r_geometry[i].FastGetSolutionStepValue(
                    ANGULAR_VELOCITY_X,
                    Step);
        }
    }

    void NonLinearBernoulliBeamElement3D::GetSecondDerivativesVector(
        Vector& rValues,
        const int Step) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        if (rValues.size() != number_of_nodes * 4) {
            rValues.resize(number_of_nodes * 4, false);
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 4;

            const auto& r_acceleration =
                r_geometry[i].FastGetSolutionStepValue(
                    ACCELERATION,
                    Step);

            rValues[index] = r_acceleration[0];
            rValues[index + 1] = r_acceleration[1];
            rValues[index + 2] = r_acceleration[2];
            rValues[index + 3] =
                r_geometry[i].FastGetSolutionStepValue(
                    ANGULAR_ACCELERATION_X,
                    Step);
        }
    }

    void NonLinearBernoulliBeamElement3D::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType number_of_dofs = number_of_nodes * 4;

        if (rMassMatrix.size1() != number_of_dofs || rMassMatrix.size2() != number_of_dofs) {
            rMassMatrix.resize(number_of_dofs, number_of_dofs, false);
        }
        noalias(rMassMatrix) = ZeroMatrix(number_of_dofs, number_of_dofs);

        const auto& r_properties = GetProperties();
        const double density = r_properties[DENSITY];
        const double mass_per_length = density * r_properties[CROSS_AREA];

        // Polar mass moment of inertia per unit length, for the rotation about the center line
        const double rotational_inertia_per_length =
            density * (r_properties[I22] + r_properties[I33]);

        const auto& r_integration_points = r_geometry.IntegrationPoints();
        const Matrix& r_shape_functions = r_geometry.ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            const Matrix& r_shape_function_derivatives =
                r_geometry.ShapeFunctionDerivatives(1, point_number);

            // Reference base vector A1, so that dL = |A1| dxi
            array_1d<double, 3> reference_base_vector = ZeroVector(3);
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                reference_base_vector +=
                    r_shape_function_derivatives(i, 0) *
                    r_geometry[i].GetInitialPosition().Coordinates();
            }

            const double integration_weight =
                r_integration_points[point_number].Weight() *
                norm_2(reference_base_vector);

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                for (IndexType j = 0; j < number_of_nodes; ++j) {
                    const double shape_function_product =
                        r_shape_functions(point_number, i) *
                        r_shape_functions(point_number, j) *
                        integration_weight;

                    for (IndexType k = 0; k < 3; ++k) {
                        rMassMatrix(i * 4 + k, j * 4 + k) +=
                            shape_function_product * mass_per_length;
                    }

                    rMassMatrix(i * 4 + 3, j * 4 + 3) +=
                        shape_function_product * rotational_inertia_per_length;
                }
            }
        }

        KRATOS_CATCH("")
    }



    void NonLinearBernoulliBeamElement3D::CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables)
    {
        const auto& r_geometry = GetGeometry();

        // Get shape functions and derivatives
        const Matrix& r_shape_functions =
            r_geometry.ShapeFunctionsValues();

        const Matrix& r_shape_function_derivatives =
            r_geometry.ShapeFunctionDerivatives(
                1,
                IntegrationPointIndex);

        const Matrix& r_shape_function_second_derivatives =
            r_geometry.ShapeFunctionDerivatives(
                2,
                IntegrationPointIndex);

        // Clear variables
        rKinematicVariables.R1.clear();
        rKinematicVariables.R2.clear();
        rKinematicVariables.r1.clear();
        rKinematicVariables.r2.clear();

        // Calculate reference and current configurations
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            const array_1d<double, 3>& r_initial_position =
                r_geometry[i].GetInitialPosition();

            const array_1d<double, 3>& r_displacement =
                r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

            // Reference configuration
            rKinematicVariables.R1 +=
                r_shape_function_derivatives(i, 0) *
                r_initial_position;

            rKinematicVariables.R2 +=
                r_shape_function_second_derivatives(i, 0) *
                r_initial_position;

            // Current configuration
            const array_1d<double, 3> current_position =
                r_initial_position + r_displacement;

            rKinematicVariables.r1 +=
                r_shape_function_derivatives(i, 0) *
                current_position;

            rKinematicVariables.r2 +=
                r_shape_function_second_derivatives(i, 0) *
                current_position;

            // Rotations
            const double rotation =
                r_geometry[i].FastGetSolutionStepValue(ROTATION_X);

            rKinematicVariables.phi +=
                r_shape_functions(IntegrationPointIndex, i) *
                rotation;

            rKinematicVariables.phi_der +=
                r_shape_function_derivatives(i, 0) *
                rotation;
        }

        // Calculate length measures
        rKinematicVariables.A =
            norm_2(rKinematicVariables.R1);

        rKinematicVariables.a =
            norm_2(rKinematicVariables.r1);

        // Calculate curvature measures
        const double reference_base_vector_product =
            inner_prod(
                rKinematicVariables.R1,
                rKinematicVariables.R2);

        const double current_base_vector_product =
            inner_prod(
                rKinematicVariables.r1,
                rKinematicVariables.r2);

        rKinematicVariables.B = std::sqrt(
            inner_prod(
                rKinematicVariables.R2,
                rKinematicVariables.R2)
            - std::pow(
                reference_base_vector_product /
                    rKinematicVariables.A,
                2));

        rKinematicVariables.b = std::sqrt(
            inner_prod(
                rKinematicVariables.r2,
                rKinematicVariables.r2)
            - std::pow(
                current_base_vector_product /
                    rKinematicVariables.a,
                2));

        // Compute reference cross-section geometry
        ComputeReferenceCrossSectionGeometry(
            rKinematicVariables);

        // Compute current cross-section geometry
        ComputeCurrentCrossSectionGeometry(
            rKinematicVariables);
    }


void NonLinearBernoulliBeamElement3D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points =
        r_geometry.IntegrationPoints();

    const int build_level = rCurrentProcessInfo.Has(BUILD_LEVEL)
        ? rCurrentProcessInfo[BUILD_LEVEL]
        : 0;

    const bool calculate_material_stiffness =
        build_level != 2;

    const bool calculate_geometric_stiffness =
        build_level != 1;

    // Body forces (e.g. self weight) are applied if the nodes provide a VOLUME_ACCELERATION
    const bool has_body_forces =
        r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION);

    for (
        IndexType point_number = 0;
        point_number < r_integration_points.size();
        ++point_number)
    {
        // Compute kinematics and metric
        KinematicVariables kinematic_variables(
            GetGeometry().WorkingSpaceDimension());

        CalculateKinematics(
            point_number,
            kinematic_variables);

        // Create constitutive law parameters
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(),
            GetProperties(),
            rCurrentProcessInfo);

        ConstitutiveVariables constitutive_variables(5);

        CalculateConstitutiveVariables(
            point_number,
            kinematic_variables,
            constitutive_variables,
            constitutive_law_parameters,
            ConstitutiveLaw::StressMeasure_PK2);

        Matrix b_axial;
        Matrix b_bending_1;
        Matrix b_bending_2;
        Matrix b_torsion_1;
        Matrix b_torsion_2;

        ComputeBMatrices(
            point_number,
            kinematic_variables,
            b_axial,
            b_bending_1,
            b_bending_2,
            b_torsion_1,
            b_torsion_2);

        // Assemble stiffness with cross-sectional scaling
        const double integration_weight =
            r_integration_points[point_number].Weight() *
            kinematic_variables.A;

        const double inverse_a2 =
            1.0 /
            (kinematic_variables.A *
             kinematic_variables.A);

        const double inverse_a4 =
            inverse_a2 * inverse_a2;

        if (CalculateStiffnessMatrixFlag) {
            if (calculate_material_stiffness) {
                // Material stiffness
                rLeftHandSideMatrix +=
                    inverse_a4 *
                    integration_weight *
                    prod(
                        trans(b_axial),
                        Matrix(prod(
                            constitutive_variables.ConstitutiveMatrix,
                            b_axial)));

                rLeftHandSideMatrix +=
                    inverse_a4 *
                    integration_weight *
                    prod(
                        trans(b_bending_1),
                        Matrix(prod(
                            constitutive_variables.ConstitutiveMatrix,
                            b_bending_1)));

                rLeftHandSideMatrix +=
                    inverse_a4 *
                    integration_weight *
                    prod(
                        trans(b_bending_2),
                        Matrix(prod(
                            constitutive_variables.ConstitutiveMatrix,
                            b_bending_2)));

                rLeftHandSideMatrix +=
                    inverse_a2 *
                    integration_weight *
                    prod(
                        trans(b_torsion_1),
                        Matrix(prod(
                            constitutive_variables.ConstitutiveMatrix,
                            b_torsion_1)));

                rLeftHandSideMatrix +=
                    inverse_a2 *
                    integration_weight *
                    prod(
                        trans(b_torsion_2),
                        Matrix(prod(
                            constitutive_variables.ConstitutiveMatrix,
                            b_torsion_2)));
            }

            if (calculate_geometric_stiffness) {
                // The second variations are only needed for the geometric stiffness
                Matrix g_axial;
                Matrix g_bending_1;
                Matrix g_bending_2;
                Matrix g_torsion_1;
                Matrix g_torsion_2;

                ComputeGMatrices(
                    point_number,
                    kinematic_variables,
                    g_axial,
                    g_bending_1,
                    g_bending_2,
                    g_torsion_1,
                    g_torsion_2);

                // Geometrical stiffness
                rLeftHandSideMatrix +=
                    inverse_a4 *
                    integration_weight *
                    g_axial *
                    constitutive_variables.StressVector[0];

                rLeftHandSideMatrix +=
                    inverse_a4 *
                    integration_weight *
                    g_bending_1 *
                    constitutive_variables.StressVector[1];

                rLeftHandSideMatrix +=
                    inverse_a4 *
                    integration_weight *
                    g_bending_2 *
                    constitutive_variables.StressVector[2];

                rLeftHandSideMatrix +=
                    inverse_a2 *
                    integration_weight *
                    g_torsion_1 *
                    constitutive_variables.StressVector[3];

                rLeftHandSideMatrix +=
                    inverse_a2 *
                    integration_weight *
                    g_torsion_2 *
                    constitutive_variables.StressVector[4];
            }
        }

        if (CalculateResidualVectorFlag) {
            // Internal forces. Each B matrix has a single non-zero row, the one of its generalized strain.
            const Vector& r_stress_vector = constitutive_variables.StressVector;

            noalias(rRightHandSideVector) -=
                inverse_a4 * integration_weight * r_stress_vector[0] * row(b_axial, 0);
            noalias(rRightHandSideVector) -=
                inverse_a4 * integration_weight * r_stress_vector[1] * row(b_bending_1, 1);
            noalias(rRightHandSideVector) -=
                inverse_a4 * integration_weight * r_stress_vector[2] * row(b_bending_2, 2);
            noalias(rRightHandSideVector) -=
                inverse_a2 * integration_weight * r_stress_vector[3] * row(b_torsion_1, 3);
            noalias(rRightHandSideVector) -=
                inverse_a2 * integration_weight * r_stress_vector[4] * row(b_torsion_2, 4);

            // External body forces (dead load): f_i = int N_i * rho * A * b dL, Bauer et al. (2016) Eq. (3)
            if (has_body_forces) {
                const Matrix& r_shape_functions = r_geometry.ShapeFunctionsValues();
                const double mass_per_length =
                    GetProperties()[DENSITY] * GetProperties()[CROSS_AREA];

                array_1d<double, 3> body_force = ZeroVector(3);
                for (IndexType i = 0; i < r_geometry.size(); ++i) {
                    body_force += r_shape_functions(point_number, i)
                        * r_geometry[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
                }

                for (IndexType i = 0; i < r_geometry.size(); ++i) {
                    for (IndexType j = 0; j < 3; ++j) {
                        rRightHandSideVector[i * mDofsPerNode + j] +=
                            r_shape_functions(point_number, i) * mass_per_length * body_force[j] * integration_weight;
                    }
                }
            }
        }
    }

    KRATOS_CATCH("")
}


    int NonLinearBernoulliBeamElement3D::Check(
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        BeamBaseElement3D::Check(rCurrentProcessInfo);

        const double numerical_limit =
            std::numeric_limits<double>::epsilon();

        // I_T is optional: the constitutive law uses I22 + I33 if it is not provided
        KRATOS_ERROR_IF(
            GetProperties().Has(I_T) &&
            GetProperties()[I_T] <= numerical_limit)
            << "Please provide a reasonable value for \"I_T\" "
            << "(torsional moment of inertia) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(
            !GetProperties().Has(I22) ||
            GetProperties()[I22] <= numerical_limit)
            << "Please provide a reasonable value for \"I22\" "
            << "(bending moment of inertia about N axis) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(
            !GetProperties().Has(I33) ||
            GetProperties()[I33] <= numerical_limit)
            << "Please provide a reasonable value for \"I33\" "
            << "(bending moment of inertia about V axis) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(
            !GetProperties().Has(LOCAL_AXIS_ORIENTATION))
            << "\"LOCAL_AXIS_ORIENTATION\" not provided for element #"
            << Id() << std::endl;

        // Each row is [parameter, nx, ny, nz]: the cross-section normal at that curve parameter.
        // The twist is interpolated between consecutive rows, so at least two rows with
        // strictly increasing parameters and non-zero normals are required.
        const Matrix& r_local_axis_orientation = GetProperties()[LOCAL_AXIS_ORIENTATION];

        KRATOS_ERROR_IF(r_local_axis_orientation.size2() != 4)
            << "\"LOCAL_AXIS_ORIENTATION\" must have 4 columns [parameter, nx, ny, nz], but has "
            << r_local_axis_orientation.size2() << " for element #" << Id() << std::endl;

        KRATOS_ERROR_IF(r_local_axis_orientation.size1() < 2)
            << "\"LOCAL_AXIS_ORIENTATION\" must have at least 2 rows, but has "
            << r_local_axis_orientation.size1() << " for element #" << Id() << std::endl;

        for (IndexType i = 0; i < r_local_axis_orientation.size1(); ++i) {
            KRATOS_ERROR_IF(i > 0 && r_local_axis_orientation(i, 0) <= r_local_axis_orientation(i - 1, 0))
                << "\"LOCAL_AXIS_ORIENTATION\" parameters must be strictly increasing, but row " << i
                << " has parameter " << r_local_axis_orientation(i, 0) << " after "
                << r_local_axis_orientation(i - 1, 0) << " for element #" << Id() << std::endl;

            const double normal_norm = std::sqrt(
                r_local_axis_orientation(i, 1) * r_local_axis_orientation(i, 1) +
                r_local_axis_orientation(i, 2) * r_local_axis_orientation(i, 2) +
                r_local_axis_orientation(i, 3) * r_local_axis_orientation(i, 3));

            KRATOS_ERROR_IF(normal_norm <= numerical_limit)
                << "\"LOCAL_AXIS_ORIENTATION\" row " << i
                << " has a zero normal vector for element #" << Id() << std::endl;
        }

        KRATOS_ERROR_IF(!Has(T_0) && !GetProperties().Has(T_0))
            << "\"T_0\" not provided for element #" << Id()
            << ". Run the ComputeBeamVectorsProcess first." << std::endl;

        KRATOS_ERROR_IF(!Has(N_0) && !GetProperties().Has(N_0))
            << "\"N_0\" not provided for element #" << Id()
            << ". Run the ComputeBeamVectorsProcess first." << std::endl;

        const Vector3d initial_tangent = Has(T_0)
            ? GetValue(T_0)
            : GetProperties()[T_0];

        KRATOS_ERROR_IF(norm_2(initial_tangent) <= numerical_limit)
            << "\"T_0\" reference director vector must have "
            << "non-zero magnitude for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(norm_2(GetInitialNormal()) <= numerical_limit)
            << "\"N_0\" reference director vector must have "
            << "non-zero magnitude for element #"
            << Id() << std::endl;

        return 0;

        KRATOS_CATCH("")
    }


    void NonLinearBernoulliBeamElement3D::CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        ConstitutiveVariables& rConstitutiveVariables,
        ConstitutiveLaw::Parameters& rConstitutiveLawParameters,
        const ConstitutiveLaw::StressMeasure StressMeasure)
    {
        rConstitutiveLawParameters.GetOptions().Set(
            ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN,
            true);

        rConstitutiveLawParameters.GetOptions().Set(
            ConstitutiveLaw::COMPUTE_STRESS);

        rConstitutiveLawParameters.GetOptions().Set(
            ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        rConstitutiveVariables.StrainVector[0] =
            0.5 * (
                rKinematicVariables.a * rKinematicVariables.a -
                rKinematicVariables.A * rKinematicVariables.A);

        rConstitutiveVariables.StrainVector[1] =
            rKinematicVariables.b_n -
            rKinematicVariables.B_n;

        rConstitutiveVariables.StrainVector[2] =
            rKinematicVariables.b_v -
            rKinematicVariables.B_v;

        rConstitutiveVariables.StrainVector[3] =
            rKinematicVariables.c_12 -
            rKinematicVariables.C_12;

        rConstitutiveVariables.StrainVector[4] =
            rKinematicVariables.c_13 -
            rKinematicVariables.C_13;

        rConstitutiveLawParameters.SetStrainVector(
            rConstitutiveVariables.StrainVector);

        rConstitutiveLawParameters.SetStressVector(
            rConstitutiveVariables.StressVector);

        rConstitutiveLawParameters.SetConstitutiveMatrix(
            rConstitutiveVariables.ConstitutiveMatrix);

        mConstitutiveLawVector[IntegrationPointIndex]
            ->CalculateMaterialResponse(
                rConstitutiveLawParameters,
                StressMeasure);
    }


    void NonLinearBernoulliBeamElement3D::ComputeReferenceTwistAngleAndDerivative(
        const KinematicVariables& rKinematicVariables,
        double& rPhi,
        double& rPhiDerivative)
    {

            const auto& r_geometry = GetGeometry();
            const double integration_point_coordinate = r_geometry.IntegrationPoints()[0].Coordinates()[0];

            double lower_twist_angle;
            double upper_twist_angle;
            double twist_angle_difference;
            double lower_parameter_coordinate;
            double upper_parameter_coordinate;
            int number_of_orientation_points;

            Matrix local_axis_orientation = GetProperties()[LOCAL_AXIS_ORIENTATION];
            number_of_orientation_points = local_axis_orientation.size1();

            lower_parameter_coordinate = local_axis_orientation(0, 0);
            upper_parameter_coordinate = local_axis_orientation(number_of_orientation_points - 1, 0);

            // Extract normal vector components for start and end
            Vector3d lower_normal;
            lower_normal[0] = local_axis_orientation(0, 1);
            lower_normal[1] = local_axis_orientation(0, 2);
            lower_normal[2] = local_axis_orientation(0, 3);

            Vector3d upper_normal;
            upper_normal[0] = local_axis_orientation(number_of_orientation_points - 1, 1);
            upper_normal[1] = local_axis_orientation(number_of_orientation_points - 1, 2);
            upper_normal[2] = local_axis_orientation(number_of_orientation_points - 1, 3);

            // Convert normal vectors to rotation angles (phi) for backward compatibility
            lower_twist_angle = CalculateDeltaPhi(rKinematicVariables, lower_normal);
            upper_twist_angle = CalculateDeltaPhi(rKinematicVariables, upper_normal);

            for (int i = 1; i < number_of_orientation_points; i++)
            {
                if (local_axis_orientation(i, 0) > integration_point_coordinate)
                {
                    lower_parameter_coordinate = local_axis_orientation(i - 1, 0);
                    Vector3d interpolated_normal(3);
                    interpolated_normal[0] = local_axis_orientation(i - 1, 1);
                    interpolated_normal[1] = local_axis_orientation(i - 1, 2);
                    interpolated_normal[2] = local_axis_orientation(i - 1, 3);
                    lower_twist_angle = CalculateDeltaPhi(rKinematicVariables, interpolated_normal);
                    break;
                }
            }

            for (int i = 1; i < number_of_orientation_points; i++)
            {
                if (local_axis_orientation(number_of_orientation_points - i - 1, 0) <= integration_point_coordinate)
                {
                    upper_parameter_coordinate = local_axis_orientation(number_of_orientation_points - i, 0);
                    Vector3d interpolated_normal(3);
                    interpolated_normal[0] = local_axis_orientation(number_of_orientation_points - i, 1);
                    interpolated_normal[1] = local_axis_orientation(number_of_orientation_points - i, 2);
                    interpolated_normal[2] = local_axis_orientation(number_of_orientation_points - i, 3);
                    upper_twist_angle = CalculateDeltaPhi(rKinematicVariables, interpolated_normal);
                    break;
                }
            }
            const double pi = 4.0 * std::atan(1.0);

            twist_angle_difference = (upper_twist_angle - lower_twist_angle);
            if (std::abs(upper_twist_angle - lower_twist_angle) > pi)
            {
                twist_angle_difference = twist_angle_difference - (twist_angle_difference) / std::abs(twist_angle_difference) * 2 * pi;
            }

            rPhi += lower_twist_angle + (integration_point_coordinate - lower_parameter_coordinate) / (upper_parameter_coordinate - lower_parameter_coordinate) * twist_angle_difference;
            rPhiDerivative += twist_angle_difference / (upper_parameter_coordinate - lower_parameter_coordinate);
        }

    double NonLinearBernoulliBeamElement3D::CalculateDeltaPhi(
        const KinematicVariables& rKinematicVariables,
        const Vector3d& rNormal)
    {
        const Vector3d reference_tangent =
            rKinematicVariables.R1 / norm_2(rKinematicVariables.R1);

        // Project the requested normal onto the plane normal to the reference tangent.
        Vector3d requested_axis =
            rNormal - inner_prod(reference_tangent, rNormal) * reference_tangent;
        requested_axis /= norm_2(requested_axis);

        // Untwisted cross-section axis Lambda(T0, T) * N0, i.e. Eq. (15) with Psi = 0.
        Matrix3d lambda_matrix;
        ComputeLambdaMatrix(
            GetInitialTangent(),
            reference_tangent,
            lambda_matrix);

        const Vector3d untwisted_axis =
            prod(lambda_matrix, rKinematicVariables.N0);

        // Signed angle about the reference tangent which rotates the untwisted axis onto the requested one.
        const double cos_phi = inner_prod(untwisted_axis, requested_axis);
        const double sin_phi = inner_prod(
            MathUtils<double>::CrossProduct(untwisted_axis, requested_axis),
            reference_tangent);

        double phi = std::atan2(sin_phi, cos_phi);

        if (phi < 0.0) {
            phi += 2.0 * Globals::Pi;
        }

        return phi;
    }

    NonLinearBernoulliBeamElement3D::Vector3d NonLinearBernoulliBeamElement3D::GetInitialTangent() const
    {
        const Vector3d initial_tangent = Has(T_0)
            ? GetValue(T_0)
            : GetProperties()[T_0];

        return initial_tangent / norm_2(initial_tangent);
    }

    NonLinearBernoulliBeamElement3D::Vector3d NonLinearBernoulliBeamElement3D::GetInitialNormal() const
    {
        const Vector3d initial_normal = Has(N_0)
            ? GetValue(N_0)
            : GetProperties()[N_0];

        return initial_normal;
    }

    void NonLinearBernoulliBeamElement3D::ComputeReferenceCrossSectionGeometry(
        KinematicVariables& rKinematicVariables)
    {
        const Vector3d reference_first_base_vector = rKinematicVariables.R1;
        const Vector3d reference_second_base_vector = rKinematicVariables.R2;
        const Vector3d initial_reference_tangent = GetInitialTangent();

        // Initial cross-section axes N0 and V0, orthonormal to T0. Needed by the twist angle computation below.
        rKinematicVariables.N0 = GetInitialNormal();
        rKinematicVariables.N0 -=
            inner_prod(initial_reference_tangent, rKinematicVariables.N0)
            * initial_reference_tangent;
        rKinematicVariables.N0 /= norm_2(rKinematicVariables.N0);
        rKinematicVariables.V0 = MathUtils<double>::CrossProduct(
            initial_reference_tangent,
            rKinematicVariables.N0);

        ComputeReferenceTwistAngleAndDerivative(
            rKinematicVariables,
            rKinematicVariables.Phi,
            rKinematicVariables.Phi_der);

        Matrix3d reference_lambda_matrix;
        Matrix3d reference_lambda_matrix_derivative;
        Matrix3d reference_rodrigues_matrix;
        Matrix3d reference_rodrigues_matrix_derivative;
        Matrix3d cross_section_transformation_derivative;
        cross_section_transformation_derivative.clear();

        const double reference_first_base_vector_norm =
            norm_2(reference_first_base_vector);

        const Vector3d reference_tangent =
            reference_first_base_vector / reference_first_base_vector_norm;

        const Vector3d reference_tangent_derivative =
            reference_second_base_vector / reference_first_base_vector_norm
            - inner_prod(reference_first_base_vector, reference_second_base_vector)
                / std::pow(reference_first_base_vector_norm, 3)
                * reference_first_base_vector;

        Vector3d initial_reference_tangent_derivative;
        initial_reference_tangent_derivative.clear();

        ComputeLambdaMatrix(
            initial_reference_tangent,
            reference_tangent,
            reference_lambda_matrix);

        ComputeLambdaMatrixDerivative(
            initial_reference_tangent,
            reference_tangent,
            initial_reference_tangent_derivative,
            reference_tangent_derivative,
            reference_lambda_matrix_derivative);

        ComputeRodriguesMatrix(
            reference_tangent,
            rKinematicVariables.Phi,
            reference_rodrigues_matrix);

        ComputeRodriguesMatrixDerivative(
            reference_tangent,
            reference_tangent_derivative,
            rKinematicVariables.Phi,
            rKinematicVariables.Phi_der,
            reference_rodrigues_matrix_derivative);

        rKinematicVariables.n.clear();

        Vector3d transformed_reference_normal;
        transformed_reference_normal.clear();

        noalias(transformed_reference_normal) += prod(reference_lambda_matrix, rKinematicVariables.N0);

        noalias(rKinematicVariables.n) += prod(reference_rodrigues_matrix, transformed_reference_normal);

        rKinematicVariables.n /=
            norm_2(rKinematicVariables.n);

        rKinematicVariables.v =
            MathUtils<double>::CrossProduct(reference_tangent, rKinematicVariables.n);

        noalias(cross_section_transformation_derivative) += prod(reference_rodrigues_matrix_derivative, reference_lambda_matrix);
        noalias(cross_section_transformation_derivative) += prod(reference_rodrigues_matrix, reference_lambda_matrix_derivative);

        Vector3d a21_vector;
        Vector3d a31_vector;
        a21_vector.clear();
        a31_vector.clear();

        noalias(a21_vector) += prod(cross_section_transformation_derivative, rKinematicVariables.N0);
        noalias(a31_vector) += prod(cross_section_transformation_derivative, rKinematicVariables.V0);

        rKinematicVariables.B_n =
            inner_prod(a21_vector, reference_first_base_vector);

        rKinematicVariables.B_v =
            inner_prod(a31_vector, reference_first_base_vector);

        rKinematicVariables.C_12 =
            inner_prod(a31_vector, rKinematicVariables.n);

        rKinematicVariables.C_13 =
            inner_prod(a21_vector, rKinematicVariables.v);
    }

    void NonLinearBernoulliBeamElement3D::ComputeCurrentCrossSectionGeometry(
        KinematicVariables& rKinematicVariables)
    {
        const Vector3d current_first_base_vector = rKinematicVariables.r1;
        const Vector3d current_second_base_vector = rKinematicVariables.r2;
        const Vector3d reference_first_base_vector = rKinematicVariables.R1;
        const Vector3d reference_second_base_vector = rKinematicVariables.R2;
        const Vector3d initial_reference_tangent = GetInitialTangent();

        rKinematicVariables.n.clear();
        rKinematicVariables.v.clear();

        Matrix3d current_lambda_matrix;
        Matrix3d current_lambda_matrix_derivative;
        Matrix3d reference_lambda_matrix;
        Matrix3d reference_lambda_matrix_derivative;
        Matrix3d current_rodrigues_matrix;
        Matrix3d current_rodrigues_matrix_derivative;
        Matrix3d reference_rodrigues_matrix;
        Matrix3d reference_rodrigues_matrix_derivative;

        const double current_first_base_vector_norm =
            norm_2(current_first_base_vector);

        const double reference_first_base_vector_norm =
            norm_2(reference_first_base_vector);

        const Vector3d current_tangent =
            current_first_base_vector / current_first_base_vector_norm;

        const Vector3d current_tangent_derivative =
            current_second_base_vector / current_first_base_vector_norm
            - inner_prod(current_first_base_vector, current_second_base_vector)
                / std::pow(current_first_base_vector_norm, 3)
                * current_first_base_vector;

        const Vector3d reference_tangent =
            reference_first_base_vector / reference_first_base_vector_norm;

        const Vector3d reference_tangent_derivative =
            reference_second_base_vector / reference_first_base_vector_norm
            - inner_prod(reference_first_base_vector, reference_second_base_vector)
                / std::pow(reference_first_base_vector_norm, 3)
                * reference_first_base_vector;

        Vector3d initial_reference_tangent_derivative;
        initial_reference_tangent_derivative.clear();

        ComputeLambdaMatrix(
            reference_tangent,
            current_tangent,
            current_lambda_matrix);

        ComputeLambdaMatrixDerivative(
            reference_tangent,
            current_tangent,
            reference_tangent_derivative,
            current_tangent_derivative,
            current_lambda_matrix_derivative);

        ComputeLambdaMatrix(
            initial_reference_tangent,
            reference_tangent,
            reference_lambda_matrix);

        ComputeLambdaMatrixDerivative(
            initial_reference_tangent,
            reference_tangent,
            initial_reference_tangent_derivative,
            reference_tangent_derivative,
            reference_lambda_matrix_derivative);

        ComputeRodriguesMatrix(
            current_tangent,
            rKinematicVariables.phi,
            current_rodrigues_matrix);

        ComputeRodriguesMatrixDerivative(
            current_tangent,
            current_tangent_derivative,
            rKinematicVariables.phi,
            rKinematicVariables.phi_der,
            current_rodrigues_matrix_derivative);

        ComputeRodriguesMatrix(
            reference_tangent,
            rKinematicVariables.Phi,
            reference_rodrigues_matrix);

        ComputeRodriguesMatrixDerivative(
            reference_tangent,
            reference_tangent_derivative,
            rKinematicVariables.Phi,
            rKinematicVariables.Phi_der,
            reference_rodrigues_matrix_derivative);

        Matrix3d reference_rodrigues_lambda_product;
        Matrix3d reference_rodrigues_lambda_derivative_product;
        Matrix3d reference_rodrigues_derivative_lambda_product;
        reference_rodrigues_lambda_product.clear();
        reference_rodrigues_lambda_derivative_product.clear();
        reference_rodrigues_derivative_lambda_product.clear();

        noalias(reference_rodrigues_lambda_product) += prod(reference_rodrigues_matrix, reference_lambda_matrix);
        noalias(reference_rodrigues_lambda_derivative_product) += prod(reference_rodrigues_matrix, reference_lambda_matrix_derivative);
        noalias(reference_rodrigues_derivative_lambda_product) += prod(reference_rodrigues_matrix_derivative, reference_lambda_matrix);

        Matrix3d current_lambda_reference_rodrigues_lambda_product;
        Matrix3d current_lambda_reference_rodrigues_lambda_derivative_product;
        Matrix3d current_lambda_reference_rodrigues_derivative_lambda_product;
        Matrix3d current_lambda_derivative_reference_rodrigues_lambda_product;
        current_lambda_reference_rodrigues_lambda_product.clear();
        current_lambda_reference_rodrigues_lambda_derivative_product.clear();
        current_lambda_reference_rodrigues_derivative_lambda_product.clear();
        current_lambda_derivative_reference_rodrigues_lambda_product.clear();

        noalias(current_lambda_reference_rodrigues_lambda_product) += prod(current_lambda_matrix, reference_rodrigues_lambda_product);
        noalias(current_lambda_reference_rodrigues_lambda_derivative_product) += prod(current_lambda_matrix, reference_rodrigues_lambda_derivative_product);
        noalias(current_lambda_reference_rodrigues_derivative_lambda_product) += prod(current_lambda_matrix, reference_rodrigues_derivative_lambda_product);
        noalias(current_lambda_derivative_reference_rodrigues_lambda_product) += prod(current_lambda_matrix_derivative, reference_rodrigues_lambda_product);

        Matrix3d current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term;
        Matrix3d current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term;
        Matrix3d current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term;
        Matrix3d current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;
        Matrix3d cross_section_transformation;

        current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term.clear();
        current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term.clear();
        current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term.clear();
        current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term.clear();
        cross_section_transformation.clear();

        noalias(current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_lambda_derivative_product);
        noalias(current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_derivative_lambda_product);
        noalias(current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term) += prod(current_rodrigues_matrix, current_lambda_derivative_reference_rodrigues_lambda_product);
        noalias(current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term) += prod(current_rodrigues_matrix_derivative, current_lambda_reference_rodrigues_lambda_product);
        noalias(cross_section_transformation) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_lambda_product);

        const Matrix3d cross_section_transformation_derivative =
            current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term
            + current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term
            + current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term
            + current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;

        Vector3d a21_vector;
        Vector3d a31_vector;
        a21_vector.clear();
        a31_vector.clear();

        noalias(a21_vector) += prod(cross_section_transformation_derivative, rKinematicVariables.N0);
        noalias(a31_vector) += prod(cross_section_transformation_derivative, rKinematicVariables.V0);
        noalias(rKinematicVariables.n) += prod(cross_section_transformation, rKinematicVariables.N0);
        noalias(rKinematicVariables.v) += prod(cross_section_transformation, rKinematicVariables.V0);

        rKinematicVariables.b_n =
            inner_prod(a21_vector, current_first_base_vector);

        rKinematicVariables.b_v =
            inner_prod(a31_vector, current_first_base_vector);

        rKinematicVariables.c_12 =
            inner_prod(a31_vector, rKinematicVariables.n);

        rKinematicVariables.c_13 =
            inner_prod(a21_vector, rKinematicVariables.v);
    }


    void NonLinearBernoulliBeamElement3D::ComputeBMatrices(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        Matrix& rBAxial,
        Matrix& rBBending1,
        Matrix& rBBending2,
        Matrix& rBTorsion1,
        Matrix& rBTorsion2)
    {
        KRATOS_TRY

        // Scratch matrices, laid out as (3 * number of dofs) x 3 or (3 * number of dofs) x (3 * number of dofs)
        Matrix lambda_derivative_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_derivative_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_variation_rodrigues_derivative_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_variation_rodrigues_lambda_derivative_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_derivative_lambda_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_derivative_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_derivative_variation_lambda_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_derivative_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_rodrigues_lambda_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_variation_rodrigues_derivative_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_variation_rodrigues_lambda_derivative_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_lambda_derivative_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_lambda_rodrigues_derivative_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_lambda_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_lambda_rodrigues_lambda_derivative_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);

        const auto& r_geometry = GetGeometry();

        // Get shape functions and derivatives at integration point
        Vector shape_functions = row(r_geometry.ShapeFunctionsValues(GetIntegrationMethod()), IntegrationPointIndex);
        Vector shape_function_derivatives = column(r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex, GetIntegrationMethod()), 0);
        Vector shape_function_second_derivatives = column(r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetIntegrationMethod()), 0);

        // Extract kinematic variables
        Vector3d current_first_base_vector = rKinematicVariables.r1;
        Vector3d current_second_base_vector = rKinematicVariables.r2;
        Vector3d reference_first_base_vector = rKinematicVariables.R1;
        Vector3d reference_second_base_vector = rKinematicVariables.R2;
        Vector3d reference_normal = rKinematicVariables.N0;
        Vector3d reference_binormal = rKinematicVariables.V0;
        double current_twist_angle = rKinematicVariables.phi;
        double current_twist_angle_derivative = rKinematicVariables.phi_der;
        double reference_twist_angle = rKinematicVariables.Phi;
        double reference_twist_angle_derivative = rKinematicVariables.Phi_der;
        Vector3d initial_reference_tangent = GetInitialTangent();

        // Compute curvature and torsion variations
        Vector axial_variation(mNumberOfDofs);
        Vector normal_curvature_variation(mNumberOfDofs);
        Vector binormal_curvature_variation(mNumberOfDofs);
        Vector normal_torsion_variation(mNumberOfDofs);
        Vector binormal_torsion_variation(mNumberOfDofs);
        axial_variation.clear();
        normal_curvature_variation.clear();
        binormal_curvature_variation.clear();
        normal_torsion_variation.clear();
        binormal_torsion_variation.clear();

        //compute axial_variation
        for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
        {
            size_t dof_component = dof_index % mDofsPerNode;
            size_t control_point_index = dof_index / mDofsPerNode;
            if (dof_component > 2)
                axial_variation[dof_index] = 0;
            else
                axial_variation[dof_index] = rKinematicVariables.r1[dof_component] * shape_function_derivatives[control_point_index];
        }
        //compute normal_curvature_variation, binormal_curvature_variation, normal_torsion_variation, binormal_torsion_variation
        Vector3d current_tangent;
        current_tangent.clear();
        Vector3d current_tangent_derivative;
        current_tangent_derivative.clear();
        Vector3d reference_tangent;
        reference_tangent.clear();
        Vector3d reference_tangent_derivative;
        reference_tangent_derivative.clear();
        Vector3d initial_reference_tangent_derivative;
        initial_reference_tangent_derivative.clear();
        Vector current_tangent_variation;
        Vector current_tangent_derivative_variation;
        Matrix3d current_lambda_matrix;
        Matrix3d current_lambda_matrix_derivative;
        Matrix3d reference_lambda_matrix;
        Matrix3d reference_lambda_matrix_derivative;
        Matrix3d current_rodrigues_matrix;
        Matrix3d current_rodrigues_matrix_derivative;
        Matrix3d reference_rodrigues_matrix;
        Matrix3d reference_rodrigues_matrix_derivative;

        current_tangent = current_first_base_vector / norm_2(current_first_base_vector);
        current_tangent_derivative = current_second_base_vector / norm_2(current_first_base_vector) - inner_prod(current_first_base_vector, current_second_base_vector) / std::pow(norm_2(current_first_base_vector), 3) * current_first_base_vector;
        reference_tangent = reference_first_base_vector / norm_2(reference_first_base_vector);
        reference_tangent_derivative = reference_second_base_vector / norm_2(reference_first_base_vector) - inner_prod(reference_first_base_vector, reference_second_base_vector) / std::pow(norm_2(reference_first_base_vector), 3) * reference_first_base_vector;
        ComputeTangentVariation(shape_function_derivatives, current_first_base_vector, current_tangent_variation);
        ComputeTangentDerivativeVariation(shape_function_derivatives, shape_function_second_derivatives, current_first_base_vector, current_second_base_vector, current_tangent_derivative_variation);

        ComputeLambdaMatrix(reference_tangent, current_tangent, current_lambda_matrix);
        ComputeLambdaMatrixDerivative(reference_tangent, current_tangent, reference_tangent_derivative, current_tangent_derivative, current_lambda_matrix_derivative);
        ComputeLambdaMatrixVariation(reference_tangent, current_tangent, current_tangent_variation, mNumberOfDofs, mDofsPerNode, lambda_variation_matrix);
        ComputeLambdaMatrixDerivativeVariation(reference_tangent, current_tangent, reference_tangent_derivative, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, mNumberOfDofs, mDofsPerNode, lambda_derivative_variation_matrix);
        ComputeLambdaMatrix(initial_reference_tangent, reference_tangent, reference_lambda_matrix);
        ComputeLambdaMatrixDerivative(initial_reference_tangent, reference_tangent, initial_reference_tangent_derivative, reference_tangent_derivative, reference_lambda_matrix_derivative);

        ComputeRodriguesMatrix(current_tangent, current_twist_angle, current_rodrigues_matrix);
        ComputeRodriguesMatrixDerivative(current_tangent, current_tangent_derivative, current_twist_angle, current_twist_angle_derivative, current_rodrigues_matrix_derivative);
        ComputeRodriguesMatrixVariation(current_tangent, current_tangent_variation, shape_functions, current_twist_angle, mNumberOfDofs, mDofsPerNode, rodrigues_variation_matrix);
        ComputeRodriguesMatrixDerivativeVariation(current_tangent, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, shape_functions, shape_function_derivatives, current_twist_angle, current_twist_angle_derivative, mNumberOfDofs, mDofsPerNode, rodrigues_derivative_variation_matrix);
        ComputeRodriguesMatrix(reference_tangent, reference_twist_angle, reference_rodrigues_matrix);
        ComputeRodriguesMatrixDerivative(reference_tangent, reference_tangent_derivative, reference_twist_angle, reference_twist_angle_derivative, reference_rodrigues_matrix_derivative);

        Matrix3d reference_rodrigues_lambda_derivative_product;
        reference_rodrigues_lambda_derivative_product.clear();
        Matrix3d reference_rodrigues_derivative_lambda_product;
        reference_rodrigues_derivative_lambda_product.clear();
        Matrix3d reference_rodrigues_lambda_product;
        reference_rodrigues_lambda_product.clear();
        Matrix3d reference_rodrigues_lambda_product_derivative;
        reference_rodrigues_lambda_product_derivative.clear();

        noalias(reference_rodrigues_lambda_product) += prod(reference_rodrigues_matrix, reference_lambda_matrix);
        noalias(reference_rodrigues_lambda_derivative_product) += prod(reference_rodrigues_matrix, reference_lambda_matrix_derivative);
        noalias(reference_rodrigues_derivative_lambda_product) += prod(reference_rodrigues_matrix_derivative, reference_lambda_matrix);
        reference_rodrigues_lambda_product_derivative = reference_rodrigues_lambda_derivative_product + reference_rodrigues_derivative_lambda_product;

        Matrix3d current_lambda_reference_rodrigues_lambda_derivative_product;
        current_lambda_reference_rodrigues_lambda_derivative_product.clear();
        Matrix3d current_lambda_reference_rodrigues_derivative_lambda_product;
        current_lambda_reference_rodrigues_derivative_lambda_product.clear();
        Matrix3d current_lambda_reference_rodrigues_lambda_product;
        current_lambda_reference_rodrigues_lambda_product.clear();
        Matrix3d current_lambda_derivative_reference_rodrigues_lambda_product;
        current_lambda_derivative_reference_rodrigues_lambda_product.clear();
        Matrix3d current_lambda_reference_rodrigues_lambda_product_derivative;
        current_lambda_reference_rodrigues_lambda_product_derivative.clear();

        noalias(current_lambda_reference_rodrigues_lambda_product) += prod(current_lambda_matrix, reference_rodrigues_lambda_product);
        noalias(current_lambda_reference_rodrigues_lambda_derivative_product) += prod(current_lambda_matrix, reference_rodrigues_lambda_derivative_product);
        noalias(current_lambda_reference_rodrigues_derivative_lambda_product) += prod(current_lambda_matrix, reference_rodrigues_derivative_lambda_product);
        noalias(current_lambda_derivative_reference_rodrigues_lambda_product) += prod(current_lambda_matrix_derivative, reference_rodrigues_lambda_product);

        current_lambda_reference_rodrigues_lambda_product_derivative = current_lambda_reference_rodrigues_lambda_derivative_product + current_lambda_reference_rodrigues_derivative_lambda_product + current_lambda_derivative_reference_rodrigues_lambda_product;

        Matrix3d current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term;
        current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term.clear();
        Matrix3d current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term;
        current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term.clear();
        Matrix3d current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term;
        current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term.clear();
        Matrix3d current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;
        current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term.clear();
        Matrix3d cross_section_transformation;
        cross_section_transformation.clear();

        noalias(current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_lambda_derivative_product);
        noalias(current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_derivative_lambda_product);
        noalias(current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term) += prod(current_rodrigues_matrix, current_lambda_derivative_reference_rodrigues_lambda_product);
        noalias(current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term) += prod(current_rodrigues_matrix_derivative, current_lambda_reference_rodrigues_lambda_product);
        noalias(cross_section_transformation) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_lambda_product);

        Matrix3d cross_section_transformation_derivative;
        cross_section_transformation_derivative.clear();
        cross_section_transformation_derivative = current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term + current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term + current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term + current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;

        lambda_variation_rodrigues_lambda_derivative_matrix.clear();
        lambda_variation_rodrigues_derivative_lambda_matrix.clear();
        lambda_variation_rodrigues_lambda_matrix.clear();
        lambda_derivative_variation_rodrigues_lambda_matrix.clear();

        noalias(lambda_variation_rodrigues_lambda_matrix) += prod(lambda_variation_matrix, reference_rodrigues_lambda_product);
        noalias(lambda_variation_rodrigues_lambda_derivative_matrix) += prod(lambda_variation_matrix, reference_rodrigues_lambda_derivative_product);
        noalias(lambda_variation_rodrigues_derivative_lambda_matrix) += prod(lambda_variation_matrix, reference_rodrigues_derivative_lambda_product);
        noalias(lambda_derivative_variation_rodrigues_lambda_matrix) += prod(lambda_derivative_variation_matrix, reference_rodrigues_lambda_product);

        Matrix current_lambda_reference_rodrigues_lambda_product_derivative_variation;
        current_lambda_reference_rodrigues_lambda_product_derivative_variation.resize(3 * mNumberOfDofs, 3);
        current_lambda_reference_rodrigues_lambda_product_derivative_variation.clear();

        current_lambda_reference_rodrigues_lambda_product_derivative_variation = lambda_variation_rodrigues_lambda_derivative_matrix + lambda_variation_rodrigues_derivative_lambda_matrix + lambda_derivative_variation_rodrigues_lambda_matrix;

        rodrigues_variation_lambda_rodrigues_lambda_derivative_matrix.clear();
        rodrigues_variation_lambda_rodrigues_derivative_lambda_matrix.clear();
        rodrigues_variation_lambda_derivative_rodrigues_lambda_matrix.clear();
        rodrigues_derivative_variation_lambda_rodrigues_lambda_matrix.clear();
        rodrigues_derivative_lambda_variation_rodrigues_lambda_matrix.clear();
        rodrigues_lambda_derivative_variation_rodrigues_lambda_matrix.clear();
        rodrigues_lambda_variation_rodrigues_derivative_lambda_matrix.clear();
        rodrigues_lambda_variation_rodrigues_lambda_derivative_matrix.clear();
        rodrigues_lambda_variation_rodrigues_lambda_matrix.clear();
        rodrigues_variation_lambda_rodrigues_lambda_matrix.clear();

        noalias(rodrigues_variation_lambda_rodrigues_lambda_derivative_matrix) += prod(rodrigues_variation_matrix, current_lambda_reference_rodrigues_lambda_derivative_product);
        noalias(rodrigues_variation_lambda_rodrigues_derivative_lambda_matrix) += prod(rodrigues_variation_matrix, current_lambda_reference_rodrigues_derivative_lambda_product);
        noalias(rodrigues_variation_lambda_derivative_rodrigues_lambda_matrix) += prod(rodrigues_variation_matrix, current_lambda_derivative_reference_rodrigues_lambda_product);
        noalias(rodrigues_derivative_variation_lambda_rodrigues_lambda_matrix) += prod(rodrigues_derivative_variation_matrix, current_lambda_reference_rodrigues_lambda_product);
        AddMatrixTimesVariation(rodrigues_lambda_variation_rodrigues_lambda_derivative_matrix, current_rodrigues_matrix, lambda_variation_rodrigues_lambda_derivative_matrix, mNumberOfDofs);
        AddMatrixTimesVariation(rodrigues_lambda_variation_rodrigues_derivative_lambda_matrix, current_rodrigues_matrix, lambda_variation_rodrigues_derivative_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesVariation(rodrigues_lambda_derivative_variation_rodrigues_lambda_matrix, current_rodrigues_matrix, lambda_derivative_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesVariation(rodrigues_derivative_lambda_variation_rodrigues_lambda_matrix, current_rodrigues_matrix_derivative, lambda_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesVariation(rodrigues_lambda_variation_rodrigues_lambda_matrix, current_rodrigues_matrix, lambda_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        noalias(rodrigues_variation_lambda_rodrigues_lambda_matrix) += prod(rodrigues_variation_matrix, current_lambda_reference_rodrigues_lambda_product);

        rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix.clear();
        rodrigues_lambda_rodrigues_lambda_variation_matrix.clear();

        rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix = rodrigues_variation_lambda_rodrigues_lambda_derivative_matrix + rodrigues_variation_lambda_rodrigues_derivative_lambda_matrix + rodrigues_variation_lambda_derivative_rodrigues_lambda_matrix + rodrigues_derivative_variation_lambda_rodrigues_lambda_matrix
            + rodrigues_lambda_variation_rodrigues_lambda_derivative_matrix + rodrigues_lambda_variation_rodrigues_derivative_lambda_matrix + rodrigues_lambda_derivative_variation_rodrigues_lambda_matrix + rodrigues_derivative_lambda_variation_rodrigues_lambda_matrix;
        rodrigues_lambda_rodrigues_lambda_variation_matrix = rodrigues_variation_lambda_rodrigues_lambda_matrix + rodrigues_lambda_variation_rodrigues_lambda_matrix;

        Vector3d normal_vector;
        normal_vector.clear();
        Vector3d binormal_vector;
        binormal_vector.clear();

        noalias(normal_vector) += prod(cross_section_transformation, reference_normal);
        noalias(binormal_vector) += prod(cross_section_transformation, reference_binormal);

        Vector normal_vector_variation;
        normal_vector_variation.resize(3 * mNumberOfDofs);
        normal_vector_variation.clear();
        Vector binormal_vector_variation;
        binormal_vector_variation.resize(3 * mNumberOfDofs);
        binormal_vector_variation.clear();

        noalias(normal_vector_variation) += prod(rodrigues_lambda_rodrigues_lambda_variation_matrix, reference_normal);
        noalias(binormal_vector_variation) += prod(rodrigues_lambda_rodrigues_lambda_variation_matrix, reference_binormal);

        Vector current_first_base_vector_variation;
        current_first_base_vector_variation.resize(3 * mNumberOfDofs);
        current_first_base_vector_variation.clear();
        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (size_t dof_index = 0; dof_index < mNumberOfDofs; dof_index++)
            {
                size_t dof_component = dof_index % mDofsPerNode;
                size_t control_point_index = dof_index / mDofsPerNode;
                if (row_component == dof_component)
                    current_first_base_vector_variation(row_component * mNumberOfDofs + dof_index) = shape_function_derivatives[control_point_index];
            }
        }

        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (int product_index = 0; product_index < 3; product_index++)
            {
                for (size_t dof_index = 0; dof_index < mNumberOfDofs; dof_index++)
                {
                    normal_curvature_variation(dof_index) += rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index) * current_first_base_vector[row_component] + cross_section_transformation_derivative(row_component, product_index) * reference_normal(product_index) * current_first_base_vector_variation[row_component * mNumberOfDofs + dof_index];
                    binormal_curvature_variation(dof_index) += rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index) * current_first_base_vector[row_component] + cross_section_transformation_derivative(row_component, product_index) * reference_binormal(product_index) * current_first_base_vector_variation[row_component * mNumberOfDofs + dof_index];
                    normal_torsion_variation(dof_index) += rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index) * normal_vector(row_component) + cross_section_transformation_derivative(row_component, product_index) * reference_binormal(product_index) * normal_vector_variation(row_component * mNumberOfDofs + dof_index);
                    binormal_torsion_variation(dof_index) += rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index) * binormal_vector(row_component) + cross_section_transformation_derivative(row_component, product_index) * reference_normal(product_index) * binormal_vector_variation(row_component * mNumberOfDofs + dof_index);
                }
            }
        }


        // Initialize B matrices
        rBAxial.resize(5, mNumberOfDofs);
        rBBending1.resize(5, mNumberOfDofs);
        rBBending2.resize(5, mNumberOfDofs);
        rBTorsion1.resize(5, mNumberOfDofs);
        rBTorsion2.resize(5, mNumberOfDofs);
        rBAxial.clear();
        rBBending1.clear();
        rBBending2.clear();
        rBTorsion1.clear();
        rBTorsion2.clear();

        noalias(row(rBAxial, 0)) = axial_variation; //Normal force
        noalias(row(rBBending1, 1)) = normal_curvature_variation;  // Bending about n-axis
        noalias(row(rBBending2, 2)) = binormal_curvature_variation;  // Bending about v-axis
        noalias(row(rBTorsion1, 3)) = normal_torsion_variation;   // Torsion n in row 1 (shear stress)
        noalias(row(rBTorsion2, 4)) = binormal_torsion_variation;   // Torsion v in row 2 (shear stress)

        KRATOS_CATCH("")
    }


    /**
     * Second variations of the generalized strains with respect to the dofs (Bauer et al. 2016, Eq. 70 and Appendix A).
     *
     * The current cross-section axes are a_alpha = T * A0_alpha with the frame transformation
     * T = R_t(psi) * Lambda(T, t) * P, where P = R_T(Psi) * Lambda(T0, T) only depends on the reference
     * configuration. Their derivative along the curve is a_alpha,1 = T,1 * A0_alpha. The second variations of
     * T and T,1 follow from the product rule and are assembled term by term below.
     *
     * Layout of the variation containers (n = number of element dofs):
     *  - 3x3 matrix M, first variation:  (3 n) x 3     with M,r(i, j)  stored at (i * n + r, j)
     *  - 3x3 matrix M, second variation: (3 n) x (3 n) with M,rs(i, j) stored at (i * n + r, j * n + s)
     *  - 3-vector v, first variation:    3 n            with v,r(i)     stored at (i * n + r)
     */
    void NonLinearBernoulliBeamElement3D::ComputeGMatrices(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        Matrix& rGAxial,
        Matrix& rGBending1,
        Matrix& rGBending2,
        Matrix& rGTorsion1,
        Matrix& rGTorsion2)
    {
        KRATOS_TRY

        // Scratch matrices, laid out as (3 * number of dofs) x 3 or (3 * number of dofs) x (3 * number of dofs)
        Matrix lambda_derivative_second_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        Matrix lambda_derivative_second_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        Matrix lambda_derivative_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_derivative_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_second_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        Matrix lambda_second_variation_rodrigues_derivative_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        Matrix lambda_second_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        Matrix lambda_second_variation_rodrigues_lambda_derivative_matrix = ZeroMatrix(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        Matrix lambda_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_variation_rodrigues_derivative_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix lambda_variation_rodrigues_lambda_derivative_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_derivative_lambda_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_derivative_second_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        Matrix rodrigues_derivative_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_derivative_variation_lambda_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_derivative_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_rodrigues_lambda_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_variation_rodrigues_derivative_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_variation_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_lambda_variation_rodrigues_lambda_derivative_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_second_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        Matrix rodrigues_variation_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_lambda_derivative_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_lambda_rodrigues_derivative_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_lambda_rodrigues_lambda_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);
        Matrix rodrigues_variation_lambda_rodrigues_lambda_derivative_matrix = ZeroMatrix(3 * mNumberOfDofs, 3);


        const auto& r_geometry = GetGeometry();

        // Get shape functions and derivatives at integration point
        Vector shape_functions = row(r_geometry.ShapeFunctionsValues(GetIntegrationMethod()), IntegrationPointIndex);
        Vector shape_function_derivatives = column(r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex, GetIntegrationMethod()), 0);
        Vector shape_function_second_derivatives = column(r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetIntegrationMethod()), 0);

        // Extract kinematic variables
        Vector3d current_first_base_vector = rKinematicVariables.r1;
        Vector3d current_second_base_vector = rKinematicVariables.r2;
        Vector3d reference_first_base_vector = rKinematicVariables.R1;
        Vector3d reference_second_base_vector = rKinematicVariables.R2;
        Vector3d reference_normal = rKinematicVariables.N0;
        Vector3d reference_binormal = rKinematicVariables.V0;
        double current_twist_angle = rKinematicVariables.phi;
        double current_twist_angle_derivative = rKinematicVariables.phi_der;
        double reference_twist_angle = rKinematicVariables.Phi;
        double reference_twist_angle_derivative = rKinematicVariables.Phi_der;
        Vector3d initial_reference_tangent = GetInitialTangent();

        // Compute curvature and torsion variations
        Matrix axial_second_variation(mNumberOfDofs, mNumberOfDofs);
        Matrix normal_curvature_second_variation(mNumberOfDofs, mNumberOfDofs);
        Matrix binormal_curvature_second_variation(mNumberOfDofs, mNumberOfDofs);
        Matrix normal_torsion_second_variation(mNumberOfDofs, mNumberOfDofs);
        Matrix binormal_torsion_second_variation(mNumberOfDofs, mNumberOfDofs);
        axial_second_variation.clear();
        normal_curvature_second_variation.clear();
        binormal_curvature_second_variation.clear();
        normal_torsion_second_variation.clear();
        binormal_torsion_second_variation.clear();

        //compute axial_variation
        for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
        {
            size_t dof_component = dof_index % mDofsPerNode;
            size_t control_point_index = dof_index / mDofsPerNode;
            if (dof_component > 2)
                for (size_t  second_dof_index  = 0;second_dof_index < mNumberOfDofs;second_dof_index++)
                    axial_second_variation(dof_index, second_dof_index) = 0.0;
            else
            {
                for (size_t  second_dof_index  = 0;second_dof_index < mNumberOfDofs;second_dof_index++)
                {
                    size_t second_dof_component = second_dof_index % mDofsPerNode;
                    int second_control_point_index = second_dof_index / mDofsPerNode;
                    if (second_dof_component > 2)
                        axial_second_variation(dof_index, second_dof_index) = 0;
                    else
                        if (dof_component == second_dof_component)
                            axial_second_variation(dof_index, second_dof_index) = shape_function_derivatives[control_point_index] * shape_function_derivatives[second_control_point_index];
                        else
                            axial_second_variation(dof_index, second_dof_index) = 0;
                }
            }
        }

        //compute normal_curvature_second_variation, binormal_curvature_second_variation, normal_torsion_second_variation, binormal_torsion_second_variation
        Vector3d current_tangent;
        current_tangent.clear();
        Vector3d current_tangent_derivative;
        current_tangent_derivative.clear();
        Vector3d reference_tangent;
        reference_tangent.clear();
        Vector3d reference_tangent_derivative;
        reference_tangent_derivative.clear();
        Vector3d initial_reference_tangent_derivative;
        initial_reference_tangent_derivative.clear();
        Vector current_tangent_variation;
        Vector current_tangent_derivative_variation;
        Matrix current_tangent_second_variation;
        Matrix current_tangent_derivative_second_variation;
        Matrix3d current_lambda_matrix;
        Matrix3d current_lambda_matrix_derivative;
        Matrix3d reference_lambda_matrix;
        Matrix3d reference_lambda_matrix_derivative;
        Matrix3d current_rodrigues_matrix;
        Matrix3d current_rodrigues_matrix_derivative;
        Matrix3d reference_rodrigues_matrix;
        Matrix3d reference_rodrigues_matrix_derivative;

        current_tangent = current_first_base_vector / norm_2(current_first_base_vector);
        current_tangent_derivative = current_second_base_vector / norm_2(current_first_base_vector) - inner_prod(current_first_base_vector, current_second_base_vector) / std::pow(norm_2(current_first_base_vector), 3) * current_first_base_vector;
        reference_tangent = reference_first_base_vector / norm_2(reference_first_base_vector);
        reference_tangent_derivative = reference_second_base_vector / norm_2(reference_first_base_vector) - inner_prod(reference_first_base_vector, reference_second_base_vector) / std::pow(norm_2(reference_first_base_vector), 3) * reference_first_base_vector;
        // 1. Current tangent t = a1 / |a1|, its derivative t,1 and their first and second variations
        ComputeTangentVariation(shape_function_derivatives, current_first_base_vector, current_tangent_variation);
        ComputeTangentDerivativeVariation(shape_function_derivatives, shape_function_second_derivatives, current_first_base_vector, current_second_base_vector, current_tangent_derivative_variation);
        ComputeTangentSecondVariation(shape_function_derivatives, current_first_base_vector, current_tangent_second_variation);
        ComputeTangentDerivativeSecondVariation(shape_function_derivatives, shape_function_second_derivatives, current_first_base_vector, current_second_base_vector, current_tangent_derivative_second_variation);

        // 2. Mapping matrices (Eq. 23): Lambda(T, t) depends on the dofs (Eq. A.13-A.19), Lambda(T0, T) does not
        ComputeLambdaMatrix(reference_tangent, current_tangent, current_lambda_matrix);
        ComputeLambdaMatrixDerivative(reference_tangent, current_tangent, reference_tangent_derivative, current_tangent_derivative, current_lambda_matrix_derivative);

        ComputeLambdaMatrix(initial_reference_tangent, reference_tangent, reference_lambda_matrix);
        ComputeLambdaMatrixDerivative(initial_reference_tangent, reference_tangent, initial_reference_tangent_derivative, reference_tangent_derivative, reference_lambda_matrix_derivative);
        ComputeLambdaMatrixVariations(reference_tangent, current_tangent, reference_tangent_derivative, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, current_tangent_second_variation, current_tangent_derivative_second_variation, mNumberOfDofs, mDofsPerNode, lambda_variation_matrix, lambda_derivative_variation_matrix, lambda_second_variation_matrix, lambda_derivative_second_variation_matrix);

        // 3. Rotation matrices (Eq. 26): R_t(psi) depends on the dofs (Eq. A.5-A.11), R_T(Psi) does not
        ComputeRodriguesMatrix(current_tangent, current_twist_angle, current_rodrigues_matrix);
        ComputeRodriguesMatrixDerivative(current_tangent, current_tangent_derivative, current_twist_angle, current_twist_angle_derivative, current_rodrigues_matrix_derivative);
        ComputeRodriguesMatrixVariation(current_tangent, current_tangent_variation, shape_functions, current_twist_angle, mNumberOfDofs, mDofsPerNode, rodrigues_variation_matrix);
        ComputeRodriguesMatrixDerivativeVariation(current_tangent, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, shape_functions, shape_function_derivatives, current_twist_angle, current_twist_angle_derivative, mNumberOfDofs, mDofsPerNode, rodrigues_derivative_variation_matrix);
        ComputeRodriguesMatrixSecondVariation(current_tangent, current_tangent_variation, current_tangent_second_variation, shape_functions, current_twist_angle, mNumberOfDofs, mDofsPerNode, rodrigues_second_variation_matrix);
        ComputeRodriguesMatrixDerivativeSecondVariation(current_tangent, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, current_tangent_second_variation, current_tangent_derivative_second_variation, shape_functions, shape_function_derivatives, current_twist_angle, current_twist_angle_derivative, mNumberOfDofs, mDofsPerNode, rodrigues_derivative_second_variation_matrix);

        ComputeRodriguesMatrix(reference_tangent, reference_twist_angle, reference_rodrigues_matrix);
        ComputeRodriguesMatrixDerivative(reference_tangent, reference_tangent_derivative, reference_twist_angle, reference_twist_angle_derivative, reference_rodrigues_matrix_derivative);

        // 4. Products of R_t, Lambda(T, t) and P and of their derivatives and variations

        Matrix3d reference_rodrigues_lambda_derivative_product;
        reference_rodrigues_lambda_derivative_product.clear();
        Matrix3d reference_rodrigues_derivative_lambda_product;
        reference_rodrigues_derivative_lambda_product.clear();
        Matrix3d reference_rodrigues_lambda_product;
        reference_rodrigues_lambda_product.clear();
        Matrix3d reference_rodrigues_lambda_product_derivative;
        reference_rodrigues_lambda_product_derivative.clear();

        noalias(reference_rodrigues_lambda_product) += prod(reference_rodrigues_matrix, reference_lambda_matrix);
        noalias(reference_rodrigues_lambda_derivative_product) += prod(reference_rodrigues_matrix, reference_lambda_matrix_derivative);
        noalias(reference_rodrigues_derivative_lambda_product) += prod(reference_rodrigues_matrix_derivative, reference_lambda_matrix);
        reference_rodrigues_lambda_product_derivative = reference_rodrigues_lambda_derivative_product + reference_rodrigues_derivative_lambda_product;

        Matrix3d current_lambda_reference_rodrigues_lambda_derivative_product;
        current_lambda_reference_rodrigues_lambda_derivative_product.clear();
        Matrix3d current_lambda_reference_rodrigues_derivative_lambda_product;
        current_lambda_reference_rodrigues_derivative_lambda_product.clear();
        Matrix3d current_lambda_reference_rodrigues_lambda_product;
        current_lambda_reference_rodrigues_lambda_product.clear();
        Matrix3d current_lambda_derivative_reference_rodrigues_lambda_product;
        current_lambda_derivative_reference_rodrigues_lambda_product.clear();
        Matrix3d current_lambda_reference_rodrigues_lambda_product_derivative;
        current_lambda_reference_rodrigues_lambda_product_derivative.clear();

        noalias(current_lambda_reference_rodrigues_lambda_product) += prod(current_lambda_matrix, reference_rodrigues_lambda_product);
        noalias(current_lambda_reference_rodrigues_lambda_derivative_product) += prod(current_lambda_matrix, reference_rodrigues_lambda_derivative_product);
        noalias(current_lambda_reference_rodrigues_derivative_lambda_product) += prod(current_lambda_matrix, reference_rodrigues_derivative_lambda_product);
        noalias(current_lambda_derivative_reference_rodrigues_lambda_product) += prod(current_lambda_matrix_derivative, reference_rodrigues_lambda_product);

        current_lambda_reference_rodrigues_lambda_product_derivative = current_lambda_reference_rodrigues_lambda_derivative_product + current_lambda_reference_rodrigues_derivative_lambda_product + current_lambda_derivative_reference_rodrigues_lambda_product;

        Matrix3d current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term;
        current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term.clear();
        Matrix3d current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term;
        current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term.clear();
        Matrix3d current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term;
        current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term.clear();
        Matrix3d current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;
        current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term.clear();
        Matrix3d cross_section_transformation;
        cross_section_transformation.clear();

        noalias(current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_lambda_derivative_product);
        noalias(current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_derivative_lambda_product);
        noalias(current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term) += prod(current_rodrigues_matrix, current_lambda_derivative_reference_rodrigues_lambda_product);
        noalias(current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term) += prod(current_rodrigues_matrix_derivative, current_lambda_reference_rodrigues_lambda_product);
        noalias(cross_section_transformation) += prod(current_rodrigues_matrix, current_lambda_reference_rodrigues_lambda_product);

        Matrix3d cross_section_transformation_derivative;
        cross_section_transformation_derivative.clear();
        cross_section_transformation_derivative = current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term + current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term + current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term + current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;


        lambda_variation_rodrigues_lambda_derivative_matrix.clear();
        lambda_variation_rodrigues_derivative_lambda_matrix.clear();
        lambda_variation_rodrigues_lambda_matrix.clear();
        lambda_derivative_variation_rodrigues_lambda_matrix.clear();

        noalias(lambda_variation_rodrigues_lambda_matrix) += prod(lambda_variation_matrix, reference_rodrigues_lambda_product);
        noalias(lambda_variation_rodrigues_lambda_derivative_matrix) += prod(lambda_variation_matrix, reference_rodrigues_lambda_derivative_product);
        noalias(lambda_variation_rodrigues_derivative_lambda_matrix) += prod(lambda_variation_matrix, reference_rodrigues_derivative_lambda_product);
        noalias(lambda_derivative_variation_rodrigues_lambda_matrix) += prod(lambda_derivative_variation_matrix, reference_rodrigues_lambda_product);

        Matrix current_lambda_reference_rodrigues_lambda_product_derivative_variation;
        current_lambda_reference_rodrigues_lambda_product_derivative_variation.resize(3 * mNumberOfDofs, 3);
        current_lambda_reference_rodrigues_lambda_product_derivative_variation.clear();

        current_lambda_reference_rodrigues_lambda_product_derivative_variation = lambda_variation_rodrigues_lambda_derivative_matrix + lambda_variation_rodrigues_derivative_lambda_matrix + lambda_derivative_variation_rodrigues_lambda_matrix;

        rodrigues_variation_lambda_rodrigues_lambda_derivative_matrix.clear();
        rodrigues_variation_lambda_rodrigues_derivative_lambda_matrix.clear();
        rodrigues_variation_lambda_derivative_rodrigues_lambda_matrix.clear();
        rodrigues_derivative_variation_lambda_rodrigues_lambda_matrix.clear();
        rodrigues_derivative_lambda_variation_rodrigues_lambda_matrix.clear();
        rodrigues_lambda_derivative_variation_rodrigues_lambda_matrix.clear();
        rodrigues_lambda_variation_rodrigues_derivative_lambda_matrix.clear();
        rodrigues_lambda_variation_rodrigues_lambda_derivative_matrix.clear();
        rodrigues_lambda_variation_rodrigues_lambda_matrix.clear();
        rodrigues_variation_lambda_rodrigues_lambda_matrix.clear();

        noalias(rodrigues_variation_lambda_rodrigues_lambda_derivative_matrix) += prod(rodrigues_variation_matrix, current_lambda_reference_rodrigues_lambda_derivative_product);
        noalias(rodrigues_variation_lambda_rodrigues_derivative_lambda_matrix) += prod(rodrigues_variation_matrix, current_lambda_reference_rodrigues_derivative_lambda_product);
        noalias(rodrigues_variation_lambda_derivative_rodrigues_lambda_matrix) += prod(rodrigues_variation_matrix, current_lambda_derivative_reference_rodrigues_lambda_product);
        noalias(rodrigues_derivative_variation_lambda_rodrigues_lambda_matrix) += prod(rodrigues_derivative_variation_matrix, current_lambda_reference_rodrigues_lambda_product);
        AddMatrixTimesVariation(rodrigues_lambda_variation_rodrigues_lambda_derivative_matrix, current_rodrigues_matrix, lambda_variation_rodrigues_lambda_derivative_matrix, mNumberOfDofs);
        AddMatrixTimesVariation(rodrigues_lambda_variation_rodrigues_derivative_lambda_matrix, current_rodrigues_matrix, lambda_variation_rodrigues_derivative_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesVariation(rodrigues_lambda_derivative_variation_rodrigues_lambda_matrix, current_rodrigues_matrix, lambda_derivative_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesVariation(rodrigues_derivative_lambda_variation_rodrigues_lambda_matrix, current_rodrigues_matrix_derivative, lambda_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesVariation(rodrigues_lambda_variation_rodrigues_lambda_matrix, current_rodrigues_matrix, lambda_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        noalias(rodrigues_variation_lambda_rodrigues_lambda_matrix) += prod(rodrigues_variation_matrix, current_lambda_reference_rodrigues_lambda_product);



        rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix.clear();
        rodrigues_lambda_rodrigues_lambda_variation_matrix.clear();
        rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix = rodrigues_variation_lambda_rodrigues_lambda_derivative_matrix + rodrigues_variation_lambda_rodrigues_derivative_lambda_matrix + rodrigues_variation_lambda_derivative_rodrigues_lambda_matrix + rodrigues_derivative_variation_lambda_rodrigues_lambda_matrix
            + rodrigues_lambda_variation_rodrigues_lambda_derivative_matrix + rodrigues_lambda_variation_rodrigues_derivative_lambda_matrix + rodrigues_lambda_derivative_variation_rodrigues_lambda_matrix + rodrigues_derivative_lambda_variation_rodrigues_lambda_matrix;
        rodrigues_lambda_rodrigues_lambda_variation_matrix = rodrigues_variation_lambda_rodrigues_lambda_matrix + rodrigues_lambda_variation_rodrigues_lambda_matrix;

        lambda_second_variation_rodrigues_lambda_matrix.clear();
        lambda_derivative_second_variation_rodrigues_lambda_matrix.clear();
        lambda_second_variation_rodrigues_derivative_lambda_matrix.clear();
        lambda_second_variation_rodrigues_lambda_derivative_matrix.clear();

        AddSecondVariationTimesMatrix(lambda_second_variation_rodrigues_lambda_matrix, lambda_second_variation_matrix, reference_rodrigues_lambda_product, mNumberOfDofs);
        AddSecondVariationTimesMatrix(lambda_derivative_second_variation_rodrigues_lambda_matrix, lambda_derivative_second_variation_matrix, reference_rodrigues_lambda_product, mNumberOfDofs);
        AddSecondVariationTimesMatrix(lambda_second_variation_rodrigues_derivative_lambda_matrix, lambda_second_variation_matrix, reference_rodrigues_derivative_lambda_product, mNumberOfDofs);
        AddSecondVariationTimesMatrix(lambda_second_variation_rodrigues_lambda_derivative_matrix, lambda_second_variation_matrix, reference_rodrigues_lambda_derivative_product, mNumberOfDofs);

        // 5. Second variation of T,1 = (R_t Lambda P),1: all product-rule terms of R_t, Lambda and P with two
        //    variations distributed over R_t and Lambda (P is independent of the dofs), Eq. (A.2) and (A.4)
        Matrix current_rodrigues_derivative_lambda_second_variation_reference_rodrigues_lambda_term;
        current_rodrigues_derivative_lambda_second_variation_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_derivative_lambda_second_variation_reference_rodrigues_lambda_term.clear();
        Matrix current_rodrigues_lambda_derivative_second_variation_reference_rodrigues_lambda_term;
        current_rodrigues_lambda_derivative_second_variation_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_lambda_derivative_second_variation_reference_rodrigues_lambda_term.clear();
        Matrix current_rodrigues_lambda_second_variation_reference_rodrigues_derivative_lambda_term;
        current_rodrigues_lambda_second_variation_reference_rodrigues_derivative_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_lambda_second_variation_reference_rodrigues_derivative_lambda_term.clear();
        Matrix current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_derivative_term;
        current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_derivative_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_derivative_term.clear();

        Matrix current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term;
        current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term.clear();
        Matrix current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term;
        current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term.clear();
        Matrix current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term;
        current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term.clear();
        Matrix current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term;
        current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term.clear();

        Matrix current_rodrigues_derivative_second_variation_lambda_reference_rodrigues_lambda_term;
        current_rodrigues_derivative_second_variation_lambda_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_derivative_second_variation_lambda_reference_rodrigues_lambda_term.clear();
        Matrix current_rodrigues_second_variation_lambda_derivative_reference_rodrigues_lambda_term;
        current_rodrigues_second_variation_lambda_derivative_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_second_variation_lambda_derivative_reference_rodrigues_lambda_term.clear();
        Matrix current_rodrigues_second_variation_lambda_reference_rodrigues_derivative_lambda_term;
        current_rodrigues_second_variation_lambda_reference_rodrigues_derivative_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_second_variation_lambda_reference_rodrigues_derivative_lambda_term.clear();
        Matrix current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_derivative_term;
        current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_derivative_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_derivative_term.clear();
        Matrix current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_term;
        current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_term.clear();
        Matrix current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term;
        current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term.clear();
        Matrix current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_term;
        current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_term.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_term.clear();

        AddMatrixTimesSecondVariation(current_rodrigues_derivative_lambda_second_variation_reference_rodrigues_lambda_term, current_rodrigues_matrix_derivative, lambda_second_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesSecondVariation(current_rodrigues_lambda_derivative_second_variation_reference_rodrigues_lambda_term, current_rodrigues_matrix, lambda_derivative_second_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesSecondVariation(current_rodrigues_lambda_second_variation_reference_rodrigues_derivative_lambda_term, current_rodrigues_matrix, lambda_second_variation_rodrigues_derivative_lambda_matrix, mNumberOfDofs);
        AddMatrixTimesSecondVariation(current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_derivative_term, current_rodrigues_matrix, lambda_second_variation_rodrigues_lambda_derivative_matrix, mNumberOfDofs);
        AddVariationTimesVariation(current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term, rodrigues_derivative_variation_matrix, lambda_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddVariationTimesVariation(current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term, rodrigues_variation_matrix, lambda_derivative_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddVariationTimesVariation(current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term, rodrigues_variation_matrix, lambda_variation_rodrigues_derivative_lambda_matrix, mNumberOfDofs);
        AddVariationTimesVariation(current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term, rodrigues_variation_matrix, lambda_variation_rodrigues_lambda_derivative_matrix, mNumberOfDofs);
        AddSecondVariationTimesMatrix(current_rodrigues_derivative_second_variation_lambda_reference_rodrigues_lambda_term, rodrigues_derivative_second_variation_matrix, current_lambda_reference_rodrigues_lambda_product, mNumberOfDofs);
        AddSecondVariationTimesMatrix(current_rodrigues_second_variation_lambda_derivative_reference_rodrigues_lambda_term, rodrigues_second_variation_matrix, current_lambda_derivative_reference_rodrigues_lambda_product, mNumberOfDofs);
        AddSecondVariationTimesMatrix(current_rodrigues_second_variation_lambda_reference_rodrigues_derivative_lambda_term, rodrigues_second_variation_matrix, current_lambda_reference_rodrigues_derivative_lambda_product, mNumberOfDofs);
        AddSecondVariationTimesMatrix(current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_derivative_term, rodrigues_second_variation_matrix, current_lambda_reference_rodrigues_lambda_derivative_product, mNumberOfDofs);
        AddMatrixTimesSecondVariation(current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_term, current_rodrigues_matrix, lambda_second_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddVariationTimesVariation(current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term, rodrigues_variation_matrix, lambda_variation_rodrigues_lambda_matrix, mNumberOfDofs);
        AddSecondVariationTimesMatrix(current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_term, rodrigues_second_variation_matrix, current_lambda_reference_rodrigues_lambda_product, mNumberOfDofs);

        Matrix cross_section_transformation_derivative_second_variation;
        cross_section_transformation_derivative_second_variation.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        cross_section_transformation_derivative_second_variation.clear();

        Vector3d normal_vector;
        normal_vector.clear();
        Vector3d binormal_vector;
        binormal_vector.clear();

        noalias(normal_vector) += prod(cross_section_transformation, reference_normal);
        noalias(binormal_vector) += prod(cross_section_transformation, reference_binormal);

        Vector normal_vector_variation;
        normal_vector_variation.resize(3 * mNumberOfDofs);
        normal_vector_variation.clear();
        Vector binormal_vector_variation;
        binormal_vector_variation.resize(3 * mNumberOfDofs);
        binormal_vector_variation.clear();

        noalias(normal_vector_variation) += prod(rodrigues_lambda_rodrigues_lambda_variation_matrix, reference_normal);
        noalias(binormal_vector_variation) += prod(rodrigues_lambda_rodrigues_lambda_variation_matrix, reference_binormal);

        // 6. Second variation of T = R_t Lambda P, needed for the second variation of the axes a_alpha
        Matrix cross_section_transformation_second_variation;
        cross_section_transformation_second_variation.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        cross_section_transformation_second_variation.clear();

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int column_component = 0;column_component < 3;column_component++)
            {
                for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
                {
                    for (size_t  second_dof_index  = 0;second_dof_index < mNumberOfDofs;second_dof_index++)
                    {
                        cross_section_transformation_second_variation(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) +=
                            current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index)
                            + current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
                        cross_section_transformation_derivative_second_variation(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) +=
                            current_rodrigues_derivative_lambda_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_lambda_derivative_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_lambda_second_variation_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index)
                            + current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index)
                            + current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index)
                            + current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index)
                            + current_rodrigues_derivative_second_variation_lambda_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_second_variation_lambda_derivative_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_second_variation_lambda_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
                    }
                }
            }
        }

        Matrix normal_vector_second_variation;
        normal_vector_second_variation.resize(3 * mNumberOfDofs, mNumberOfDofs);
        normal_vector_second_variation.clear();
        Matrix binormal_vector_second_variation;
        binormal_vector_second_variation.resize(3 * mNumberOfDofs, mNumberOfDofs);
        binormal_vector_second_variation.clear();

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
            {
                for (size_t  second_dof_index  = 0;second_dof_index < mNumberOfDofs;second_dof_index++)
                {
                    for (size_t  product_index = 0;product_index < 3;product_index++)
                    {
                        normal_vector_second_variation(row_component * mNumberOfDofs + dof_index, second_dof_index) += cross_section_transformation_second_variation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_normal(product_index);
                        binormal_vector_second_variation(row_component * mNumberOfDofs + dof_index, second_dof_index) += cross_section_transformation_second_variation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_binormal(product_index);
                    }
                }
            }
        }

        // 7. Second variations of the generalized strains, Eq. (70):
        //    bending b_alpha = a_alpha,1 . a1 and torsion t_alpha = a_beta,1 . a_alpha
        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (size_t  product_index = 0;product_index < 3;product_index++)
            {
                for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
                {
                    for (size_t  second_dof_index  = 0;second_dof_index < mNumberOfDofs;second_dof_index++)
                    {
                        const double first_base_vector_variation_dof =
                            row_component == dof_index % mDofsPerNode
                                ? shape_function_derivatives[dof_index / mDofsPerNode]
                                : 0.0;

                        const double first_base_vector_variation_second_dof =
                            row_component == second_dof_index % mDofsPerNode
                                ? shape_function_derivatives[second_dof_index / mDofsPerNode]
                                : 0.0;

                        normal_curvature_second_variation(dof_index, second_dof_index) += cross_section_transformation_derivative_second_variation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_normal(product_index) * current_first_base_vector(row_component)
                            + rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index) * first_base_vector_variation_second_dof
                            + rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + second_dof_index, product_index) * reference_normal(product_index) * first_base_vector_variation_dof;

                        binormal_curvature_second_variation(dof_index, second_dof_index) += cross_section_transformation_derivative_second_variation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_binormal(product_index) * current_first_base_vector(row_component)
                            + rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index) * first_base_vector_variation_second_dof
                            + rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + second_dof_index, product_index) * reference_binormal(product_index) * first_base_vector_variation_dof;

                        normal_torsion_second_variation(dof_index, second_dof_index) += cross_section_transformation_derivative_second_variation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_binormal(product_index) * normal_vector(row_component)
                            + rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index) * normal_vector_variation(row_component * mNumberOfDofs + second_dof_index)
                            + rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + second_dof_index, product_index) * reference_binormal(product_index) * normal_vector_variation(row_component * mNumberOfDofs + dof_index)
                            + cross_section_transformation_derivative(row_component, product_index) * reference_binormal(product_index) * normal_vector_second_variation(row_component * mNumberOfDofs + dof_index, second_dof_index);

                        binormal_torsion_second_variation(dof_index, second_dof_index) += cross_section_transformation_derivative_second_variation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_normal(product_index) * binormal_vector(row_component)
                            + rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index) * binormal_vector_variation(row_component * mNumberOfDofs + second_dof_index)
                            + rodrigues_lambda_rodrigues_lambda_derivative_variation_matrix(row_component * mNumberOfDofs + second_dof_index, product_index) * reference_normal(product_index) * binormal_vector_variation(row_component * mNumberOfDofs + dof_index)
                            + cross_section_transformation_derivative(row_component, product_index) * reference_normal(product_index) * binormal_vector_second_variation(row_component * mNumberOfDofs + dof_index, second_dof_index);
                    }
                }
            }
        }

        // Initialize G matrices
        rGAxial.resize(mNumberOfDofs, mNumberOfDofs);
        rGBending1.resize(mNumberOfDofs, mNumberOfDofs);
        rGBending2.resize(mNumberOfDofs, mNumberOfDofs);
        rGTorsion1.resize(mNumberOfDofs, mNumberOfDofs);
        rGTorsion2.resize(mNumberOfDofs, mNumberOfDofs);
        rGAxial.clear();
        rGBending1.clear();
        rGBending2.clear();
        rGTorsion1.clear();
        rGTorsion2.clear();

        noalias(rGAxial) = axial_second_variation;
        noalias(rGBending1) = normal_curvature_second_variation;
        noalias(rGBending2) = binormal_curvature_second_variation;
        noalias(rGTorsion1) = normal_torsion_second_variation;
        noalias(rGTorsion2) = binormal_torsion_second_variation;

        KRATOS_CATCH("")

    }
}
