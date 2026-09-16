//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti, Ricky Aristio
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/non_linear_bernoulli_beam_element_3D.h"
#include "utilities/math_utils.h"

namespace Kratos
{
    void NonLinearBernoulliBeamElement3D::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        InitializeMaterial();


        mDofsPerNode = 4;
        mNumberOfDofs  = this->GetGeometry().size() * mDofsPerNode;
        mSMatrixRodriguesVariation.resize(mNumberOfDofs * 3, 3);
        mSMatrixLambdaVariation.resize(mNumberOfDofs * 3, 3);
        mSMatrixLambdaVariationRodriguesLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesLambdaVariationRodriguesLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesVariationLambdaRodriguesLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesLambdaRodriguesLambdaVariation.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesDerivativeVariation.resize(mNumberOfDofs * 3, 3);
        mSMatrixLambdaDerivativeVariation.resize(mNumberOfDofs * 3, 3);
        mSMatrixLambdaDerivativeVariationRodriguesLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixLambdaVariationRodriguesDerivativeLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixLambdaVariationRodriguesLambdaDerivative.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesDerivativeLambdaVariationRodriguesLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesLambdaDerivativeVariationRodriguesLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesLambdaVariationRodriguesDerivativeLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesLambdaVariationRodriguesLambdaDerivative.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesDerivativeVariationLambdaRodriguesLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesVariationLambdaDerivativeRodriguesLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesVariationLambdaRodriguesDerivativeLambda.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesVariationLambdaRodriguesLambdaDerivative.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation.resize(mNumberOfDofs * 3, 3);
        mSMatrixRodriguesSecondVariation.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatrixLambdaSecondVariation.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatrixLambdaSecondVariationRodriguesLambda.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatrixRodriguesDerivativeSecondVariation.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatrixLambdaDerivativeSecondVariation.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatrixLambdaDerivativeSecondVariationRodriguesLambda.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatrixLambdaSecondVariationRodriguesDerivativeLambda.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatrixLambdaSecondVariationRodriguesLambdaDerivative.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
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

        rKinematicVariables.B = sqrt(
            inner_prod(
                rKinematicVariables.R2,
                rKinematicVariables.R2)
            - pow(
                reference_base_vector_product /
                    rKinematicVariables.A,
                2));

        rKinematicVariables.b = sqrt(
            inner_prod(
                rKinematicVariables.r2,
                rKinematicVariables.r2)
            - pow(
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

        Vector stress_axial = ZeroVector(5);
        Vector stress_bending_1 = ZeroVector(5);
        Vector stress_bending_2 = ZeroVector(5);
        Vector stress_torsion_1 = ZeroVector(5);
        Vector stress_torsion_2 = ZeroVector(5);

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

            if (CalculateStiffnessMatrixFlag == true) {
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

            if (CalculateResidualVectorFlag == true) {
                // Map stress components to different behaviors
                stress_axial[0] =
                    constitutive_variables.StressVector[0];

                stress_bending_1[1] =
                    constitutive_variables.StressVector[1];

                stress_bending_2[2] =
                    constitutive_variables.StressVector[2];

                stress_torsion_1[3] =
                    constitutive_variables.StressVector[3];

                stress_torsion_2[4] =
                    constitutive_variables.StressVector[4];

                // Assemble internal forces
                rRightHandSideVector -=
                    inverse_a4 *
                    integration_weight *
                    prod(
                        trans(b_axial),
                        stress_axial);

                rRightHandSideVector -=
                    inverse_a4 *
                    integration_weight *
                    prod(
                        trans(b_bending_1),
                        stress_bending_1);

                rRightHandSideVector -=
                    inverse_a4 *
                    integration_weight *
                    prod(
                        trans(b_bending_2),
                        stress_bending_2);

                rRightHandSideVector -=
                    inverse_a2 *
                    integration_weight *
                    prod(
                        trans(b_torsion_1),
                        stress_torsion_1);

                rRightHandSideVector -=
                    inverse_a2 *
                    integration_weight *
                    prod(
                        trans(b_torsion_2),
                        stress_torsion_2);
            }
        }

        KRATOS_CATCH("")
    }

        
    int NonLinearBernoulliBeamElement3D::Check(
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        Check(rCurrentProcessInfo);

        const double numerical_limit =
            std::numeric_limits<double>::epsilon();

        KRATOS_ERROR_IF(
            !GetProperties().Has(I_T) ||
            GetProperties()[I_T] <= numerical_limit)
            << "Please provide a reasonable value for \"I_T\" "
            << "(torsional moment of inertia) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(
            !GetProperties().Has(I_N) ||
            GetProperties()[I_N] <= numerical_limit)
            << "Please provide a reasonable value for \"I_N\" "
            << "(bending moment of inertia about N axis) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(
            !GetProperties().Has(I_V) ||
            GetProperties()[I_V] <= numerical_limit)
            << "Please provide a reasonable value for \"I_V\" "
            << "(bending moment of inertia about V axis) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(
            !GetProperties().Has(LOCAL_AXIS_ORIENTATION))
            << "\"LOCAL_AXIS_ORIENTATION\" not provided for element #"
            << Id() << std::endl;

        if (GetProperties().Has(T_0)) {
            const Vector3d& r_reference_tangent =
                GetProperties()[T_0];

            KRATOS_ERROR_IF(
                norm_2(r_reference_tangent) <= numerical_limit)
                << "\"T_0\" reference director vector must have "
                << "non-zero magnitude for element #"
                << Id() << std::endl;
        }

        if (GetProperties().Has(N_0)) {
            const Vector3d& r_reference_normal =
                GetProperties()[N_0];

            KRATOS_ERROR_IF(
                norm_2(r_reference_normal) <= numerical_limit)
                << "\"N_0\" reference director vector must have "
                << "non-zero magnitude for element #"
                << Id() << std::endl;
        }

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

        noalias(rConstitutiveVariables.StressVector) =
            prod(
                trans(rConstitutiveVariables.ConstitutiveMatrix),
                rConstitutiveVariables.StrainVector);
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
        const Vector3d initial_tangent = GetProperties()[T_0];
        const Vector3d reference_base_vector = rKinematicVariables.R1;

        const Vector3d current_reference_tangent =
            reference_base_vector / norm_2(reference_base_vector);

        const Vector3d normalized_initial_tangent =
            initial_tangent / norm_2(initial_tangent);

        Vector3d normal = rNormal;

        const double tangent_normal_projection =
            inner_prod(current_reference_tangent, normal);

        // Project the normal onto the plane normal to the reference tangent.
        normal -= tangent_normal_projection * current_reference_tangent;
        normal /= norm_2(normal);

        const Vector3d reference_principal_axis = normal;

        Matrix3d lambda_matrix;
        ComputeLambdaMatrix(
            normalized_initial_tangent,
            current_reference_tangent,
            lambda_matrix);

        const Vector3d rotated_reference_axis =
            prod(lambda_matrix, reference_principal_axis);

        const Vector3d orientation_reference =
            MathUtils<double>::CrossProduct(reference_base_vector, rotated_reference_axis);

        const double axis_dot_product =
            inner_prod(rotated_reference_axis, normal);

        const double orientation_dot_product =
            inner_prod(normal, orientation_reference);

        double cos_theta =
            axis_dot_product /
            (norm_2(rotated_reference_axis) * norm_2(normal));

        if (std::abs(1.0 - std::abs(cos_theta)) < 1e-9) {
            cos_theta = MathUtils<double>::Sign(cos_theta);
        }

        double phi = std::acos(cos_theta);

        const double pi = 4.0 * std::atan(1.0);

        if (orientation_dot_product < -1e-12) {
            phi = 2.0 * pi - phi;
        }

        return phi;
    }
            

    void NonLinearBernoulliBeamElement3D::ComputeReferenceCrossSectionGeometry(
        KinematicVariables& rKinematicVariables)
    {
        const Vector3d reference_first_base_vector = rKinematicVariables.R1;
        const Vector3d reference_second_base_vector = rKinematicVariables.R2;
        const Vector3d initial_reference_tangent = GetProperties()[T_0];

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
                / pow(reference_first_base_vector_norm, 3)
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
        rKinematicVariables.N0 = GetProperties()[N_0];

        const double initial_reference_tangent_norm =
            norm_2(initial_reference_tangent);

        const Vector3d normalized_initial_reference_tangent =
            initial_reference_tangent / initial_reference_tangent_norm;

        rKinematicVariables.N0 -=
            inner_prod(
                normalized_initial_reference_tangent,
                rKinematicVariables.N0)
            * normalized_initial_reference_tangent;

        rKinematicVariables.V0 =
            MathUtils<double>::CrossProduct(
                normalized_initial_reference_tangent,
                rKinematicVariables.N0);

        rKinematicVariables.N0 /=
            norm_2(rKinematicVariables.N0);

        rKinematicVariables.V0 /=
            norm_2(rKinematicVariables.V0);

        Vector3d transformed_reference_normal;
        transformed_reference_normal.clear();

        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                transformed_reference_normal[i] +=
                    reference_lambda_matrix(i, j)
                    * rKinematicVariables.N0[j];
            }
        }

        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                rKinematicVariables.n(i) +=
                    reference_rodrigues_matrix(i, j)
                    * transformed_reference_normal[j];
            }
        }

        rKinematicVariables.n /=
            norm_2(rKinematicVariables.n);

        rKinematicVariables.v =
            MathUtils<double>::CrossProduct(reference_tangent, rKinematicVariables.n);

        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                for (IndexType k = 0; k < 3; ++k) {
                    cross_section_transformation_derivative(i, j) +=
                        reference_rodrigues_matrix_derivative(i, k)
                        * reference_lambda_matrix(k, j);

                    cross_section_transformation_derivative(i, j) +=
                        reference_rodrigues_matrix(i, k)
                        * reference_lambda_matrix_derivative(k, j);
                }
            }
        }

        Vector3d a21_vector;
        Vector3d a31_vector;
        a21_vector.clear();
        a31_vector.clear();

        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                a21_vector[i] +=
                    cross_section_transformation_derivative(i, j)
                    * rKinematicVariables.N0[j];

                a31_vector[i] +=
                    cross_section_transformation_derivative(i, j)
                    * rKinematicVariables.V0[j];
            }
        }

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
        const Vector3d initial_reference_tangent = GetProperties()[T_0];

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
                / pow(current_first_base_vector_norm, 3)
                * current_first_base_vector;

        const Vector3d reference_tangent =
            reference_first_base_vector / reference_first_base_vector_norm;

        const Vector3d reference_tangent_derivative =
            reference_second_base_vector / reference_first_base_vector_norm
            - inner_prod(reference_first_base_vector, reference_second_base_vector)
                / pow(reference_first_base_vector_norm, 3)
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

        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                for (IndexType k = 0; k < 3; ++k) {
                    reference_rodrigues_lambda_product(i, j) +=
                        reference_rodrigues_matrix(i, k)
                        * reference_lambda_matrix(k, j);

                    reference_rodrigues_lambda_derivative_product(i, j) +=
                        reference_rodrigues_matrix(i, k)
                        * reference_lambda_matrix_derivative(k, j);

                    reference_rodrigues_derivative_lambda_product(i, j) +=
                        reference_rodrigues_matrix_derivative(i, k)
                        * reference_lambda_matrix(k, j);
                }
            }
        }

        Matrix3d current_lambda_reference_rodrigues_lambda_product;
        Matrix3d current_lambda_reference_rodrigues_lambda_derivative_product;
        Matrix3d current_lambda_reference_rodrigues_derivative_lambda_product;
        Matrix3d current_lambda_derivative_reference_rodrigues_lambda_product;
        current_lambda_reference_rodrigues_lambda_product.clear();
        current_lambda_reference_rodrigues_lambda_derivative_product.clear();
        current_lambda_reference_rodrigues_derivative_lambda_product.clear();
        current_lambda_derivative_reference_rodrigues_lambda_product.clear();

        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                for (IndexType k = 0; k < 3; ++k) {
                    current_lambda_reference_rodrigues_lambda_product(i, j) +=
                        current_lambda_matrix(i, k)
                        * reference_rodrigues_lambda_product(k, j);

                    current_lambda_reference_rodrigues_lambda_derivative_product(i, j) +=
                        current_lambda_matrix(i, k)
                        * reference_rodrigues_lambda_derivative_product(k, j);

                    current_lambda_reference_rodrigues_derivative_lambda_product(i, j) +=
                        current_lambda_matrix(i, k)
                        * reference_rodrigues_derivative_lambda_product(k, j);

                    current_lambda_derivative_reference_rodrigues_lambda_product(i, j) +=
                        current_lambda_matrix_derivative(i, k)
                        * reference_rodrigues_lambda_product(k, j);
                }
            }
        }

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

        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                for (IndexType k = 0; k < 3; ++k) {
                    current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term(i, j) +=
                        current_rodrigues_matrix(i, k)
                        * current_lambda_reference_rodrigues_lambda_derivative_product(k, j);

                    current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term(i, j) +=
                        current_rodrigues_matrix(i, k)
                        * current_lambda_reference_rodrigues_derivative_lambda_product(k, j);

                    current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term(i, j) +=
                        current_rodrigues_matrix(i, k)
                        * current_lambda_derivative_reference_rodrigues_lambda_product(k, j);

                    current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term(i, j) +=
                        current_rodrigues_matrix_derivative(i, k)
                        * current_lambda_reference_rodrigues_lambda_product(k, j);

                    cross_section_transformation(i, j) +=
                        current_rodrigues_matrix(i, k)
                        * current_lambda_reference_rodrigues_lambda_product(k, j);
                }
            }
        }

        const Matrix3d cross_section_transformation_derivative =
            current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term
            + current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term
            + current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term
            + current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;

        Vector3d a21_vector;
        Vector3d a31_vector;
        a21_vector.clear();
        a31_vector.clear();

        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                a21_vector[i] +=
                    cross_section_transformation_derivative(i, j)
                    * rKinematicVariables.N0[j];

                a31_vector[i] +=
                    cross_section_transformation_derivative(i, j)
                    * rKinematicVariables.V0[j];

                rKinematicVariables.n[i] +=
                    cross_section_transformation(i, j)
                    * rKinematicVariables.N0[j];

                rKinematicVariables.v[i] +=
                    cross_section_transformation(i, j)
                    * rKinematicVariables.V0[j];
            }
        }

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
        Vector3d initial_reference_tangent = GetProperties()[T_0];

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
        current_tangent_derivative = current_second_base_vector / norm_2(current_first_base_vector) - inner_prod(current_first_base_vector, current_second_base_vector) / pow(norm_2(current_first_base_vector), 3) * current_first_base_vector;
        reference_tangent = reference_first_base_vector / norm_2(reference_first_base_vector);
        reference_tangent_derivative = reference_second_base_vector / norm_2(reference_first_base_vector) - inner_prod(reference_first_base_vector, reference_second_base_vector) / pow(norm_2(reference_first_base_vector), 3) * reference_first_base_vector;
        ComputeTangentVariation(shape_function_derivatives, current_first_base_vector, current_tangent_variation);
        ComputeTangentDerivativeVariation(shape_function_derivatives, shape_function_second_derivatives, current_first_base_vector, current_second_base_vector, current_tangent_derivative_variation);

        ComputeLambdaMatrix(reference_tangent, current_tangent, current_lambda_matrix);
        ComputeLambdaMatrixDerivative(reference_tangent, current_tangent, reference_tangent_derivative, current_tangent_derivative, current_lambda_matrix_derivative);
        ComputeLambdaMatrixVariation(reference_tangent, current_tangent, current_tangent_variation, mNumberOfDofs, mDofsPerNode, mSMatrixLambdaVariation);
        ComputeLambdaMatrixDerivativeVariation(reference_tangent, current_tangent, reference_tangent_derivative, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, mNumberOfDofs, mDofsPerNode, mSMatrixLambdaDerivativeVariation);
        ComputeLambdaMatrix(initial_reference_tangent, reference_tangent, reference_lambda_matrix);
        ComputeLambdaMatrixDerivative(initial_reference_tangent, reference_tangent, initial_reference_tangent_derivative, reference_tangent_derivative, reference_lambda_matrix_derivative);

        ComputeRodriguesMatrix(current_tangent, current_twist_angle, current_rodrigues_matrix);
        ComputeRodriguesMatrixDerivative(current_tangent, current_tangent_derivative, current_twist_angle, current_twist_angle_derivative, current_rodrigues_matrix_derivative);
        ComputeRodriguesMatrixVariation(current_tangent, current_tangent_variation, shape_functions, current_twist_angle, mNumberOfDofs, mDofsPerNode, mSMatrixRodriguesVariation);
        ComputeRodriguesMatrixDerivativeVariation(current_tangent, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, shape_functions, shape_function_derivatives, current_twist_angle, current_twist_angle_derivative, mNumberOfDofs, mDofsPerNode, mSMatrixRodriguesDerivativeVariation);
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

        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (int column_component = 0; column_component < 3; column_component++)
            {
                for (int product_index = 0; product_index < 3; product_index++)
                {
                    reference_rodrigues_lambda_product(row_component, column_component) += reference_rodrigues_matrix(row_component, product_index) * reference_lambda_matrix(product_index, column_component);
                    reference_rodrigues_lambda_derivative_product(row_component, column_component) += reference_rodrigues_matrix(row_component, product_index) * reference_lambda_matrix_derivative(product_index, column_component);
                    reference_rodrigues_derivative_lambda_product(row_component, column_component) += reference_rodrigues_matrix_derivative(row_component, product_index) * reference_lambda_matrix(product_index, column_component);
                }
            }
        }
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

        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (int column_component = 0; column_component < 3; column_component++)
            {
                for (int product_index = 0; product_index < 3; product_index++)
                {
                    current_lambda_reference_rodrigues_lambda_product(row_component, column_component) += current_lambda_matrix(row_component, product_index) * reference_rodrigues_lambda_product(product_index, column_component);
                    current_lambda_reference_rodrigues_lambda_derivative_product(row_component, column_component) += current_lambda_matrix(row_component, product_index) * reference_rodrigues_lambda_derivative_product(product_index, column_component);
                    current_lambda_reference_rodrigues_derivative_lambda_product(row_component, column_component) += current_lambda_matrix(row_component, product_index) * reference_rodrigues_derivative_lambda_product(product_index, column_component);
                    current_lambda_derivative_reference_rodrigues_lambda_product(row_component, column_component) += current_lambda_matrix_derivative(row_component, product_index) * reference_rodrigues_lambda_product(product_index, column_component);
                }
            }
        }

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

        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (int column_component = 0; column_component < 3; column_component++)
            {
                for (int product_index = 0; product_index < 3; product_index++)
                {
                    current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term(row_component, column_component) += current_rodrigues_matrix(row_component, product_index) * current_lambda_reference_rodrigues_lambda_derivative_product(product_index, column_component);
                    current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term(row_component, column_component) += current_rodrigues_matrix(row_component, product_index) * current_lambda_reference_rodrigues_derivative_lambda_product(product_index, column_component);
                    current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term(row_component, column_component) += current_rodrigues_matrix(row_component, product_index) * current_lambda_derivative_reference_rodrigues_lambda_product(product_index, column_component);
                    current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term(row_component, column_component) += current_rodrigues_matrix_derivative(row_component, product_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                    cross_section_transformation(row_component, column_component) += current_rodrigues_matrix(row_component, product_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                }
            }
        }

        Matrix3d cross_section_transformation_derivative;
        cross_section_transformation_derivative.clear();
        cross_section_transformation_derivative = current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term + current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term + current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term + current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;

        mSMatrixLambdaVariationRodriguesLambdaDerivative.clear();
        mSMatrixLambdaVariationRodriguesDerivativeLambda.clear();
        mSMatrixLambdaVariationRodriguesLambda.clear();
        mSMatrixLambdaDerivativeVariationRodriguesLambda.clear();

        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (int column_component = 0; column_component < 3; column_component++)
            {
                for (int product_index = 0; product_index < 3; product_index++)
                {
                    for (size_t dof_index = 0; dof_index < mNumberOfDofs; dof_index++)
                    {
                        mSMatrixLambdaVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_rodrigues_lambda_product(product_index, column_component);
                        mSMatrixLambdaVariationRodriguesLambdaDerivative(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_rodrigues_lambda_derivative_product(product_index, column_component);
                        mSMatrixLambdaVariationRodriguesDerivativeLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_rodrigues_derivative_lambda_product(product_index, column_component);
                        mSMatrixLambdaDerivativeVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixLambdaDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_rodrigues_lambda_product(product_index, column_component);
                    }
                }
            }
        }

        Matrix current_lambda_reference_rodrigues_lambda_product_derivative_variation;
        current_lambda_reference_rodrigues_lambda_product_derivative_variation.resize(3 * mNumberOfDofs, 3);
        current_lambda_reference_rodrigues_lambda_product_derivative_variation.clear();

        current_lambda_reference_rodrigues_lambda_product_derivative_variation = mSMatrixLambdaVariationRodriguesLambdaDerivative + mSMatrixLambdaVariationRodriguesDerivativeLambda + mSMatrixLambdaDerivativeVariationRodriguesLambda;

        mSMatrixRodriguesVariationLambdaRodriguesLambdaDerivative.clear();
        mSMatrixRodriguesVariationLambdaRodriguesDerivativeLambda.clear();
        mSMatrixRodriguesVariationLambdaDerivativeRodriguesLambda.clear();
        mSMatrixRodriguesDerivativeVariationLambdaRodriguesLambda.clear();
        mSMatrixRodriguesDerivativeLambdaVariationRodriguesLambda.clear();
        mSMatrixRodriguesLambdaDerivativeVariationRodriguesLambda.clear();
        mSMatrixRodriguesLambdaVariationRodriguesDerivativeLambda.clear();
        mSMatrixRodriguesLambdaVariationRodriguesLambdaDerivative.clear();
        mSMatrixRodriguesLambdaVariationRodriguesLambda.clear();
        mSMatrixRodriguesVariationLambdaRodriguesLambda.clear();

        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (int column_component = 0; column_component < 3; column_component++)
            {
                for (int product_index = 0; product_index < 3; product_index++)
                {
                    for (size_t dof_index = 0; dof_index < mNumberOfDofs; dof_index++)
                    {
                        mSMatrixRodriguesVariationLambdaRodriguesLambdaDerivative(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_reference_rodrigues_lambda_derivative_product(product_index, column_component);
                        mSMatrixRodriguesVariationLambdaRodriguesDerivativeLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_reference_rodrigues_derivative_lambda_product(product_index, column_component);
                        mSMatrixRodriguesVariationLambdaDerivativeRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_derivative_reference_rodrigues_lambda_product(product_index, column_component);
                        mSMatrixRodriguesDerivativeVariationLambdaRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                        mSMatrixRodriguesLambdaVariationRodriguesLambdaDerivative(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaVariationRodriguesLambdaDerivative(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesLambdaVariationRodriguesDerivativeLambda(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaVariationRodriguesDerivativeLambda(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesLambdaDerivativeVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaDerivativeVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesDerivativeLambdaVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix_derivative(row_component, product_index) * mSMatrixLambdaVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesLambdaVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesVariationLambdaRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                    }
                }
            }
        }

        mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation.clear();
        mSMatrixRodriguesLambdaRodriguesLambdaVariation.clear();

        mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation = mSMatrixRodriguesVariationLambdaRodriguesLambdaDerivative + mSMatrixRodriguesVariationLambdaRodriguesDerivativeLambda + mSMatrixRodriguesVariationLambdaDerivativeRodriguesLambda + mSMatrixRodriguesDerivativeVariationLambdaRodriguesLambda
            + mSMatrixRodriguesLambdaVariationRodriguesLambdaDerivative + mSMatrixRodriguesLambdaVariationRodriguesDerivativeLambda + mSMatrixRodriguesLambdaDerivativeVariationRodriguesLambda + mSMatrixRodriguesDerivativeLambdaVariationRodriguesLambda;
        mSMatrixRodriguesLambdaRodriguesLambdaVariation = mSMatrixRodriguesVariationLambdaRodriguesLambda + mSMatrixRodriguesLambdaVariationRodriguesLambda;

        Vector3d normal_vector;
        normal_vector.clear();
        Vector3d binormal_vector;
        binormal_vector.clear();

        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (int product_index = 0; product_index < 3; product_index++)
            {
                normal_vector(row_component) += cross_section_transformation(row_component, product_index) * reference_normal(product_index);
                binormal_vector(row_component) += cross_section_transformation(row_component, product_index) * reference_binormal(product_index);
            }
        }

        Vector normal_vector_variation;
        normal_vector_variation.resize(3 * mNumberOfDofs);
        normal_vector_variation.clear();
        Vector binormal_vector_variation;
        binormal_vector_variation.resize(3 * mNumberOfDofs);
        binormal_vector_variation.clear();

        for (size_t row_component = 0; row_component < 3; row_component++)
        {
            for (size_t dof_index = 0; dof_index < mNumberOfDofs; dof_index++)
            {
                for (int product_index = 0; product_index < 3; product_index++)
                {
                    normal_vector_variation(row_component * mNumberOfDofs + dof_index) += mSMatrixRodriguesLambdaRodriguesLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index);
                    binormal_vector_variation(row_component * mNumberOfDofs + dof_index) += mSMatrixRodriguesLambdaRodriguesLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index);
                }
            }
        }

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
                    normal_curvature_variation(dof_index) += mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index) * current_first_base_vector[row_component] + cross_section_transformation_derivative(row_component, product_index) * reference_normal(product_index) * current_first_base_vector_variation[row_component * mNumberOfDofs + dof_index];
                    binormal_curvature_variation(dof_index) += mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index) * current_first_base_vector[row_component] + cross_section_transformation_derivative(row_component, product_index) * reference_binormal(product_index) * current_first_base_vector_variation[row_component * mNumberOfDofs + dof_index];
                    normal_torsion_variation(dof_index) += mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index) * normal_vector(row_component) + cross_section_transformation_derivative(row_component, product_index) * reference_binormal(product_index) * normal_vector_variation(row_component * mNumberOfDofs + dof_index);
                    binormal_torsion_variation(dof_index) += mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index) * binormal_vector(row_component) + cross_section_transformation_derivative(row_component, product_index) * reference_normal(product_index) * binormal_vector_variation(row_component * mNumberOfDofs + dof_index);
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
        Vector3d initial_reference_tangent = GetProperties()[T_0];

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
        current_tangent_derivative = current_second_base_vector / norm_2(current_first_base_vector) - inner_prod(current_first_base_vector, current_second_base_vector) / pow(norm_2(current_first_base_vector), 3) * current_first_base_vector;
        reference_tangent = reference_first_base_vector / norm_2(reference_first_base_vector);
        reference_tangent_derivative = reference_second_base_vector / norm_2(reference_first_base_vector) - inner_prod(reference_first_base_vector, reference_second_base_vector) / pow(norm_2(reference_first_base_vector), 3) * reference_first_base_vector;
        ComputeTangentVariation(shape_function_derivatives, current_first_base_vector, current_tangent_variation);
        ComputeTangentDerivativeVariation(shape_function_derivatives, shape_function_second_derivatives, current_first_base_vector, current_second_base_vector, current_tangent_derivative_variation);
        ComputeTangentSecondVariation(shape_function_derivatives, current_first_base_vector, current_tangent_second_variation);
        ComputeTangentDerivativeSecondVariation(shape_function_derivatives, shape_function_second_derivatives, current_first_base_vector, current_second_base_vector, current_tangent_derivative_second_variation);

        ComputeLambdaMatrix(reference_tangent, current_tangent, current_lambda_matrix);
        ComputeLambdaMatrixDerivative(reference_tangent, current_tangent, reference_tangent_derivative, current_tangent_derivative, current_lambda_matrix_derivative);

        ComputeLambdaMatrix(initial_reference_tangent, reference_tangent, reference_lambda_matrix);
        ComputeLambdaMatrixDerivative(initial_reference_tangent, reference_tangent, initial_reference_tangent_derivative, reference_tangent_derivative, reference_lambda_matrix_derivative);
        ComputeLambdaMatrixVariations(reference_tangent, current_tangent, reference_tangent_derivative, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, current_tangent_second_variation, current_tangent_derivative_second_variation, mNumberOfDofs, mDofsPerNode, mSMatrixLambdaVariation, mSMatrixLambdaDerivativeVariation, mSMatrixLambdaSecondVariation, mSMatrixLambdaDerivativeSecondVariation);

        ComputeRodriguesMatrix(current_tangent, current_twist_angle, current_rodrigues_matrix);
        ComputeRodriguesMatrixDerivative(current_tangent, current_tangent_derivative, current_twist_angle, current_twist_angle_derivative, current_rodrigues_matrix_derivative);
        ComputeRodriguesMatrixVariation(current_tangent, current_tangent_variation, shape_functions, current_twist_angle, mNumberOfDofs, mDofsPerNode, mSMatrixRodriguesVariation);
        ComputeRodriguesMatrixDerivativeVariation(current_tangent, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, shape_functions, shape_function_derivatives, current_twist_angle, current_twist_angle_derivative, mNumberOfDofs, mDofsPerNode, mSMatrixRodriguesDerivativeVariation);
        ComputeRodriguesMatrixSecondVariation(current_tangent, current_tangent_variation, current_tangent_second_variation, shape_functions, current_twist_angle, mNumberOfDofs, mDofsPerNode, mSMatrixRodriguesSecondVariation);
        ComputeRodriguesMatrixDerivativeSecondVariation(current_tangent, current_tangent_variation, current_tangent_derivative, current_tangent_derivative_variation, current_tangent_second_variation, current_tangent_derivative_second_variation, shape_functions, shape_function_derivatives, current_twist_angle, current_twist_angle_derivative, mNumberOfDofs, mDofsPerNode, mSMatrixRodriguesDerivativeSecondVariation);

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

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int column_component = 0;column_component < 3;column_component++)
            {
                for (int product_index = 0;product_index < 3;product_index++)
                {
                    reference_rodrigues_lambda_product(row_component, column_component) += reference_rodrigues_matrix(row_component, product_index) * reference_lambda_matrix(product_index, column_component);
                    reference_rodrigues_lambda_derivative_product(row_component, column_component) += reference_rodrigues_matrix(row_component, product_index) * reference_lambda_matrix_derivative(product_index, column_component);
                    reference_rodrigues_derivative_lambda_product(row_component, column_component) += reference_rodrigues_matrix_derivative(row_component, product_index) * reference_lambda_matrix(product_index, column_component);
                }
            }
        }
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

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int column_component = 0;column_component < 3;column_component++)
            {
                for (int product_index = 0;product_index < 3;product_index++)
                {
                    current_lambda_reference_rodrigues_lambda_product(row_component, column_component) += current_lambda_matrix(row_component, product_index) * reference_rodrigues_lambda_product(product_index, column_component);
                    current_lambda_reference_rodrigues_lambda_derivative_product(row_component, column_component) += current_lambda_matrix(row_component, product_index) * reference_rodrigues_lambda_derivative_product(product_index, column_component);
                    current_lambda_reference_rodrigues_derivative_lambda_product(row_component, column_component) += current_lambda_matrix(row_component, product_index) * reference_rodrigues_derivative_lambda_product(product_index, column_component);
                    current_lambda_derivative_reference_rodrigues_lambda_product(row_component, column_component) += current_lambda_matrix_derivative(row_component, product_index) * reference_rodrigues_lambda_product(product_index, column_component);
                }
            }
        }

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

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int column_component = 0;column_component < 3;column_component++)
            {
                for (int product_index = 0;product_index < 3;product_index++)
                {
                    current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term(row_component, column_component) += current_rodrigues_matrix(row_component, product_index) * current_lambda_reference_rodrigues_lambda_derivative_product(product_index, column_component);
                    current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term(row_component, column_component) += current_rodrigues_matrix(row_component, product_index) * current_lambda_reference_rodrigues_derivative_lambda_product(product_index, column_component);
                    current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term(row_component, column_component) += current_rodrigues_matrix(row_component, product_index) * current_lambda_derivative_reference_rodrigues_lambda_product(product_index, column_component);
                    current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term(row_component, column_component) += current_rodrigues_matrix_derivative(row_component, product_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                    cross_section_transformation(row_component, column_component) += current_rodrigues_matrix(row_component, product_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                }
            }
        }

        Matrix3d cross_section_transformation_derivative;
        cross_section_transformation_derivative.clear();
        cross_section_transformation_derivative = current_rodrigues_current_lambda_reference_rodrigues_lambda_derivative_term + current_rodrigues_current_lambda_reference_rodrigues_derivative_lambda_term + current_rodrigues_current_lambda_derivative_reference_rodrigues_lambda_term + current_rodrigues_derivative_current_lambda_reference_rodrigues_lambda_term;

        
        mSMatrixLambdaVariationRodriguesLambdaDerivative.clear();
        mSMatrixLambdaVariationRodriguesDerivativeLambda.clear();
        mSMatrixLambdaVariationRodriguesLambda.clear();
        mSMatrixLambdaDerivativeVariationRodriguesLambda.clear();

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int column_component = 0;column_component < 3;column_component++)
            {
                for (int product_index = 0;product_index < 3;product_index++)
                {
                    for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
                    {
                        mSMatrixLambdaVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_rodrigues_lambda_product(product_index, column_component);
                        mSMatrixLambdaVariationRodriguesLambdaDerivative(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_rodrigues_lambda_derivative_product(product_index, column_component);
                        mSMatrixLambdaVariationRodriguesDerivativeLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_rodrigues_derivative_lambda_product(product_index, column_component);
                        mSMatrixLambdaDerivativeVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixLambdaDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_rodrigues_lambda_product(product_index, column_component);
                    }
                }
            }
        }

        Matrix current_lambda_reference_rodrigues_lambda_product_derivative_variation;
        current_lambda_reference_rodrigues_lambda_product_derivative_variation.resize(3 * mNumberOfDofs, 3);
        current_lambda_reference_rodrigues_lambda_product_derivative_variation.clear();

        current_lambda_reference_rodrigues_lambda_product_derivative_variation = mSMatrixLambdaVariationRodriguesLambdaDerivative + mSMatrixLambdaVariationRodriguesDerivativeLambda + mSMatrixLambdaDerivativeVariationRodriguesLambda;

        mSMatrixRodriguesVariationLambdaRodriguesLambdaDerivative.clear();
        mSMatrixRodriguesVariationLambdaRodriguesDerivativeLambda.clear();
        mSMatrixRodriguesVariationLambdaDerivativeRodriguesLambda.clear();
        mSMatrixRodriguesDerivativeVariationLambdaRodriguesLambda.clear();
        mSMatrixRodriguesDerivativeLambdaVariationRodriguesLambda.clear();
        mSMatrixRodriguesLambdaDerivativeVariationRodriguesLambda.clear();
        mSMatrixRodriguesLambdaVariationRodriguesDerivativeLambda.clear();
        mSMatrixRodriguesLambdaVariationRodriguesLambdaDerivative.clear();
        mSMatrixRodriguesLambdaVariationRodriguesLambda.clear();
        mSMatrixRodriguesVariationLambdaRodriguesLambda.clear();

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int column_component = 0;column_component < 3;column_component++)
            {
                for (int product_index = 0;product_index < 3;product_index++)
                {
                    for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
                    {
                        mSMatrixRodriguesVariationLambdaRodriguesLambdaDerivative(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_reference_rodrigues_lambda_derivative_product(product_index, column_component);
                        mSMatrixRodriguesVariationLambdaRodriguesDerivativeLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_reference_rodrigues_derivative_lambda_product(product_index, column_component);
                        mSMatrixRodriguesVariationLambdaDerivativeRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_derivative_reference_rodrigues_lambda_product(product_index, column_component);
                        mSMatrixRodriguesDerivativeVariationLambdaRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                        mSMatrixRodriguesLambdaVariationRodriguesLambdaDerivative(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaVariationRodriguesLambdaDerivative(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesLambdaVariationRodriguesDerivativeLambda(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaVariationRodriguesDerivativeLambda(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesLambdaDerivativeVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaDerivativeVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesDerivativeLambdaVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix_derivative(row_component, product_index) * mSMatrixLambdaVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesLambdaVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component);
                        mSMatrixRodriguesVariationLambdaRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                    }
                }
            }
        }

        
        
        mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation.clear();
        mSMatrixRodriguesLambdaRodriguesLambdaVariation.clear();
        mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation = mSMatrixRodriguesVariationLambdaRodriguesLambdaDerivative + mSMatrixRodriguesVariationLambdaRodriguesDerivativeLambda + mSMatrixRodriguesVariationLambdaDerivativeRodriguesLambda + mSMatrixRodriguesDerivativeVariationLambdaRodriguesLambda
            + mSMatrixRodriguesLambdaVariationRodriguesLambdaDerivative + mSMatrixRodriguesLambdaVariationRodriguesDerivativeLambda + mSMatrixRodriguesLambdaDerivativeVariationRodriguesLambda + mSMatrixRodriguesDerivativeLambdaVariationRodriguesLambda;
        mSMatrixRodriguesLambdaRodriguesLambdaVariation = mSMatrixRodriguesVariationLambdaRodriguesLambda + mSMatrixRodriguesLambdaVariationRodriguesLambda;

        mSMatrixLambdaSecondVariationRodriguesLambda.clear();
        mSMatrixLambdaDerivativeSecondVariationRodriguesLambda.clear();
        mSMatrixLambdaSecondVariationRodriguesDerivativeLambda.clear();
        mSMatrixLambdaSecondVariationRodriguesLambdaDerivative.clear();

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int column_component = 0;column_component < 3;column_component++)
            {
                for (int product_index = 0;product_index < 3;product_index++)
                {
                    for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
                    {
                        for (size_t  second_dof_index  = 0;second_dof_index < mNumberOfDofs;second_dof_index++)
                        {
                            mSMatrixLambdaSecondVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixLambdaSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_rodrigues_lambda_product(product_index, column_component);
                            mSMatrixLambdaDerivativeSecondVariationRodriguesLambda(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixLambdaDerivativeSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_rodrigues_lambda_product(product_index, column_component);
                            mSMatrixLambdaSecondVariationRodriguesDerivativeLambda(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixLambdaSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_rodrigues_derivative_lambda_product(product_index, column_component);
                            mSMatrixLambdaSecondVariationRodriguesLambdaDerivative(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixLambdaSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_rodrigues_lambda_derivative_product(product_index, column_component);
                        }
                    }
                }
            }
        }

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

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int column_component = 0;column_component < 3;column_component++)
            {
                for (int product_index = 0;product_index < 3;product_index++)
                {
                    for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
                    {
                        for (size_t  second_dof_index  = 0;second_dof_index < mNumberOfDofs;second_dof_index++)
                        {
                            current_rodrigues_derivative_lambda_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += current_rodrigues_matrix_derivative(row_component, product_index) * mSMatrixLambdaSecondVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
                            current_rodrigues_lambda_derivative_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaDerivativeSecondVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
                            current_rodrigues_lambda_second_variation_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaSecondVariationRodriguesDerivativeLambda(product_index * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
                            current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaSecondVariationRodriguesLambdaDerivative(product_index * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
                            current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * mSMatrixLambdaVariationRodriguesLambda(product_index * mNumberOfDofs + second_dof_index, column_component);
                            current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * mSMatrixLambdaDerivativeVariationRodriguesLambda(product_index * mNumberOfDofs + second_dof_index, column_component);
                            current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * mSMatrixLambdaVariationRodriguesDerivativeLambda(product_index * mNumberOfDofs + second_dof_index, column_component);
                            current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * mSMatrixLambdaVariationRodriguesLambdaDerivative(product_index * mNumberOfDofs + second_dof_index, column_component);
                            current_rodrigues_derivative_second_variation_lambda_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesDerivativeSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                            current_rodrigues_second_variation_lambda_derivative_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * current_lambda_derivative_reference_rodrigues_lambda_product(product_index, column_component);
                            current_rodrigues_second_variation_lambda_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * current_lambda_reference_rodrigues_derivative_lambda_product(product_index, column_component);
                            current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * current_lambda_reference_rodrigues_lambda_derivative_product(product_index, column_component);
                            current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += current_rodrigues_matrix(row_component, product_index) * mSMatrixLambdaSecondVariationRodriguesLambda(product_index * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
                            current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesVariation(row_component * mNumberOfDofs + dof_index, product_index) * mSMatrixLambdaVariationRodriguesLambda(product_index * mNumberOfDofs + second_dof_index, column_component);
                            current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += mSMatrixRodriguesSecondVariation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * current_lambda_reference_rodrigues_lambda_product(product_index, column_component);
                        }
                    }
                }
            }
        }

        Matrix cross_section_transformation_derivative_second_variation;
        cross_section_transformation_derivative_second_variation.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
        cross_section_transformation_derivative_second_variation.clear();

        Vector3d normal_vector;
        normal_vector.clear();
        Vector3d binormal_vector;
        binormal_vector.clear();

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (int product_index = 0;product_index < 3;product_index++)
            {
                normal_vector(row_component) += cross_section_transformation(row_component, product_index) * reference_normal(product_index);
                binormal_vector(row_component) += cross_section_transformation(row_component, product_index) * reference_binormal(product_index);
            }
        }

        Vector normal_vector_variation;
        normal_vector_variation.resize(3 * mNumberOfDofs);
        normal_vector_variation.clear();
        Vector binormal_vector_variation;
        binormal_vector_variation.resize(3 * mNumberOfDofs);
        binormal_vector_variation.clear();

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
            {
                for (int product_index = 0;product_index < 3;product_index++)
                {
                    normal_vector_variation(row_component * mNumberOfDofs + dof_index) += mSMatrixRodriguesLambdaRodriguesLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index);
                    binormal_vector_variation(row_component * mNumberOfDofs + dof_index) += mSMatrixRodriguesLambdaRodriguesLambdaVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index);
                }
            }
        }

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
                        cross_section_transformation_second_variation(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index) + current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
                        cross_section_transformation_derivative_second_variation(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) += current_rodrigues_derivative_lambda_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_lambda_derivative_second_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_lambda_second_variation_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_lambda_second_variation_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index)
                            + current_rodrigues_derivative_variation_lambda_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index) + current_rodrigues_variation_lambda_derivative_variation_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index) + current_rodrigues_variation_lambda_variation_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index) + current_rodrigues_variation_lambda_variation_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + second_dof_index, column_component * mNumberOfDofs + dof_index)
                            + current_rodrigues_derivative_second_variation_lambda_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_second_variation_lambda_derivative_reference_rodrigues_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_second_variation_lambda_reference_rodrigues_derivative_lambda_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index) + current_rodrigues_second_variation_lambda_reference_rodrigues_lambda_derivative_term(row_component * mNumberOfDofs + dof_index, column_component * mNumberOfDofs + second_dof_index);
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

        for (size_t  row_component   = 0;row_component < 3;row_component++)
        {
            for (size_t  product_index = 0;product_index < 3;product_index++)
            {
                for (size_t  dof_index  = 0;dof_index < mNumberOfDofs;dof_index++)
                {
                    for (size_t  second_dof_index  = 0;second_dof_index < mNumberOfDofs;second_dof_index++)
                    {
                        normal_curvature_second_variation(dof_index, second_dof_index) += 0;
                        normal_curvature_second_variation(dof_index, second_dof_index) += 0;

                        normal_torsion_second_variation(dof_index, second_dof_index) += cross_section_transformation_derivative_second_variation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_binormal(product_index) * normal_vector(row_component)
                            + mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_binormal(product_index) * normal_vector_variation(row_component * mNumberOfDofs + second_dof_index)
                            + mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation(row_component * mNumberOfDofs + second_dof_index, product_index) * reference_binormal(product_index) * normal_vector_variation(row_component * mNumberOfDofs + dof_index)
                            + cross_section_transformation_derivative(row_component, product_index) * reference_binormal(product_index) * normal_vector_second_variation(row_component * mNumberOfDofs + dof_index, second_dof_index);

                        binormal_torsion_second_variation(dof_index, second_dof_index) += cross_section_transformation_derivative_second_variation(row_component * mNumberOfDofs + dof_index, product_index * mNumberOfDofs + second_dof_index) * reference_normal(product_index) * binormal_vector(row_component)
                            + mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation(row_component * mNumberOfDofs + dof_index, product_index) * reference_normal(product_index) * binormal_vector_variation(row_component * mNumberOfDofs + second_dof_index)
                            + mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation(row_component * mNumberOfDofs + second_dof_index, product_index) * reference_normal(product_index) * binormal_vector_variation(row_component * mNumberOfDofs + dof_index)
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
        noalias(rGTorsion1) = binormal_torsion_second_variation;   
            
        KRATOS_CATCH("")
        
    }
}