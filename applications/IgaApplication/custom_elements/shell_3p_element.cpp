//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//


// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/shell_3p_element.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void Shell3pElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (m_A_ab_covariant_vector.size() != r_number_of_integration_points)
            m_A_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_B_ab_covariant_vector.size() != r_number_of_integration_points)
            m_B_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_dA_vector.size() != r_number_of_integration_points)
            m_dA_vector.resize(r_number_of_integration_points);
        if (m_T_vector.size() != r_number_of_integration_points)
            m_T_vector.resize(r_number_of_integration_points);

        KinematicVariables kinematic_variables(
            GetGeometry().WorkingSpaceDimension());

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
        {
            CalculateKinematics(
                point_number,
                kinematic_variables);

            m_A_ab_covariant_vector[point_number] = kinematic_variables.a_ab_covariant;
            m_B_ab_covariant_vector[point_number] = kinematic_variables.b_ab_covariant;

            m_dA_vector[point_number] = kinematic_variables.dA;

            CalculateTransformation(kinematic_variables, m_T_vector[point_number]);
        }

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void Shell3pElement::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_number_of_integration_points)
            mConstitutiveLawVector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    void Shell3pElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(
                constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);
        }
    }

    ///@}
    ///@name Results on Gauss Points
    ///@{

    void Shell3pElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

        if(rVariable==SHEAR_FORCE_1 || rVariable==SHEAR_FORCE_2)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

                array_1d<double, 2> q = ZeroVector(2);
                CalculateShearForce(point_number, q, rCurrentProcessInfo);

                if (rVariable==SHEAR_FORCE_1)
                {
                    rOutput[point_number] = q[0];
                }
                else if (rVariable==SHEAR_FORCE_2)
                {
                    rOutput[point_number] = q[1];
                }
            }
        }
        else if (rVariable==PK2_STRESS_XX || rVariable==PK2_STRESS_YY || rVariable==PK2_STRESS_XY)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                
                array_1d<double, 3> membrane_stress_pk2_car;
                array_1d<double, 3> bending_stress_pk2_car;

                CalculatePK2Stress(point_number, membrane_stress_pk2_car, bending_stress_pk2_car, rCurrentProcessInfo);

                if (rVariable==PK2_STRESS_XX)
                {
                    rOutput[point_number] = membrane_stress_pk2_car[0];
                }
                else if (rVariable==PK2_STRESS_YY)
                {
                    rOutput[point_number] = membrane_stress_pk2_car[1];
                }
                else if (rVariable==PK2_STRESS_XY)
                {
                    rOutput[point_number] = membrane_stress_pk2_car[2];
                }
            }
        }
        else if(rVariable==CAUCHY_STRESS_XX || rVariable==CAUCHY_STRESS_YY || rVariable==CAUCHY_STRESS_XY 
            || rVariable==CAUCHY_STRESS_TOP_XX || rVariable==CAUCHY_STRESS_TOP_YY || rVariable==CAUCHY_STRESS_TOP_XY  
            || rVariable==CAUCHY_STRESS_BOTTOM_XX || rVariable==CAUCHY_STRESS_BOTTOM_YY || rVariable==CAUCHY_STRESS_BOTTOM_XY
            || rVariable==MEMBRANE_FORCE_XX || rVariable==MEMBRANE_FORCE_YY || rVariable==MEMBRANE_FORCE_XY 
            || rVariable==INTERNAL_MOMENT_XX || rVariable==INTERNAL_MOMENT_YY || rVariable==INTERNAL_MOMENT_XY)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                
                array_1d<double, 3> membrane_stress_cau_car;
                array_1d<double, 3> bending_stress_cau_car;

                CalculateCauchyStress(point_number, membrane_stress_cau_car, bending_stress_cau_car, rCurrentProcessInfo);
                double thickness = this->GetProperties().GetValue(THICKNESS);

                if (rVariable==CAUCHY_STRESS_XX)
                {
                    rOutput[point_number] = membrane_stress_cau_car[0];
                }
                else if (rVariable==CAUCHY_STRESS_YY)
                {
                    rOutput[point_number] = membrane_stress_cau_car[1];
                }
                else if (rVariable==CAUCHY_STRESS_XY)
                {
                    rOutput[point_number] = membrane_stress_cau_car[2];
                }
                else if (rVariable==CAUCHY_STRESS_TOP_XX) 
                {
                    rOutput[point_number] = membrane_stress_cau_car[0] + thickness / 2 * bending_stress_cau_car[0];
                }
                else if (rVariable==CAUCHY_STRESS_TOP_YY) 
                {
                    rOutput[point_number] = membrane_stress_cau_car[1] + thickness / 2 * bending_stress_cau_car[1];
                }
                else if (rVariable==CAUCHY_STRESS_TOP_XY) 
                {
                    rOutput[point_number] = membrane_stress_cau_car[2] + thickness / 2 * bending_stress_cau_car[2];
                }
                else if (rVariable==CAUCHY_STRESS_BOTTOM_XX)
                {
                    rOutput[point_number] = membrane_stress_cau_car[0] - thickness / 2 * bending_stress_cau_car[0];
                }
                else if (rVariable==CAUCHY_STRESS_BOTTOM_YY)
                {
                    rOutput[point_number] = membrane_stress_cau_car[1] - thickness / 2 * bending_stress_cau_car[1];
                }
                else if (rVariable==CAUCHY_STRESS_BOTTOM_XY)
                {
                    rOutput[point_number] = membrane_stress_cau_car[2] - thickness / 2 * bending_stress_cau_car[2];
                }
                else if (rVariable==MEMBRANE_FORCE_XX)
                {
                    rOutput[point_number] = membrane_stress_cau_car[0] * thickness;
                }
                else if (rVariable==MEMBRANE_FORCE_YY)
                {
                    rOutput[point_number] = membrane_stress_cau_car[1] * thickness;
                }
                else if (rVariable==MEMBRANE_FORCE_XY)
                {
                    rOutput[point_number] = membrane_stress_cau_car[2] * thickness;
                }
                else if (rVariable==INTERNAL_MOMENT_XX)
                {
                    rOutput[point_number] = bending_stress_cau_car[0] * pow(thickness, 3) / 12;
                }
                else if (rVariable==INTERNAL_MOMENT_XX)
                {
                    rOutput[point_number] = bending_stress_cau_car[1] * pow(thickness, 3) / 12;
                }
                else if (rVariable==INTERNAL_MOMENT_XY)
                {
                    rOutput[point_number] = bending_stress_cau_car[2] * pow(thickness, 3) / 12;
                }
            } 
        }
        else if (mConstitutiveLawVector[0]->Has(rVariable)) {
            GetValueOnConstitutiveLaw(rVariable, rOutput);
        }
    }

    void Shell3pElement::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 >>& rVariable,
        std::vector<array_1d<double, 3 >>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

        if (rVariable==PK2_STRESS)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

                array_1d<double, 3> membrane_stress_pk2_car;
                array_1d<double, 3> bending_stress_pk2_car;

                CalculatePK2Stress(point_number, membrane_stress_pk2_car, bending_stress_pk2_car, rCurrentProcessInfo);
                rOutput[point_number] = membrane_stress_pk2_car;
            }
        }
        else if(rVariable==CAUCHY_STRESS || rVariable==CAUCHY_STRESS_TOP || rVariable==CAUCHY_STRESS_BOTTOM 
                || rVariable==MEMBRANE_FORCE ||  rVariable==INTERNAL_MOMENT)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

                array_1d<double, 3> membrane_stress_cau_car;
                array_1d<double, 3> bending_stress_cau_car;

                CalculateCauchyStress(point_number, membrane_stress_cau_car, bending_stress_cau_car, rCurrentProcessInfo);
                double thickness = this->GetProperties().GetValue(THICKNESS);

                if (rVariable==CAUCHY_STRESS)
                {
                    rOutput[point_number] = membrane_stress_cau_car;
                }
                else if (rVariable==CAUCHY_STRESS_TOP) 
                {
                    rOutput[point_number] = membrane_stress_cau_car + thickness / 2 * bending_stress_cau_car;
                }
                else if (rVariable==CAUCHY_STRESS_BOTTOM)
                {
                    rOutput[point_number] = membrane_stress_cau_car - thickness / 2 * bending_stress_cau_car;
                }
                else if (rVariable==MEMBRANE_FORCE)
                {
                    rOutput[point_number] = membrane_stress_cau_car * thickness;
                }
                else if (rVariable==INTERNAL_MOMENT)
                {
                    rOutput[point_number] = bending_stress_cau_car * pow(thickness, 3) / 12;
                }
            }
        }
    }


    ///@}
    ///@name Assembly
    ///@{

    void Shell3pElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    ) const
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(
                GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane(3);
            ConstitutiveVariables constitutive_variables_curvature(3);
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_variables_curvature,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            Matrix BMembrane = ZeroMatrix(3, mat_size);
            Matrix BCurvature = ZeroMatrix(3, mat_size);
            CalculateBMembrane(
                point_number,
                BMembrane,
                kinematic_variables);
            CalculateBCurvature(
                point_number,
                BCurvature,
                kinematic_variables);

            // Nonlinear Deformation
            SecondVariations second_variations_strain(mat_size);
            SecondVariations second_variations_curvature(mat_size);
            CalculateSecondVariationStrainCurvature(
                point_number,
                second_variations_strain,
                second_variations_curvature,
                kinematic_variables);

            double integration_weight =
                r_integration_points[point_number].Weight()
                * m_dA_vector[point_number]
                * GetProperties()[THICKNESS];

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                //adding membrane contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BMembrane,
                    constitutive_variables_membrane.ConstitutiveMatrix,
                    integration_weight);
                //adding curvature contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BCurvature,
                    constitutive_variables_curvature.ConstitutiveMatrix,
                    integration_weight);

                // adding  non-linear-contribution to Stiffness-Matrix
                CalculateAndAddNonlinearKm(
                    rLeftHandSideMatrix,
                    second_variations_strain,
                    constitutive_variables_membrane.StressVector,
                    integration_weight);

                CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                    second_variations_curvature,
                    constitutive_variables_curvature.StressVector,
                    integration_weight);
            }
            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
            {
                // operation performed: rRightHandSideVector -= Weight*IntForce
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.StressVector);
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.StressVector);
            }
        }
        KRATOS_CATCH("");
    }

    ///@}
     ///@name Explicit dynamic functions
     ///@{
 
     void Shell3pElement::AddExplicitContribution(
         const VectorType& rRHSVector,
         const Variable<VectorType>& rRHSVariable,
         const Variable<double >& rDestinationVariable,
         const ProcessInfo& rCurrentProcessInfo
         )
     {
         auto& r_geometry = GetGeometry();
         const IndexType nb_nodes = r_geometry.size();
         const IndexType nb_dofs = r_geometry.size() * 3;
 
         if (rDestinationVariable == NODAL_MASS) {
             VectorType element_mass_vector(nb_dofs);
             CalculateLumpedMassVector(element_mass_vector, rCurrentProcessInfo);
 
             for (IndexType i = 0; i < nb_nodes; ++i) {
                 double& r_nodal_mass = r_geometry[i].GetValue(NODAL_MASS);
                 IndexType index = i * 3;
 
                 #pragma omp atomic
                 r_nodal_mass += element_mass_vector(index);
 
             }
         }
     }
 
     void Shell3pElement::AddExplicitContribution(
         const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
         const Variable<array_1d<double, 3>>& rDestinationVariable,
         const ProcessInfo& rCurrentProcessInfo
         )
     {
         auto& r_geometry = GetGeometry();
         const IndexType nb_nodes = r_geometry.size();
         const IndexType nb_dofs = nb_nodes * 3;
 
         #pragma omp critical
         {
             //KRATOS_WATCH("fint")
             //KRATOS_WATCH(rRHSVector)
      
 
         if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
 
             Vector damping_residual_contribution = ZeroVector(nb_dofs);
             Vector current_nodal_velocities = ZeroVector(nb_dofs);
             GetFirstDerivativesVector(current_nodal_velocities, 0);
             Matrix damping_matrix;
             ProcessInfo temp_process_information; // cant pass const ProcessInfo
             CalculateDampingMatrix(damping_matrix, temp_process_information);
             // current residual contribution due to damping
             noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
 
             for (IndexType i = 0; i < nb_nodes; ++i) {
                 IndexType index = 3 * i;
                 array_1d<double, 3>& r_force_residual = r_geometry[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
 
 
                 for (IndexType j = 0; j < 3; ++j) {
                 #pragma omp atomic
                     r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
                 }
             }
         }
         else if (rDestinationVariable == NODAL_INERTIA) {
 
             // Getting the vector mass
             VectorType mass_vector(nb_dofs);
             CalculateLumpedMassVector(mass_vector, rCurrentProcessInfo);
 
             for (IndexType i = 0; i < nb_nodes; ++i) {
                 double& r_nodal_mass = r_geometry[i].GetValue(NODAL_MASS);
                 array_1d<double, 3>& r_nodal_inertia = r_geometry[i].GetValue(NODAL_INERTIA);
                 IndexType index = i * 3;
 
                 #pragma omp atomic
                 r_nodal_mass += mass_vector[index];
 
                 for (IndexType k = 0; k < 3; ++k) {
                     #pragma omp atomic
                     r_nodal_inertia[k] += 0.0;
                 }
             }
         }
         }   
     }

    ///@}
    ///@name Implicit
    ///@{

    void Shell3pElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;
        // Rayleigh Damping Matrix: alpha*M + beta*K

        // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
        double alpha = 0.0;
        if (GetProperties().Has(RAYLEIGH_ALPHA))
            alpha = GetProperties()[RAYLEIGH_ALPHA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA))
            alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

        double beta = 0.0;
        if (GetProperties().Has(RAYLEIGH_BETA))
            beta = GetProperties()[RAYLEIGH_BETA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
            beta = rCurrentProcessInfo[RAYLEIGH_BETA];

        // 2.-Calculate StiffnessMatrix and MassMatrix:
        if (std::abs(alpha) < 1E-12 && std::abs(beta) < 1E-12) {
            // no damping specified, only setting the matrix to zero
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType mat_size = number_of_nodes * 3;
            if (rDampingMatrix.size1() != mat_size || rDampingMatrix.size2() != mat_size) {
                rDampingMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);
        } else if (std::abs(alpha) > 1E-12 && std::abs(beta) < 1E-12) {
            // damping only required with the mass matrix
            CalculateMassMatrix(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= alpha;
        } else if (std::abs(alpha) < 1E-12 && std::abs(beta) > 1E-12) {
            // damping only required with the stiffness matrix
            CalculateLeftHandSide(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= beta;
        } else {
            // damping with both mass matrix and stiffness matrix required
            CalculateLeftHandSide(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= beta;

            Matrix mass_matrix;
            CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);
            noalias(rDampingMatrix) += alpha  * mass_matrix;
        }

        KRATOS_CATCH("")
    }

    void Shell3pElement::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        // Shape function values for all integration points
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            double integration_weight = r_integration_points[point_number].Weight();

            double thickness = this->GetProperties().GetValue(THICKNESS);
            double density = this->GetProperties().GetValue(DENSITY);
            double mass = thickness * density * m_dA_vector[point_number] * integration_weight;

            if (rMassMatrix.size1() != mat_size)
                rMassMatrix.resize(mat_size, mat_size, false);

            rMassMatrix = ZeroMatrix(mat_size, mat_size);

            for (unsigned int r = 0; r<number_of_nodes; r++)
            {
                for (unsigned int s = 0; s<number_of_nodes; s++)
                {
                    rMassMatrix(3 * s, 3 * r) = r_N(point_number, s)*r_N(point_number, r) * mass;
                    rMassMatrix(3 * s + 1, 3 * r + 1) = rMassMatrix(3 * s, 3 * r);
                    rMassMatrix(3 * s + 2, 3 * r + 2) = rMassMatrix(3 * s, 3 * r);
                }
            }
        }
        KRATOS_CATCH("")
    }


     void Shell3pElement::CalculateLumpedMassVector(
         Vector& rLumpedMassVector,
         const ProcessInfo& rCurrentProcessInfo
         ) const
     {
         const auto& r_geometry = GetGeometry();
         const IndexType nb_nodes = r_geometry.size();
         const Matrix& r_N = r_geometry.ShapeFunctionsValues();
         auto& r_integration_points = r_geometry.IntegrationPoints();
         const double num_integration_points = r_integration_points.size();
         // Clear matrix
         if (rLumpedMassVector.size() != nb_nodes * 3) {
             rLumpedMassVector.resize(nb_nodes * 3, false);
         }
 
         for (IndexType point_number = 0; point_number < num_integration_points; ++point_number) {
             
              double integration_weight = r_integration_points[point_number].Weight();
 
             double thickness = this->GetProperties().GetValue(THICKNESS);
             double density = this->GetProperties().GetValue(DENSITY);
             double mass = thickness * density * m_dA_vector[point_number] * integration_weight;
 
             // #pragma omp critical
             // {
             //     KRATOS_WATCH(mass)
             // }
 
             for (IndexType i = 0; i < nb_nodes; ++i) {
                 for (IndexType j = 0; j < 3; ++j) {
                     IndexType index = i * 3 + j;
 
                     rLumpedMassVector[index] = mass * r_N(point_number, i);
                 }
             }
         }
         #pragma omp critical
         {
             //KRATOS_WATCH(rLumpedMassVector)
         }
     }

    ///@}
    ///@name Kinematics
    ///@{

    void Shell3pElement::CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables
    ) const
    {
        Matrix J;
        GetGeometry().Jacobian(J, IntegrationPointIndex);

        rKinematicVariables.a1 = column(J, 0);
        rKinematicVariables.a2 = column(J, 1);

        //not-normalized base vector 3
        MathUtils<double>::CrossProduct(rKinematicVariables.a3_tilde, rKinematicVariables.a1, rKinematicVariables.a2);

        //differential area dA
        rKinematicVariables.dA = norm_2(rKinematicVariables.a3_tilde);

        //base vector 3 normalized
        noalias(rKinematicVariables.a3) = rKinematicVariables.a3_tilde / rKinematicVariables.dA;

        //GetCovariantMetric
        rKinematicVariables.a_ab_covariant[0] = pow(rKinematicVariables.a1[0], 2) + pow(rKinematicVariables.a1[1], 2) + pow(rKinematicVariables.a1[2], 2);
        rKinematicVariables.a_ab_covariant[1] = pow(rKinematicVariables.a2[0], 2) + pow(rKinematicVariables.a2[1], 2) + pow(rKinematicVariables.a2[2], 2);
        rKinematicVariables.a_ab_covariant[2] = rKinematicVariables.a1[0] * rKinematicVariables.a2[0] + rKinematicVariables.a1[1] * rKinematicVariables.a2[1] + rKinematicVariables.a1[2] * rKinematicVariables.a2[2];

        Matrix H = ZeroMatrix(3, 3);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        rKinematicVariables.b_ab_covariant[0] = H(0, 0) * rKinematicVariables.a3[0] + H(1, 0) * rKinematicVariables.a3[1] + H(2, 0) * rKinematicVariables.a3[2];
        rKinematicVariables.b_ab_covariant[1] = H(0, 1) * rKinematicVariables.a3[0] + H(1, 1) * rKinematicVariables.a3[1] + H(2, 1) * rKinematicVariables.a3[2];
        rKinematicVariables.b_ab_covariant[2] = H(0, 2) * rKinematicVariables.a3[0] + H(1, 2) * rKinematicVariables.a3[1] + H(2, 2) * rKinematicVariables.a3[2];
    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void Shell3pElement::CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& rT
    ) const
    {
        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rKinematicVariables.a_ab_covariant[0] * rKinematicVariables.a_ab_covariant[1]
                - rKinematicVariables.a_ab_covariant[2] * rKinematicVariables.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] = inv_det_g_ab * rKinematicVariables.a_ab_covariant[1];
        a_ab_contravariant[1] = inv_det_g_ab * rKinematicVariables.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rKinematicVariables.a_ab_covariant[2];

        //Contravariant base vectors
        array_1d<double, 3> a_contravariant_1 = rKinematicVariables.a1 * a_ab_contravariant[0] + rKinematicVariables.a2 * a_ab_contravariant[2];
        array_1d<double, 3> a_contravariant_2 = rKinematicVariables.a1 * a_ab_contravariant[2] + rKinematicVariables.a2 * a_ab_contravariant[1];


        //Local cartesian coordinates
        double l_a1 = norm_2(rKinematicVariables.a1);
        array_1d<double, 3> e1 = rKinematicVariables.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        // e * a_contravariant
        Matrix G = ZeroMatrix(2, 2);
        G(0, 0) = inner_prod(e1, a_contravariant_1);
        G(0, 1) = inner_prod(e1, a_contravariant_2);
        G(1, 0) = inner_prod(e2, a_contravariant_1);
        G(1, 1) = inner_prod(e2, a_contravariant_2);

        //Transformation matrix T
        if (rT.size1() != 3 && rT.size2() != 3)
            rT.resize(3, 3);
        noalias(rT) = ZeroMatrix(3, 3);

        rT(0, 0) = pow(G(0, 0), 2);
        rT(0, 1) = pow(G(0, 1), 2);
        rT(0, 2) = 2 * G(0, 0) * G(0, 1);

        rT(1, 0) = pow(G(1, 0), 2);
        rT(1, 1) = pow(G(1, 1), 2);
        rT(1, 2) = 2 * G(1, 0) * G(1, 1);

        rT(2, 0) = 2 * G(0, 0) * G(1, 0);
        rT(2, 1) = 2 * G(0, 1) * G(1, 1);
        rT(2, 2) = 2 * (G(0, 0) * G(1, 1) + G(0, 1) * G(1, 0));
    }

    /* Computes the transformation matrix T from the local cartesian basis to
    *  the local cartesian basis.
    */
    void Shell3pElement::CalculateTransformationFromCovariantToCartesian(
        const KinematicVariables& rKinematicVariables,
        Matrix& rTCovToCar
    ) const
    {
        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rKinematicVariables.a_ab_covariant[0] * rKinematicVariables.a_ab_covariant[1]
                - rKinematicVariables.a_ab_covariant[2] * rKinematicVariables.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] =  inv_det_g_ab * rKinematicVariables.a_ab_covariant[1];
        a_ab_contravariant[1] =  inv_det_g_ab * rKinematicVariables.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rKinematicVariables.a_ab_covariant[2];

        //Contravariant base vectors
        array_1d<double, 3> a_contravariant_2 = rKinematicVariables.a1*a_ab_contravariant[2] + rKinematicVariables.a2*a_ab_contravariant[1];

        //Local cartesian coordinates
        double l_a1 = norm_2(rKinematicVariables.a1);
        array_1d<double, 3> e1 = rKinematicVariables.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        // e * a_covariant
        double G_00 = inner_prod(e1, rKinematicVariables.a1);
        double G_01 = inner_prod(e1, rKinematicVariables.a2);
        double G_11 = inner_prod(e2, rKinematicVariables.a2);

        //Transformation matrix T from covariant to local cartesian coordinate system
        rTCovToCar(0, 0) = pow(G_00, 2);
        rTCovToCar(0, 1) = pow(G_01, 2);
        rTCovToCar(0, 2) = 2 * G_00 * G_01;
        rTCovToCar(1, 1) = pow(G_11, 2);
        rTCovToCar(2, 1) = G_01 * G_11;
        rTCovToCar(2, 2) = G_00 * G_11;
    }

    void Shell3pElement::CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector[IntegrationPointIndex]);
        noalias(rThisConstitutiveVariablesMembrane.StrainVector) = prod(m_T_vector[IntegrationPointIndex], strain_vector);

        array_1d<double, 3> curvature_vector = rActualKinematic.b_ab_covariant - m_B_ab_covariant_vector[IntegrationPointIndex];
        noalias(rThisConstitutiveVariablesCurvature.StrainVector) = prod(m_T_vector[IntegrationPointIndex], curvature_vector);

        // Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector); //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

        mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);

        double thickness = this->GetProperties().GetValue(THICKNESS);
        noalias(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix) = rThisConstitutiveVariablesMembrane.ConstitutiveMatrix * (pow(thickness, 2) / 12);

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
        noalias(rThisConstitutiveVariablesCurvature.StressVector) = prod(
            trans(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix), rThisConstitutiveVariablesCurvature.StrainVector);
    }

    void Shell3pElement::CalculateBMembrane(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        const KinematicVariables& rActualKinematic) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        if (rB.size1() != 3 || rB.size2() != mat_size)
            rB.resize(3, mat_size);
        noalias(rB) = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            array_1d<double, 3> dE_curvilinear;
            // strain
            dE_curvilinear[0] = r_DN_De(kr, 0) * rActualKinematic.a1(dirr);
            dE_curvilinear[1] = r_DN_De(kr, 1) * rActualKinematic.a2(dirr);
            dE_curvilinear[2] = 0.5 * (r_DN_De(kr, 0) * rActualKinematic.a2(dirr) + rActualKinematic.a1(dirr) * r_DN_De(kr, 1));

            rB(0, r) = m_T_vector[IntegrationPointIndex](0, 0) * dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](0, 1) * dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](0, 2) * dE_curvilinear[2];
            rB(1, r) = m_T_vector[IntegrationPointIndex](1, 0) * dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](1, 1) * dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](1, 2) * dE_curvilinear[2];
            rB(2, r) = m_T_vector[IntegrationPointIndex](2, 0) * dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](2, 1) * dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](2, 2) * dE_curvilinear[2];
        }
    }

    void Shell3pElement::CalculateBCurvature(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        const KinematicVariables& rActualKinematic) const
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        Matrix da3 = ZeroMatrix(3, 3);
        Matrix dn = ZeroMatrix(3, 3);
        Matrix b = ZeroMatrix(3, mat_size);

        double inv_dA = 1 / rActualKinematic.dA;
        double inv_dA3 = 1 / std::pow(rActualKinematic.dA, 3);

        Matrix H = ZeroMatrix(3, 3);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex));

        for (IndexType i = 0; i < number_of_control_points; i++)
        {
            IndexType index = 3 * i;
            //first line
            da3(0, 0) = 0;
            da3(0, 1) = -r_DN_De(i, 0) * rActualKinematic.a2[2] + r_DN_De(i, 1) * rActualKinematic.a1[2];
            da3(0, 2) = r_DN_De(i, 0) * rActualKinematic.a2[1] - r_DN_De(i, 1) * rActualKinematic.a1[1];

            //second line
            da3(1, 0) = r_DN_De(i, 0) * rActualKinematic.a2[2] - r_DN_De(i, 1) * rActualKinematic.a1[2];
            da3(1, 1) = 0;
            da3(1, 2) = -r_DN_De(i, 0) * rActualKinematic.a2[0] + r_DN_De(i, 1) * rActualKinematic.a1[0];

            //third line
            da3(2, 0) = -r_DN_De(i, 0) * rActualKinematic.a2[1] + r_DN_De(i, 1) * rActualKinematic.a1[1];
            da3(2, 1) = r_DN_De(i, 0) * rActualKinematic.a2[0] - r_DN_De(i, 1) * rActualKinematic.a1[0];
            da3(2, 2) = 0;

            for (IndexType j = 0; j < 3; j++)
            {
                double a3da3la3 = (rActualKinematic.a3_tilde[0] * da3(j, 0) + rActualKinematic.a3_tilde[1] * da3(j, 1) + rActualKinematic.a3_tilde[2] * da3(j, 2)) * inv_dA3;

                dn(j, 0) = da3(j, 0) * inv_dA - rActualKinematic.a3_tilde[0] * a3da3la3;
                dn(j, 1) = da3(j, 1) * inv_dA - rActualKinematic.a3_tilde[1] * a3da3la3;
                dn(j, 2) = da3(j, 2) * inv_dA - rActualKinematic.a3_tilde[2] * a3da3la3;
            }

            // curvature vector [K11,K22,K12] referred to curvilinear coordinate system
            b(0, index) = 0 - (r_DDN_DDe(i, 0) * rActualKinematic.a3[0] + H(0, 0) * dn(0, 0) + H(1, 0) * dn(0, 1) + H(2, 0) * dn(0, 2));
            b(0, index + 1) = 0 - (r_DDN_DDe(i, 0) * rActualKinematic.a3[1] + H(0, 0) * dn(1, 0) + H(1, 0) * dn(1, 1) + H(2, 0) * dn(1, 2));
            b(0, index + 2) = 0 - (r_DDN_DDe(i, 0) * rActualKinematic.a3[2] + H(0, 0) * dn(2, 0) + H(1, 0) * dn(2, 1) + H(2, 0) * dn(2, 2));

            //second line
            b(1, index) = 0 - (r_DDN_DDe(i, 2) * rActualKinematic.a3[0] + H(0, 1) * dn(0, 0) + H(1, 1) * dn(0, 1) + H(2, 1) * dn(0, 2));
            b(1, index + 1) = 0 - (r_DDN_DDe(i, 2) * rActualKinematic.a3[1] + H(0, 1) * dn(1, 0) + H(1, 1) * dn(1, 1) + H(2, 1) * dn(1, 2));
            b(1, index + 2) = 0 - (r_DDN_DDe(i, 2) * rActualKinematic.a3[2] + H(0, 1) * dn(2, 0) + H(1, 1) * dn(2, 1) + H(2, 1) * dn(2, 2));

            //third line
            b(2, index) = 0 - (r_DDN_DDe(i, 1) * rActualKinematic.a3[0] + H(0, 2) * dn(0, 0) + H(1, 2) * dn(0, 1) + H(2, 2) * dn(0, 2));
            b(2, index + 1) = 0 - (r_DDN_DDe(i, 1) * rActualKinematic.a3[1] + H(0, 2) * dn(1, 0) + H(1, 2) * dn(1, 1) + H(2, 2) * dn(1, 2));
            b(2, index + 2) = 0 - (r_DDN_DDe(i, 1) * rActualKinematic.a3[2] + H(0, 2) * dn(2, 0) + H(1, 2) * dn(2, 1) + H(2, 2) * dn(2, 2));
        }

        noalias(rB) = -prod(m_T_vector[IntegrationPointIndex], b);

        KRATOS_CATCH("")
    }

    void Shell3pElement::CalculateSecondVariationStrainCurvature(
        const IndexType IntegrationPointIndex,
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const KinematicVariables& rActualKinematic) const
    {
        const auto& r_geometry = GetGeometry();

        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        double l_a3 = norm_2(rActualKinematic.a3_tilde);
        double l_a3_3 = pow(l_a3, 3);
        double l_a3_5 = pow(l_a3, 5);
        double inv_l_a3 = 1 / l_a3;
        double inv_l_a3_3 = 1 / l_a3_3;
        double inv_l_a3_5 = 1 / l_a3_5;

        Matrix S_da3 = ZeroMatrix(3, mat_size);
        Vector S_a3_da3 = ZeroVector(mat_size);
        Vector S_a3_da3_l_a3_3 = ZeroVector(mat_size);
        Matrix S_dn = ZeroMatrix(3, mat_size);

        Matrix H = ZeroMatrix(3, 3);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        // first variation of strain and curvature w.r.t. dof
        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            array_1d<double, 3> S_dg_1 = ZeroVector(3);
            array_1d<double, 3> S_dg_2 = ZeroVector(3);
            S_dg_1(dirr) = r_DN_De(kr, 0);
            S_dg_2(dirr) = r_DN_De(kr, 1);

            // curvature
            S_da3(0, r) = S_dg_1(1) * rActualKinematic.a2(2) - S_dg_1(2) * rActualKinematic.a2(1) + rActualKinematic.a1(1) * S_dg_2(2) - rActualKinematic.a1(2) * S_dg_2(1);
            S_da3(1, r) = S_dg_1(2) * rActualKinematic.a2(0) - S_dg_1(0) * rActualKinematic.a2(2) + rActualKinematic.a1(2) * S_dg_2(0) - rActualKinematic.a1(0) * S_dg_2(2);
            S_da3(2, r) = S_dg_1(0) * rActualKinematic.a2(1) - S_dg_1(1) * rActualKinematic.a2(0) + rActualKinematic.a1(0) * S_dg_2(1) - rActualKinematic.a1(1) * S_dg_2(0);

            S_a3_da3[r] = rActualKinematic.a3_tilde[0] * S_da3(0, r) + rActualKinematic.a3_tilde[1] * S_da3(1, r) + rActualKinematic.a3_tilde[2] * S_da3(2, r);
            S_a3_da3_l_a3_3[r] = S_a3_da3[r] * inv_l_a3_3;

            S_dn(0, r) = S_da3(0, r) * inv_l_a3 - rActualKinematic.a3_tilde[0] * S_a3_da3_l_a3_3[r];
            S_dn(1, r) = S_da3(1, r) * inv_l_a3 - rActualKinematic.a3_tilde[1] * S_a3_da3_l_a3_3[r];
            S_dn(2, r) = S_da3(2, r) * inv_l_a3 - rActualKinematic.a3_tilde[2] * S_a3_da3_l_a3_3[r];
        }

        // second variation of strain and curvature w.r.t. dofs
        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            for (IndexType s = 0; s <= r; s++)
            {
                // local node number ks and dof direction dirs
                IndexType ks = s / 3;
                IndexType dirs = s % 3;

                // strain
                array_1d<double, 3> ddE_cu = ZeroVector(3);
                if (dirr == dirs)
                {
                    ddE_cu[0] = r_DN_De(kr, 0) * r_DN_De(ks, 0);
                    ddE_cu[1] = r_DN_De(kr, 1) * r_DN_De(ks, 1);
                    ddE_cu[2] = 0.5 * (r_DN_De(kr, 0) * r_DN_De(ks, 1) + r_DN_De(kr, 1) * r_DN_De(ks, 0));

                    rSecondVariationsStrain.B11(r, s) = m_T_vector[IntegrationPointIndex](0, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](0, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](0, 2) * ddE_cu[2];
                    rSecondVariationsStrain.B22(r, s) = m_T_vector[IntegrationPointIndex](1, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](1, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](1, 2) * ddE_cu[2];
                    rSecondVariationsStrain.B12(r, s) = m_T_vector[IntegrationPointIndex](2, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](2, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](2, 2) * ddE_cu[2];
                }

                // curvature
                array_1d<double, 3> dda3 = ZeroVector(3);
                int dirt = 4 - static_cast<int>(dirr) - static_cast<int>(dirs);
                int ddir = static_cast<int>(dirr) - static_cast<int>(dirs);
                if (ddir == -1) dda3(dirt - 1) = r_DN_De(kr, 0) * r_DN_De(ks, 1) - r_DN_De(ks, 0) * r_DN_De(kr, 1);
                else if (ddir == 2) dda3(dirt - 1) = r_DN_De(kr, 0) * r_DN_De(ks, 1) - r_DN_De(ks, 0) * r_DN_De(kr, 1);
                else if (ddir == 1) dda3(dirt - 1) = -r_DN_De(kr, 0) * r_DN_De(ks, 1) + r_DN_De(ks, 0) * r_DN_De(kr, 1);
                else if (ddir == -2) dda3(dirt - 1) = -r_DN_De(kr, 0) * r_DN_De(ks, 1) + r_DN_De(ks, 0) * r_DN_De(kr, 1);

                double c = -(dda3[0] * rActualKinematic.a3_tilde[0] + dda3[1] * rActualKinematic.a3_tilde[1] + dda3[2] * rActualKinematic.a3_tilde[2]
                    + S_da3(0, r) * S_da3(0, s) + S_da3(1, r) * S_da3(1, s) + S_da3(2, r) * S_da3(2, s)
                    ) * inv_l_a3_3;

                double d = 3.0 * S_a3_da3[r] * S_a3_da3[s] * inv_l_a3_5;

                array_1d<double, 3> ddn = ZeroVector(3);
                ddn[0] = dda3[0] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(0, r) - S_a3_da3_l_a3_3[r] * S_da3(0, s) + (c + d) * rActualKinematic.a3_tilde[0];
                ddn[1] = dda3[1] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(1, r) - S_a3_da3_l_a3_3[r] * S_da3(1, s) + (c + d) * rActualKinematic.a3_tilde[1];
                ddn[2] = dda3[2] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(2, r) - S_a3_da3_l_a3_3[r] * S_da3(2, s) + (c + d) * rActualKinematic.a3_tilde[2];

                array_1d<double, 3> ddK_cu = ZeroVector(3);
                ddK_cu[0] = r_DDN_DDe(kr, 0) * S_dn(dirr, s) + r_DDN_DDe(ks, 0) * S_dn(dirs, r)
                    + H(0, 0) * ddn[0] + H(1, 0) * ddn[1] + H(2, 0) * ddn[2];
                ddK_cu[1] = r_DDN_DDe(kr, 2) * S_dn(dirr, s) + r_DDN_DDe(ks, 2) * S_dn(dirs, r)
                    + H(0, 1) * ddn[0] + H(1, 1) * ddn[1] + H(2, 1) * ddn[2];
                ddK_cu[2] = r_DDN_DDe(kr, 1) * S_dn(dirr, s) + r_DDN_DDe(ks, 1) * S_dn(dirs, r)
                    + H(0, 2) * ddn[0] + H(1, 2) * ddn[1] + H(2, 2) * ddn[2];

                rSecondVariationsCurvature.B11(r, s) = m_T_vector[IntegrationPointIndex](0, 0) * ddK_cu[0] + m_T_vector[IntegrationPointIndex](0, 1) * ddK_cu[1] + m_T_vector[IntegrationPointIndex](0, 2) * ddK_cu[2];
                rSecondVariationsCurvature.B11(s, r) = rSecondVariationsCurvature.B11(r, s);
                rSecondVariationsCurvature.B22(r, s) = m_T_vector[IntegrationPointIndex](1, 0) * ddK_cu[0] + m_T_vector[IntegrationPointIndex](1, 1) * ddK_cu[1] + m_T_vector[IntegrationPointIndex](1, 2) * ddK_cu[2];
                rSecondVariationsCurvature.B22(s, r) = rSecondVariationsCurvature.B22(r, s);
                rSecondVariationsCurvature.B12(r, s) = m_T_vector[IntegrationPointIndex](2, 0) * ddK_cu[0] + m_T_vector[IntegrationPointIndex](2, 1) * ddK_cu[1] + m_T_vector[IntegrationPointIndex](2, 2) * ddK_cu[2];
                rSecondVariationsCurvature.B12(s, r) = rSecondVariationsCurvature.B12(r, s);
            }
        }
    }

    ///@}
    ///@name Stiffness matrix assembly
    ///@{

    inline void Shell3pElement::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rB,
        const Matrix& rD,
        const double IntegrationWeight
    ) const
    {
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(rB), Matrix(prod(rD, rB)));
    }

    inline void Shell3pElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& rSecondVariationsStrain,
        const Vector& rSD,
        const double IntegrationWeight) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        for (IndexType n = 0; n < mat_size; n++)
        {
            for (IndexType m = 0; m <= n; m++)
            {
                double nm = (rSD[0] * rSecondVariationsStrain.B11(n, m)
                    + rSD[1] * rSecondVariationsStrain.B22(n, m)
                    + rSD[2] * rSecondVariationsStrain.B12(n, m)) * IntegrationWeight;

                rLeftHandSideMatrix(n, m) += nm;
                if (n != m)
                    rLeftHandSideMatrix(m, n) += nm;
            }
        }
    }

    ///@}
    ///@name Stress recovery
    ///@{
    
    void Shell3pElement::CalculatePK2Stress(
        const IndexType IntegrationPointIndex,
        array_1d<double, 3>& rPK2MembraneStressCartesian,
        array_1d<double, 3>& rPK2BendingStressCartesian,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        // Compute Kinematics and Metric
        KinematicVariables kinematic_variables(
            GetGeometry().WorkingSpaceDimension());
        CalculateKinematics(
            IntegrationPointIndex,
            kinematic_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        ConstitutiveVariables constitutive_variables_membrane(3);
        ConstitutiveVariables constitutive_variables_curvature(3);
        CalculateConstitutiveVariables(
            IntegrationPointIndex,
            kinematic_variables,
            constitutive_variables_membrane,
            constitutive_variables_curvature,
            constitutive_law_parameters,
            ConstitutiveLaw::StressMeasure_PK2);
        
        double thickness = this->GetProperties().GetValue(THICKNESS);

        rPK2MembraneStressCartesian = constitutive_variables_membrane.StressVector;
        rPK2BendingStressCartesian = -1.0 * constitutive_variables_curvature.StressVector / pow(thickness, 2) * 12;
    }
    
    void Shell3pElement::CalculateCauchyStress(
        const IndexType IntegrationPointIndex,
        array_1d<double, 3>& rCauchyMembraneStressesCartesian, 
        array_1d<double, 3>& rCauchyBendingStressesCartesian, 
        const ProcessInfo& rCurrentProcessInfo) const
    {
        // Compute PK2 stress
        array_1d<double, 3> membrane_stress_pk2_car;
        array_1d<double, 3> bending_stress_pk2_car;

        CalculatePK2Stress(IntegrationPointIndex, membrane_stress_pk2_car, bending_stress_pk2_car, rCurrentProcessInfo);

        // Compute Kinematics and Metric
        KinematicVariables kinematic_variables(
            GetGeometry().WorkingSpaceDimension());
        CalculateKinematics(
            IntegrationPointIndex,
            kinematic_variables);

        double detF = kinematic_variables.dA / m_dA_vector[IntegrationPointIndex];

        // Compute Transformation matrix T from local cartesian coordinate system to covariant basis
        Matrix T_car_to_cov = ZeroMatrix(3,3);
        T_car_to_cov = trans(m_T_vector[IntegrationPointIndex]);
        for (IndexType i = 0; i < 3; i++)
        {
            T_car_to_cov(2, i) = T_car_to_cov(i, 2) / 2;
        }

        // Compute Transformation matrix T from covariant basis to local cartesian coordinate system            
        Matrix T_cov_to_car = ZeroMatrix(3,3);
        CalculateTransformationFromCovariantToCartesian(kinematic_variables, T_cov_to_car);

        // Compute Cauchy stress
        array_1d<double, 3> membrane_stress_pk2_cov = prod(T_car_to_cov, membrane_stress_pk2_car);
        array_1d<double, 3> membrane_stress_cau_cov = membrane_stress_pk2_cov / detF;
        array_1d<double, 3> membrane_stress_cau_car = prod(T_cov_to_car, membrane_stress_cau_cov);

        array_1d<double, 3> bending_stress_pk2_cov = prod(T_car_to_cov, bending_stress_pk2_car);
        array_1d<double, 3> bending_stress_cau_cov = bending_stress_pk2_cov / detF;
        array_1d<double, 3> bending_stress_cau_car = prod(T_cov_to_car, bending_stress_cau_cov);

        rCauchyMembraneStressesCartesian = membrane_stress_cau_car;
        rCauchyBendingStressesCartesian = bending_stress_cau_car;
    }

    void Shell3pElement::CalculateShearForce(
        const IndexType IntegrationPointIndex,
        array_1d<double, 2>& rq, 
        const ProcessInfo& rCurrentProcessInfo) const
    {
        // Compute Kinematics and Metric
        KinematicVariables kinematic_variables(
            GetGeometry().WorkingSpaceDimension());
        CalculateKinematics(
            IntegrationPointIndex,
            kinematic_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        ConstitutiveVariables constitutive_variables_membrane(3);
        ConstitutiveVariables constitutive_variables_curvature(3);
        CalculateConstitutiveVariables(
            IntegrationPointIndex,
            kinematic_variables,
            constitutive_variables_membrane,
            constitutive_variables_curvature,
            constitutive_law_parameters,
            ConstitutiveLaw::StressMeasure_PK2);
        
        double thickness = this->GetProperties().GetValue(THICKNESS);

        // Calculate Hessian Matrix at initial configuration
        const SizeType number_of_points = GetGeometry().size();
        const SizeType working_space_dimension = 3;
        const Matrix& rDDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());
        
        Matrix H_initial;
        H_initial.resize(working_space_dimension, working_space_dimension);
        H_initial = ZeroMatrix(working_space_dimension, working_space_dimension);

        for (IndexType k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3>& coords = GetGeometry()[k].GetInitialPosition();

            H_initial(0, 0) += rDDN_DDe(k, 0)*coords[0];
            H_initial(0, 1) += rDDN_DDe(k, 2)*coords[0];
            H_initial(0, 2) += rDDN_DDe(k, 1)*coords[0];

            H_initial(1, 0) += rDDN_DDe(k, 0)*coords[1];
            H_initial(1, 1) += rDDN_DDe(k, 2)*coords[1];
            H_initial(1, 2) += rDDN_DDe(k, 1)*coords[1];

            H_initial(2, 0) += rDDN_DDe(k, 0)*coords[2];
            H_initial(2, 1) += rDDN_DDe(k, 2)*coords[2];
            H_initial(2, 2) += rDDN_DDe(k, 1)*coords[2];
        }

        // Calculate Hessian Matrix at actual configuration
        Matrix H_actual;
        H_actual.resize(working_space_dimension, working_space_dimension);
        H_actual = ZeroMatrix(working_space_dimension, working_space_dimension);
        CalculateHessian(H_actual, rDDN_DDe);

        // Calculate derivative of curvature w.r.t. theta1 and theta2 at initial configuration
        array_1d<double, 3> DCurvature_D1_initial = ZeroVector(3);
        array_1d<double, 3> DCurvature_D2_initial = ZeroVector(3);
        CalculateDerivativeOfCurvatureInitial(IntegrationPointIndex, DCurvature_D1_initial, DCurvature_D2_initial, H_initial);
        
        // Calculate derivative of curvature w.r.t. theta1 and theta2 at actual configuration
        array_1d<double, 3> DCurvature_D1_actual = ZeroVector(3);
        array_1d<double, 3> DCurvature_D2_actual = ZeroVector(3);
        CalculateDerivativeOfCurvatureActual(IntegrationPointIndex, DCurvature_D1_actual, DCurvature_D2_actual, H_actual, kinematic_variables);

        std::vector<array_1d<double, 3>> Dk_con_Dalpha(2);
        Dk_con_Dalpha[0] = DCurvature_D1_initial - DCurvature_D1_actual;
        Dk_con_Dalpha[1] = DCurvature_D2_initial - DCurvature_D2_actual;

        array_1d<double, 3> k_con = -1.0 * (kinematic_variables.b_ab_covariant - m_B_ab_covariant_vector[IntegrationPointIndex]);
        array_1d<double, 3> k_car = prod(m_T_vector[IntegrationPointIndex], k_con);
        array_1d<double, 3> m_car = prod(constitutive_variables_curvature.ConstitutiveMatrix, k_car)*thickness;

        // derivative of the transformation matrix T_con_to_car (contravariant to local Cartesian basis)
        std::vector<Matrix> DT_con_to_car_init_Dalpha(2, ZeroMatrix(3, 3));
        std::vector<Matrix> DT_car_to_cov_init_Dalpha(2, ZeroMatrix(3, 3));
        CalculateDerivativeTransformationMatrices(IntegrationPointIndex, DT_con_to_car_init_Dalpha, DT_car_to_cov_init_Dalpha, H_initial);

        // derivative of the moment
        std::vector<array_1d<double, 3>> Dm_car_Dalpha(2);
        std::vector<array_1d<double, 3>> Dm_cov_Dalpha(2);
        array_1d<double, 3> Dk_car_D1 = prod(m_T_vector[IntegrationPointIndex], Dk_con_Dalpha[0]) + prod(DT_con_to_car_init_Dalpha[0], k_con);
        array_1d<double, 3> Dk_car_D2 = prod(m_T_vector[IntegrationPointIndex], Dk_con_Dalpha[1]) + prod(DT_con_to_car_init_Dalpha[1], k_con);
        Dm_car_Dalpha[0] = prod(constitutive_variables_curvature.ConstitutiveMatrix, Dk_car_D1)*thickness;
        Dm_car_Dalpha[1] = prod(constitutive_variables_curvature.ConstitutiveMatrix, Dk_car_D2)*thickness;

        //Transformation matrix T from local cartesian coordinate system to covariant basis
        Matrix T_car_to_cov = ZeroMatrix(3,3);
        T_car_to_cov = trans(m_T_vector[IntegrationPointIndex]);
        for (IndexType i = 0; i < 3; i++)
        {
            T_car_to_cov(2, i) = T_car_to_cov(i, 2) / 2;
        }

        Dm_cov_Dalpha[0] = prod(T_car_to_cov, Dm_car_Dalpha[0]) + prod(DT_car_to_cov_init_Dalpha[0], m_car);
        Dm_cov_Dalpha[1] = prod(T_car_to_cov, Dm_car_Dalpha[1]) + prod(DT_car_to_cov_init_Dalpha[1], m_car);

        array_1d<double, 2> q_pk2_cov;
        q_pk2_cov[0] = Dm_cov_Dalpha[0](0) / std::sqrt(kinematic_variables.a_ab_covariant[0]) + Dm_cov_Dalpha[1](2) / std::sqrt(kinematic_variables.a_ab_covariant[1]);
        q_pk2_cov[1] = Dm_cov_Dalpha[1](1) / std::sqrt(kinematic_variables.a_ab_covariant[1]) + Dm_cov_Dalpha[0](2) / std::sqrt(kinematic_variables.a_ab_covariant[0]);
        double detF = kinematic_variables.dA / m_dA_vector[IntegrationPointIndex];
        array_1d<double, 2> q_cau_cov = q_pk2_cov / detF;

        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (kinematic_variables.a_ab_covariant[0] * kinematic_variables.a_ab_covariant[1]
                - kinematic_variables.a_ab_covariant[2] * kinematic_variables.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] =  inv_det_g_ab * kinematic_variables.a_ab_covariant[1];
        a_ab_contravariant[1] =  inv_det_g_ab * kinematic_variables.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * kinematic_variables.a_ab_covariant[2];

        //Contravariant base vectors
        array_1d<double, 3> a_contravariant_2 = kinematic_variables.a1*a_ab_contravariant[2] + kinematic_variables.a2*a_ab_contravariant[1];

        //Local cartesian coordinates
        double l_a1 = norm_2(kinematic_variables.a1);
        array_1d<double, 3> e1 = kinematic_variables.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        rq[0] = inner_prod(e1, kinematic_variables.a1) * q_cau_cov[0];
        rq[1] = inner_prod(e2, kinematic_variables.a2) * q_cau_cov[1];
    }

    void Shell3pElement::CalculateDerivativeOfCurvatureInitial(
        const IndexType IntegrationPointIndex,
        array_1d<double, 3>& rDCurvature_D1,
        array_1d<double, 3>& rDCurvature_D2,
        const Matrix& rHessian) const
    {
        // base vector at initial configuration
        const SizeType number_of_nodes = GetGeometry().size();
        const Matrix& rDN_De = GetGeometry().ShapeFunctionDerivatives(1, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        array_1d<double, 3> a1 = ZeroVector(3);
        array_1d<double, 3> a2 = ZeroVector(3);
        array_1d<double, 3> a3_tilde = ZeroVector(3);
        array_1d<double, 3> a3 = ZeroVector(3);

        for (SizeType i = 0; i < number_of_nodes; i++)
        {
            a1[0] += GetGeometry().GetPoint( i ).X0() * rDN_De(i, 0);
            a1[1] += GetGeometry().GetPoint( i ).Y0() * rDN_De(i, 0);
            a1[2] += GetGeometry().GetPoint( i ).Z0() * rDN_De(i, 0);

            a2[0] += GetGeometry().GetPoint( i ).X0() * rDN_De(i, 1);
            a2[1] += GetGeometry().GetPoint( i ).Y0() * rDN_De(i, 1);
            a2[2] += GetGeometry().GetPoint( i ).Z0() * rDN_De(i, 1);
        }

        MathUtils<double>::CrossProduct(a3_tilde, a1, a2);
        noalias(a3) = a3_tilde / norm_2(a3_tilde);

        // second derivatives of the base vectors w.r.t. the curvilinear coords at initial config.
        const Matrix& rDDDN_DDDe = GetGeometry().ShapeFunctionDerivatives(3, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        array_1d<double, 3> DDa1_DD11 = ZeroVector(3);
        array_1d<double, 3> DDa1_DD12 = ZeroVector(3);
        array_1d<double, 3> DDa2_DD21 = ZeroVector(3);
        array_1d<double, 3> DDa2_DD22 = ZeroVector(3);

        const SizeType number_of_points = GetGeometry().size();
    
        for (IndexType k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3> coords = GetGeometry()[k].GetInitialPosition();

            DDa1_DD11 += rDDDN_DDDe(k, 0) * coords;
            DDa1_DD12 += rDDDN_DDDe(k, 1) * coords;
            DDa2_DD21 += rDDDN_DDDe(k, 2) * coords;
            DDa2_DD22 += rDDDN_DDDe(k, 3) * coords;
        } 
        
        // derivative of the curvature
        array_1d<double, 3> cross1;
        array_1d<double, 3> cross2;
        array_1d<double, 3> Da1_D1;
        array_1d<double, 3> Da2_D2;
        array_1d<double, 3> Da1_D2;

        for (unsigned int i = 0; i < 3; i++)
        {
            Da1_D1[i] = rHessian(i, 0);
            Da2_D2[i] = rHessian(i, 1);
            Da1_D2[i] = rHessian(i, 2);
        }

        MathUtils<double>::CrossProduct(cross1, Da1_D1, a2);
        MathUtils<double>::CrossProduct(cross2, a1, Da1_D2);
        array_1d<double, 3> Da3_D1 = (m_dA_vector[IntegrationPointIndex] * (cross1 + cross2) - a3_tilde * inner_prod(a3,(cross1 + cross2))/m_dA_vector[IntegrationPointIndex]) 
            / pow(m_dA_vector[IntegrationPointIndex], 2);
        MathUtils<double>::CrossProduct(cross1, Da1_D2, a2);
        MathUtils<double>::CrossProduct(cross2, a1, Da2_D2);
        array_1d<double, 3> Da3_D2 = (m_dA_vector[IntegrationPointIndex] * (cross1 + cross2) - a3_tilde * inner_prod(a3,(cross1 + cross2))/m_dA_vector[IntegrationPointIndex])
            / pow(m_dA_vector[IntegrationPointIndex], 2);
        rDCurvature_D1[0] = inner_prod(DDa1_DD11, a3) + inner_prod(Da1_D1, Da3_D1);
        rDCurvature_D1[1] = inner_prod(DDa2_DD21, a3) + inner_prod(Da2_D2, Da3_D1);
        rDCurvature_D1[2] = inner_prod(DDa1_DD12, a3) + inner_prod(Da1_D2, Da3_D1);
        rDCurvature_D2[0] = inner_prod(DDa1_DD12, a3) + inner_prod(Da1_D1, Da3_D2);
        rDCurvature_D2[1] = inner_prod(DDa2_DD22, a3) + inner_prod(Da2_D2, Da3_D2);
        rDCurvature_D2[2] = inner_prod(DDa2_DD21, a3) + inner_prod(Da1_D2, Da3_D2);
    }
    
    void Shell3pElement::CalculateDerivativeOfCurvatureActual(
        const IndexType IntegrationPointIndex,
        array_1d<double, 3>& rDCurvature_D1,
        array_1d<double, 3>& rDCurvature_D2,
        const Matrix& rHessian,
        const KinematicVariables& rKinematicVariables) const
    {
        // second derivatives of the base vectors w.r.t. the curvilinear coords.
        const Matrix& rDDDN_DDDe = GetGeometry().ShapeFunctionDerivatives(3, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        array_1d<double, 3> DDa1_DD11 = ZeroVector(3);
        array_1d<double, 3> DDa1_DD12 = ZeroVector(3);
        array_1d<double, 3> DDa2_DD21 = ZeroVector(3);
        array_1d<double, 3> DDa2_DD22 = ZeroVector(3);

        CalculateSecondDerivativesOfBaseVectors(rDDDN_DDDe, DDa1_DD11, DDa1_DD12, DDa2_DD21, DDa2_DD22); 
        
        // derivative of the curvature
        array_1d<double, 3> cross1;
        array_1d<double, 3> cross2;
        array_1d<double, 3> Da1_D1;
        array_1d<double, 3> Da2_D2;
        array_1d<double, 3> Da1_D2;

        for (unsigned int i = 0; i < 3; i++)
        {
            Da1_D1[i] = rHessian(i, 0);
            Da2_D2[i] = rHessian(i, 1);
            Da1_D2[i] = rHessian(i, 2);
        }

        MathUtils<double>::CrossProduct(cross1, Da1_D1, rKinematicVariables.a2);
        MathUtils<double>::CrossProduct(cross2, rKinematicVariables.a1, Da1_D2);
        array_1d<double, 3> Da3_D1 = (rKinematicVariables.dA * (cross1 + cross2) - rKinematicVariables.a3_tilde * inner_prod(rKinematicVariables.a3_tilde,(cross1 + cross2))/rKinematicVariables.dA) 
            / pow(rKinematicVariables.dA, 2);
        MathUtils<double>::CrossProduct(cross1, Da1_D2, rKinematicVariables.a2);
        MathUtils<double>::CrossProduct(cross2, rKinematicVariables.a1, Da2_D2);
        array_1d<double, 3> Da3_D2 = (rKinematicVariables.dA * (cross1 + cross2) - rKinematicVariables.a3_tilde * inner_prod(rKinematicVariables.a3_tilde,(cross1 + cross2))/rKinematicVariables.dA)
            / pow(rKinematicVariables.dA, 2);
        rDCurvature_D1[0] = inner_prod(DDa1_DD11, rKinematicVariables.a3) + inner_prod(Da1_D1, Da3_D1);
        rDCurvature_D1[1] = inner_prod(DDa2_DD21, rKinematicVariables.a3) + inner_prod(Da2_D2, Da3_D1);
        rDCurvature_D1[2] = inner_prod(DDa1_DD12, rKinematicVariables.a3) + inner_prod(Da1_D2, Da3_D1);
        rDCurvature_D2[0] = inner_prod(DDa1_DD12, rKinematicVariables.a3) + inner_prod(Da1_D1, Da3_D2);
        rDCurvature_D2[1] = inner_prod(DDa2_DD22, rKinematicVariables.a3) + inner_prod(Da2_D2, Da3_D2);
        rDCurvature_D2[2] = inner_prod(DDa2_DD21, rKinematicVariables.a3) + inner_prod(Da1_D2, Da3_D2);
    }

    void Shell3pElement::CalculateDerivativeTransformationMatrices(
        const IndexType IntegrationPointIndex,
        std::vector<Matrix>& rDQ_Dalpha_init,
        std::vector<Matrix>& rDTransCartToCov_Dalpha_init,
        const Matrix& rHessian) const
    {
        // base vector at initial configuration
        const SizeType number_of_nodes = GetGeometry().size();
        const Matrix& rDN_De = GetGeometry().ShapeFunctionDerivatives(1, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        array_1d<double, 3> a1 = ZeroVector(3);
        array_1d<double, 3> a2 = ZeroVector(3);
        array_1d<double, 3> a3_tilde = ZeroVector(3);
        array_1d<double, 3> a3 = ZeroVector(3);

        for (SizeType i = 0; i < number_of_nodes; i++)
        {
            a1[0] += GetGeometry().GetPoint( i ).X0() * rDN_De(i, 0);
            a1[1] += GetGeometry().GetPoint( i ).Y0() * rDN_De(i, 0);
            a1[2] += GetGeometry().GetPoint( i ).Z0() * rDN_De(i, 0);

            a2[0] += GetGeometry().GetPoint( i ).X0() * rDN_De(i, 1);
            a2[1] += GetGeometry().GetPoint( i ).Y0() * rDN_De(i, 1);
            a2[2] += GetGeometry().GetPoint( i ).Z0() * rDN_De(i, 1);
        }

        MathUtils<double>::CrossProduct(a3_tilde, a1, a2);
        noalias(a3) = a3_tilde / norm_2(a3_tilde);

        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (m_A_ab_covariant_vector[IntegrationPointIndex][0] * m_A_ab_covariant_vector[IntegrationPointIndex][1]
                - m_A_ab_covariant_vector[IntegrationPointIndex][2] * m_A_ab_covariant_vector[IntegrationPointIndex][2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] =  inv_det_g_ab * m_A_ab_covariant_vector[IntegrationPointIndex][1];
        a_ab_contravariant[1] =  inv_det_g_ab * m_A_ab_covariant_vector[IntegrationPointIndex][0];
        a_ab_contravariant[2] = -inv_det_g_ab * m_A_ab_covariant_vector[IntegrationPointIndex][2];

        //Contravariant base vectors
        array_1d<double, 3> a_contravariant_1 = a1*a_ab_contravariant[0] + a2*a_ab_contravariant[2];
        array_1d<double, 3> a_contravariant_2 = a1*a_ab_contravariant[2] + a2*a_ab_contravariant[1];

        //Local cartesian coordinates
        double l_a1 = norm_2(a1);
        array_1d<double, 3> e1 = a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        // in 1.direction
        array_1d<double, 3> Da1_D1_init;
        array_1d<double, 3> Da1_D2_init;
        array_1d<double, 3> Da2_D2_init;

        for (unsigned int i = 0; i < 3; i++)
        {
            Da1_D1_init[i] = rHessian(i, 0);
            Da2_D2_init[i] = rHessian(i, 1);
            Da1_D2_init[i] = rHessian(i, 2);
        }
        array_1d<double, 3> De1_D1 = Da1_D1_init / l_a1 - inner_prod(e1, Da1_D1_init) * e1 / l_a1;
        array_1d<double, 3> tilde_t2 = a2 - inner_prod(a2, e1) * e1;
        double bar_tilde_t2 = norm_2(tilde_t2);
        array_1d<double, 3> tilde_t2_1 = Da1_D2_init - (inner_prod(Da1_D2_init, e1) + inner_prod(a2,De1_D1)) * e1 
            - inner_prod(a2, e1) * De1_D1;
        array_1d<double, 3> De2_D1 = tilde_t2_1 / bar_tilde_t2 + inner_prod(e2, tilde_t2_1) * e2 / bar_tilde_t2;


        //derivative of covariant base vectors
        double A_1 = 2.0 * inner_prod(Da1_D1_init, a1) * m_A_ab_covariant_vector[IntegrationPointIndex][1] 
            + 2.0 * m_A_ab_covariant_vector[IntegrationPointIndex][0] * inner_prod(Da1_D2_init, a2) 
            - 2.0 * m_A_ab_covariant_vector[IntegrationPointIndex][2] * (inner_prod(Da1_D1_init, a2) + inner_prod(a1, Da1_D2_init)); // check (stands in Carat Code)
        array_1d<double, 3> Da1_con_D1 = inv_det_g_ab 
            * (2.0 * inner_prod(Da1_D2_init, a2) * a1 + m_A_ab_covariant_vector[IntegrationPointIndex][1] * Da1_D1_init 
            - (inner_prod(Da1_D1_init, a2) + inner_prod(a1, Da1_D2_init)) 
            * a2 - m_A_ab_covariant_vector[IntegrationPointIndex][2] * Da1_D2_init) 
            - pow(inv_det_g_ab, 2) * (m_A_ab_covariant_vector[IntegrationPointIndex][1] * a1 - m_A_ab_covariant_vector[IntegrationPointIndex][2] * a2) * A_1;
        array_1d<double, 3> Da2_con_D1 = inv_det_g_ab
            * (-(inner_prod(Da1_D2_init, a1) + inner_prod(a2, Da1_D1_init)) * a1 
            + m_A_ab_covariant_vector[IntegrationPointIndex][2] * Da1_D1_init + 2.0 * inner_prod(Da1_D1_init, a1) * a2 
            - m_A_ab_covariant_vector[IntegrationPointIndex][0] * Da1_D2_init) - pow(inv_det_g_ab, 2) 
            * (-m_A_ab_covariant_vector[IntegrationPointIndex][2] * a1 + m_A_ab_covariant_vector[IntegrationPointIndex][0] * a2) * A_1;
  
        double eG11 = inner_prod(e1, a_contravariant_1);
        double eG12 = inner_prod(e1, a_contravariant_2);
        double eG21 = inner_prod(e2, a_contravariant_1);
        double eG22 = inner_prod(e2, a_contravariant_2);
        double eG11_d1 = inner_prod(De1_D1, a_contravariant_1) + inner_prod(e1, Da1_con_D1);
        double eG12_d1 = inner_prod(De1_D1, a_contravariant_2) + inner_prod(e1, Da2_con_D1);
        double eG21_d1 = inner_prod(De2_D1, a_contravariant_1) + inner_prod(e2, Da1_con_D1);
        double eG22_d1 = inner_prod(De2_D1, a_contravariant_2) + inner_prod(e2, Da2_con_D1);
        
        // derivative of the transformation matrix T_con_to_car (contravariant to local Cartesian basis) of the initial configuration w.r.t. theta1
        rDQ_Dalpha_init[0](0,0) = eG11_d1 * eG11 + eG11 * eG11_d1;
        rDQ_Dalpha_init[0](0,1) = eG12_d1 * eG12 + eG12 * eG12_d1;
        rDQ_Dalpha_init[0](0,2) = 2.0 * (eG11_d1 * eG12 + eG11 * eG12_d1);
        rDQ_Dalpha_init[0](1,0) = eG21_d1 * eG21 + eG21 * eG21_d1;  // should be always zero (ML)
        rDQ_Dalpha_init[0](1,1) = eG22_d1 * eG22 + eG22 * eG22_d1;
        rDQ_Dalpha_init[0](1,2) = 2.0 * (eG21_d1 * eG22 + eG21 * eG22_d1);  // should be always zero (ML)
        rDQ_Dalpha_init[0](2,0) = 2.0 * (eG11_d1 * eG21 + eG11 * eG21_d1);  // should be always zero (ML)
        rDQ_Dalpha_init[0](2,1) = 2.0 * (eG12_d1 * eG22 + eG12 * eG22_d1);
        rDQ_Dalpha_init[0](2,2) = 2.0 * (eG11_d1 * eG22 + eG12_d1 * eG21 + eG11 * eG22_d1 + eG12 * eG21_d1);
        // derivative of the transformation matrix T_car_to_cov (local Cartesian to covariant basis) of the initial configuration 
        // w.r.t. theta1
        rDTransCartToCov_Dalpha_init[0](0,0) = eG11_d1 * eG11 + eG11 * eG11_d1;
        rDTransCartToCov_Dalpha_init[0](0,1) = eG21_d1 * eG21 + eG21 * eG21_d1;
        rDTransCartToCov_Dalpha_init[0](0,2) = 2.0 * (eG11_d1 * eG21 + eG11 * eG21_d1);
        rDTransCartToCov_Dalpha_init[0](1,0) = eG12_d1 * eG12 + eG12 * eG12_d1;
        rDTransCartToCov_Dalpha_init[0](1,1) = eG22_d1 * eG22 + eG22 * eG22_d1;
        rDTransCartToCov_Dalpha_init[0](1,2) = 2.0 * (eG12_d1 * eG22 + eG12 * eG22_d1);
        rDTransCartToCov_Dalpha_init[0](2,0) = eG11_d1 * eG12 + eG11 * eG12_d1;
        rDTransCartToCov_Dalpha_init[0](2,1) = eG21_d1 * eG22 + eG21 * eG22_d1;
        rDTransCartToCov_Dalpha_init[0](2,2) = eG11_d1 * eG22 + eG21_d1 * eG12 + eG11 * eG22_d1 + eG21 * eG12_d1;

        // in 2.direction
        array_1d<double, 3> De1_D2 = Da1_D2_init / l_a1 + inner_prod(e1, Da1_D2_init) * e1 / l_a1;
        array_1d<double, 3> tilde_t2_2 = Da2_D2_init - (inner_prod(Da2_D2_init, e1) + inner_prod(a2, De1_D2)) * e1 
            - inner_prod(a2, e1) * De1_D2;
        array_1d<double, 3> De2_D2 = tilde_t2_2 / bar_tilde_t2 + inner_prod(e2, tilde_t2_2) * e2 / bar_tilde_t2;
  
        //derivative of covariant base vectors
        double A_2 = 2.0 * inner_prod(Da1_D2_init, a1) * m_A_ab_covariant_vector[IntegrationPointIndex][1] 
            + 2.0 * m_A_ab_covariant_vector[IntegrationPointIndex][0] * inner_prod(Da2_D2_init, a2) 
            - 2.0 * m_A_ab_covariant_vector[IntegrationPointIndex][2] * (inner_prod(Da1_D2_init, a2) + inner_prod(a1, Da2_D2_init));
        array_1d<double, 3> Da1_con_D2 = inv_det_g_ab * (2.0 * inner_prod(Da2_D2_init, a2) * a1 
            + m_A_ab_covariant_vector[IntegrationPointIndex][1] * Da1_D2_init - (inner_prod(Da1_D2_init, a2) 
            + inner_prod(a1, Da2_D2_init)) * a2 - m_A_ab_covariant_vector[IntegrationPointIndex][2] * Da2_D2_init) 
            - pow(inv_det_g_ab, 2) * (m_A_ab_covariant_vector[IntegrationPointIndex][1] * a1 - m_A_ab_covariant_vector[IntegrationPointIndex][2] * a2) * A_2;
        array_1d<double, 3> Da2_con_D2 = inv_det_g_ab * (-(inner_prod(Da2_D2_init, a1) 
            + inner_prod(a2, Da1_D2_init)) * a1 + m_A_ab_covariant_vector[IntegrationPointIndex][2] * Da1_D2_init 
            + 2.0 * inner_prod(Da1_D2_init, a1) * a2 - m_A_ab_covariant_vector[IntegrationPointIndex][0] * Da2_D2_init) 
            - pow(inv_det_g_ab, 2) * (-m_A_ab_covariant_vector[IntegrationPointIndex][2] * a1 + m_A_ab_covariant_vector[IntegrationPointIndex][0] * a2) * A_2;

        double eG11_d2 = inner_prod(De1_D2, a_contravariant_1) + inner_prod(e1, Da1_con_D2);
        double eG12_d2 = inner_prod(De1_D2, a_contravariant_2) + inner_prod(e1, Da2_con_D2);
        double eG21_d2 = inner_prod(De2_D2, a_contravariant_1) + inner_prod(e2, Da1_con_D2);
        double eG22_d2 = inner_prod(De2_D2, a_contravariant_2) + inner_prod(e2, Da2_con_D2);

        // derivative of the transformation matrix T_con_to_car (contravariant to local Cartesian basis) of the initial configuration w.r.t. theta2
        rDQ_Dalpha_init[1](0,0) = eG11_d2 * eG11 + eG11 * eG11_d2;
        rDQ_Dalpha_init[1](0,1) = eG12_d2 * eG12 + eG12 * eG12_d2;
        rDQ_Dalpha_init[1](0,2) = 2.0 * (eG11_d2 * eG12 + eG11 * eG12_d2);
        rDQ_Dalpha_init[1](1,0) = eG21_d2 * eG21 + eG21 * eG21_d2;
        rDQ_Dalpha_init[1](1,1) = eG22_d2 * eG22 + eG22 * eG22_d2;
        rDQ_Dalpha_init[1](1,2) = 2.0 * (eG21_d2 * eG22 + eG21 * eG22_d2);
        rDQ_Dalpha_init[1](2,0) = 2.0 * (eG11_d2 * eG21 + eG11 * eG21_d2);
        rDQ_Dalpha_init[1](2,1) = 2.0 * (eG12_d2 * eG22 + eG12 * eG22_d2);
        rDQ_Dalpha_init[1](2,2) = 2.0 * (eG11_d2 * eG22 + eG12_d2 * eG21 + eG11 * eG22_d2 + eG12 * eG21_d2);
        
        // derivative of the transformation matrix T_car_to_cov (local Cartesian to covariant basis) of the initial configuration 
        // w.r.t. theta2
        rDTransCartToCov_Dalpha_init[1](0,0) = eG11_d2 * eG11 + eG11 * eG11_d2;
        rDTransCartToCov_Dalpha_init[1](0,1) = eG21_d2*eG21 + eG21*eG21_d2;
        rDTransCartToCov_Dalpha_init[1](0,2) = 2.0*(eG11_d2*eG21 + eG11*eG21_d2);
        rDTransCartToCov_Dalpha_init[1](1,0) = eG12_d2*eG12 + eG12*eG12_d2;
        rDTransCartToCov_Dalpha_init[1](1,1) = eG22_d2*eG22 + eG22*eG22_d2;
        rDTransCartToCov_Dalpha_init[1](1,2) = 2.0*(eG12_d2*eG22 + eG12*eG22_d2);
        rDTransCartToCov_Dalpha_init[1](2,0) = eG11_d2*eG12 + eG11*eG12_d2;
        rDTransCartToCov_Dalpha_init[1](2,1) = eG21_d2*eG22 + eG21*eG22_d2;
        rDTransCartToCov_Dalpha_init[1](2,2) = eG11_d2*eG22 + eG21_d2*eG12 + eG11*eG22_d2 + eG21*eG12_d2;
    }

    ///@}
    ///@name Dynamic Functions
    ///@{

    void Shell3pElement::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void Shell3pElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const IndexType index = i * 3;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void Shell3pElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const IndexType index = i * 3;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void Shell3pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    };

    void Shell3pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    };

    ///@}
    ///@name Check
    ///@{

    int Shell3pElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        // Verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }
        else
        {
            // Verify that the constitutive law has the correct dimension
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
                << "THICKNESS not provided for element " << this->Id() << std::endl;

            // Check strain size
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
                << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
                << this->Id() << std::endl;
        }

        return 0;
    }

    void Shell3pElement::CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe) const
    {
        const SizeType number_of_points = GetGeometry().size();
        const SizeType working_space_dimension = 3;
        Hessian.resize(working_space_dimension, working_space_dimension);
        Hessian = ZeroMatrix(working_space_dimension, working_space_dimension);

        for (IndexType k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

            Hessian(0, 0) += rDDN_DDe(k, 0)*coords[0];
            Hessian(0, 1) += rDDN_DDe(k, 2)*coords[0];
            Hessian(0, 2) += rDDN_DDe(k, 1)*coords[0];

            Hessian(1, 0) += rDDN_DDe(k, 0)*coords[1];
            Hessian(1, 1) += rDDN_DDe(k, 2)*coords[1];
            Hessian(1, 2) += rDDN_DDe(k, 1)*coords[1];

            Hessian(2, 0) += rDDN_DDe(k, 0)*coords[2];
            Hessian(2, 1) += rDDN_DDe(k, 2)*coords[2];
            Hessian(2, 2) += rDDN_DDe(k, 1)*coords[2];
        }
    }

    void Shell3pElement::CalculateSecondDerivativesOfBaseVectors(
        const Matrix& rDDDN_DDDe,
        array_1d<double, 3>& rDDa1_DD11,
        array_1d<double, 3>& rDDa1_DD12,
        array_1d<double, 3>& rDDa2_DD21,
        array_1d<double, 3>& rDDa2_DD22) const
    {
        const SizeType number_of_points = GetGeometry().size();
    
        for (IndexType k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

            rDDa1_DD11 += rDDDN_DDDe(k, 0) * coords;
            rDDa1_DD12 += rDDDN_DDDe(k, 1) * coords;
            rDDa2_DD21 += rDDDN_DDDe(k, 2) * coords;
            rDDa2_DD22 += rDDDN_DDDe(k, 3) * coords;
        }
    }

    ///@}

} // Namespace Kratos


