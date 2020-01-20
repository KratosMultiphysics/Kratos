//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Alexander Müller
//                   Tobias Teschemacher
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/shell_5p_element.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void Shell5pElement::Initialize()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (m_B_ab_covariant_vector.size() != r_number_of_integration_points)
            m_B_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_dA_vector.size() != r_number_of_integration_points)
            m_dA_vector.resize(r_number_of_integration_points);
        if (m_cart_deriv.size() != r_number_of_integration_points)
            m_cart_deriv.resize(r_number_of_integration_points);

        KinematicVariables kinematic_variables(WorkingSpaceDimension());
        VariationVariables variation_variables;

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
        {
            CalculateCartesianDerivatives(point_number, kinematic_variables, m_cart_deriv[point_number]);

            CalculateKinematics(
                point_number,
                kinematic_variables,variation_variables);

            m_B_ab_covariant_vector[point_number] = kinematic_variables.b_ab_covariant;
            m_S_a_covariant_vector [point_number] = kinematic_variables.s_a_covariant;
        }
        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void Shell5pElement::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_number_of_integration_points)
            mConstitutiveLawVector.resize(r_number_of_integration_points);


        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Assembly
    ///@{

    void Shell5pElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 5;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(
                WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane(3);
            ConstitutiveVariables constitutive_variables_curvature(3);
            ConstitutiveVariables constitutive_variables_shear(2);
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_variables_curvature,
                constitutive_variables_shear,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            Matrix BMembrane = ZeroMatrix(8, mat_size);
            Matrix BCurvature = ZeroMatrix(8, mat_size);
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
    ///@name Kinematics
    ///@{

    void Shell5pElement::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKin,
        VariationVariables& rVar
    )
    {
        Matrix J;
        GetGeometry().Jacobian(J, IntegrationPointIndex);
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);

        rKin.a1 = column(J, 0);
        rKin.a2 = column(J, 1);

        //GetCovariantMetric
        rKin.a_ab_covariant[0] = norm_2_square(rKin.a1); //TODO calculate metric from displacements and not from geometry. This is more accurate for small strains due cancelation with reference metric since 1.0-1.0 ~ loses digits 
        rKin.a_ab_covariant[1] = norm_2_square(rKin.a2);
        rKin.a_ab_covariant[2] = inner_prod(rKin.a1, rKin.a2);

        rKin.t = ZeroVector(3);
        rKin.dtd1 = ZeroVector(3);
        rKin.dtd2 = ZeroVector(3);
        for (size_t i = 0; i < N.size; i++) 
        {
            rKin.t += N[i]* GetGeometry()[i].FastGetSolutionStepValue(DIRECTOR); //TODO does this really do what i want? t=sum_I N_I *t_I ? and is this the best way to obtain this?
            rKin.dtd1 += m_cart_deriv[IntegrationPointIndex](0,i)* GetGeometry()[i].FastGetSolutionStepValue(DIRECTOR);
            rKin.dtd2 += m_cart_deriv[IntegrationPointIndex](1,i)* GetGeometry()[i].FastGetSolutionStepValue(DIRECTOR);
        }

        double invL_t = 1.0 / norm_2(rKin.t);
        rKin.t *= invL_t;

        array_1d<double, 3>& t = rKin.t; //define alias for cleaner code
        array_1d<double, 3>& dtd1 = rKin.dtd1; //define alias for cleaner code
        array_1d<double, 3>& dtd2 = rKin.dtd2; //define alias for cleaner code
        array_1d<double, 3>& a1 = rKin.a1; //define alias for cleaner code
        array_1d<double, 3>& a2 = rKin.a2; //define alias for cleaner code

        rVar.P = outer_prod(t, t) - IdentityMatrix(3);
        rVar.P *= invL_t;

        invL_t *= invL_t;
        const Matrix tdyadt = outer_prod(t, t);
        rVar.Q1 = invL_t * ((inner_prod(t, dtd1))*(3.0* tdyadt -IdentityMatrix(3)) - outer_prod(dtd1, t)- outer_prod(t, dtd1));
        rVar.Q2 = invL_t * ((inner_prod(t, dtd2))*(3.0* tdyadt -IdentityMatrix(3)) - outer_prod(dtd2, t)- outer_prod(t, dtd2));
        rVar.S1 = invL_t * ((inner_prod(t,   a1))*(3.0* tdyadt -IdentityMatrix(3)) - outer_prod(  a1, t)- outer_prod(t,   a1));
        rVar.S2 = invL_t * ((inner_prod(t,   a2))*(3.0* tdyadt -IdentityMatrix(3)) - outer_prod(  a2, t)- outer_prod(t,   a2));

        invL_t *= invL_t;
        rVar.Chi11 = invL_t * (3.0 * inner_prod(t, dtd1) * (outer_prod(a1, t) + 0.5 * inner_prod(a1, t) * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a1, dtd1) * tdyadt + inner_prod(a1, t) * outer_prod(dtd1, t)) - outer_prod(a1, dtd1) - inner_prod(a1, dtd1) * 0.5 * IdentityMatrix(3));
        rVar.Chi21 = invL_t * (3.0 * inner_prod(t, dtd1) * (outer_prod(a2, t) + 0.5 * inner_prod(a2, t) * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a2, dtd1) * tdyadt + inner_prod(a2, t) * outer_prod(dtd1, t)) - outer_prod(a2, dtd1) - inner_prod(a2, dtd1) * 0.5 * IdentityMatrix(3));
        rVar.Chi12 = invL_t * (3.0 * inner_prod(t, dtd2) * (outer_prod(a1, t) + 0.5 * inner_prod(a1, t) * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a1, dtd2) * tdyadt + inner_prod(a1, t) * outer_prod(dtd2, t)) - outer_prod(a1, dtd2) - inner_prod(a1, dtd2) * 0.5 * IdentityMatrix(3));
        rVar.Chi22 = invL_t * (3.0 * inner_prod(t, dtd2) * (outer_prod(a2, t) + 0.5 * inner_prod(a2, t) * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a2, dtd2) * tdyadt + inner_prod(a2, t) * outer_prod(dtd2, t)) - outer_prod(a2, dtd2) - inner_prod(a2, dtd2) * 0.5 * IdentityMatrix(3));

        rKin.dtd1 = prod(rVar.P, rKin.dtd1);
        rKin.dtd2 = prod(rVar.P, rKin.dtd2);

        rKin.b_ab_covariant[0] = inner_prod(rKin.a1, rKin.dtd1);
        rKin.b_ab_covariant[1] = inner_prod(rKin.a2, rKin.dtd2);
        rKin.b_ab_covariant[2] = inner_prod(rKin.a1, rKin.dtd2) + inner_prod(rKin.a2, rKin.dtd1);

        rKin.s_a_covariant[0] = inner_prod(rKin.t, rKin.a1);
        rKin.s_a_covariant[1] = inner_prod(rKin.t, rKin.a2);
    }

    /* Transforms derivatives to obtain cartesian quantities  */
    void Shell5pElement::CalculateCartesianDerivatives(
        IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        Matrix& mCartD
    )
    {
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        array_1d<double, 3> a3;
        MathUtils<double>::CrossProduct(a3, rKinematicVariables.a1, rKinematicVariables.a2);
        m_dA_vector[IntegrationPointIndex] = norm_2(a3);

        array_1d<double, 3> axi1 = rKinematicVariables.a1/ norm_2(rKinematicVariables.a1);
        array_1d<double, 3> axi2 = rKinematicVariables.a2 / norm_2(rKinematicVariables.a2);

        MathUtils<double>::CrossProduct(a3, axi1, axi2);

        array_1d<double, 3> axi1bar = 0.5 * (axi1 + axi2);
        axi1bar = axi1bar / norm_2(axi1bar);
        array_1d<double, 3> axi2bar;
        MathUtils<double>::CrossProduct(axi2bar, a3, axi1bar);

        array_1d<double, 3> a1Cart = .70710678118654752440 * (axi1bar - axi2bar);
        array_1d<double, 3> a2Cart = .70710678118654752440 * (axi1bar + axi2bar);

        Matrix J = ZeroMatrix(2, 2);

        J(0, 0) = inner_prod(rKinematicVariables.a1, a1Cart);
        J(1, 0) = inner_prod(rKinematicVariables.a2, a1Cart);
        J(0, 1) = inner_prod(rKinematicVariables.a1, a2Cart);
        J(1, 1) = inner_prod(rKinematicVariables.a2, a2Cart);

        double detJ;
        Matrix invJ;
        MathUtils<double>::InvertMatrix2(J,invJ, detJ);
  
        mCartD = prod(invJ, r_DN_De);
    }

    void Shell5pElement::CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveVariables& rThisConstitutiveVariablesShear,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        array_1d<double, 3> strain_vector =  (rActualKinematic.a_ab_covariant);
        strain_vector(0) -= 1.0; //TODO this is unnecessary if TODO from callculatekinematics is done
        strain_vector(1) -= 1.0;

        strain_vector(0) *= 0.5;
        strain_vector(1) *= 0.5;

        noalias(rThisConstitutiveVariablesMembrane.StrainVector) =  strain_vector;

        array_1d<double, 3> curvature_vector = rActualKinematic.b_ab_covariant - m_B_ab_covariant_vector[IntegrationPointIndex];
        noalias(rThisConstitutiveVariablesCurvature.StrainVector) =  curvature_vector;

        // Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

        mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);

        double thickness = this->GetProperties().GetValue(THICKNESS);
        noalias(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix) = rThisConstitutiveVariablesMembrane.ConstitutiveMatrix*(pow(thickness, 2) / 12);   //TODO this does not work for general material laws especially including shear

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
        noalias(rThisConstitutiveVariablesCurvature.StressVector) = prod(
            trans(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix), rThisConstitutiveVariablesCurvature.StrainVector);

        //TODO add shear constiutive relatons
    }

    void Shell5pElement::CalculateBOperator(
        IndexType IntegrationPointIndex,
        Matrix& rB,
        const KinematicVariables& rActualKinematic,
        const VariationVariables& rVariations )
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType num_dof = number_of_control_points * 5;

        const Vector& r_N = this->GetValue(SHAPE_FUNCTION_VALUES);

        noalias(rB) = ZeroMatrix(8, num_dof);

        Matrix WI1(3, 3);
        Matrix WI2(3, 3);

        for (IndexType r = 0; r < number_of_control_points; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = 5 * r;

            WI1 = rVariations.Q1 * r_N[r] + rVariations.P * m_cart_deriv[IntegrationPointIndex](0, r);
            WI2 = rVariations.Q2 * r_N[r] + rVariations.P * m_cart_deriv[IntegrationPointIndex](1, r);

            //membrane
            subrange(rB, 0, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.a1); //TODO WHY DOES THIS NOT WORK is  the equal sign not overloaded for array_1d?
            subrange(rB, 1, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.a2);
            subrange(rB, 2, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.a2) + m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.a1);

            //bending
            subrange(rB, 3, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.dtd1); //TODO WHY DOES THIS NOT WORK
            subrange(rB, 4, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.dtd2);
            subrange(rB, 5, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.dtd2) + m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.dtd1);

            subrange(rB, 3, 1, kr + 3, 2) = prod(trans(rActualKinematic.a1), WI1); //TODO WHY DOES THIS NOT WORK
            subrange(rB, 4, 1, kr + 3, 2) = prod(trans(rActualKinematic.a2), WI2); //TODO WHY DOES THIS NOT WORK
            subrange(rB, 5, 1, kr + 3, 2) = prod(trans(rActualKinematic.a2), WI1) + prod(trans(rActualKinematic.a1), WI2); //TODO WHY DOES THIS NOT WORK


            //shear
            subrange(rB, 6, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.t); //TODO WHY DOES THIS NOT WORK
            subrange(rB, 7, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.t);

            subrange(rB, 6, 1, kr + 3, 2) = prod(trans(rActualKinematic.a1), rVariations.P) * r_N[r];
            subrange(rB, 7, 1, kr + 3, 2) = prod(trans(rActualKinematic.a2), rVariations.P) * r_N[r];
        };


        }
    }

    

    void Shell5pElement::CalculateGeometricStiffness(
        IndexType IntegrationPointIndex,
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const KinematicVariables& rActualKinematic)
    {
        const auto& r_geometry = GetGeometry();

        const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;


        Matrix S_da3 = ZeroMatrix(3, mat_size);
        Vector S_a3_da3 = ZeroVector(mat_size);
        Vector S_a3_da3_l_a3_3 = ZeroVector(mat_size);
        Matrix S_dn = ZeroMatrix(3, mat_size);

        // first variation of strain and curvature w.r.t. dof
        for (int r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 3;
            int dirr = r % 3;

            array_1d<double, 3> S_dg_1 = ZeroVector(3);
            array_1d<double, 3> S_dg_2 = ZeroVector(3);
            S_dg_1(dirr) = r_DN_De(0, kr);
            S_dg_2(dirr) = r_DN_De(1, kr);

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
        for (int r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 3;
            int dirr = r % 3;

            for (int s = 0; s <= r; s++)
            {
                // local node number ks and dof direction dirs
                int ks = s / 3;
                int dirs = s % 3;

                // strain
                array_1d<double, 3> ddE_cu = ZeroVector(3);
                if (dirr == dirs)
                {
                    ddE_cu[0] = r_DN_De(0, kr) * r_DN_De(0, ks);
                    ddE_cu[1] = r_DN_De(1, kr) * r_DN_De(1, ks);
                    ddE_cu[2] = 0.5 * (r_DN_De(0, kr) * r_DN_De(1, ks) + r_DN_De(1, kr) * r_DN_De(0, ks));

                    rSecondVariationsStrain.B11(r, s) = m_T_vector[IntegrationPointIndex](0, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](0, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](0, 2) * ddE_cu[2];
                    rSecondVariationsStrain.B22(r, s) = m_T_vector[IntegrationPointIndex](1, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](1, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](1, 2) * ddE_cu[2];
                    rSecondVariationsStrain.B12(r, s) = m_T_vector[IntegrationPointIndex](2, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](2, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](2, 2) * ddE_cu[2];
                }

                // curvature
                array_1d<double, 3> dda3 = ZeroVector(3);
                int dirt = 4 - dirr - dirs;
                int ddir = dirr - dirs;
                     if (ddir == -1) dda3(dirt - 1) =  r_DN_De(0, kr) * r_DN_De(1, ks) - r_DN_De(0, ks) * r_DN_De(1, kr);
                else if (ddir ==  2) dda3(dirt - 1) =  r_DN_De(0, kr) * r_DN_De(1, ks) - r_DN_De(0, ks) * r_DN_De(1, kr);
                else if (ddir ==  1) dda3(dirt - 1) = -r_DN_De(0, kr) * r_DN_De(1, ks) + r_DN_De(0, ks) * r_DN_De(1, kr);
                else if (ddir == -2) dda3(dirt - 1) = -r_DN_De(0, kr) * r_DN_De(1, ks) + r_DN_De(0, ks) * r_DN_De(1, kr);

                double c = -(dda3[0] * rActualKinematic.a3_tilde[0] + dda3[1] * rActualKinematic.a3_tilde[1] + dda3[2] * rActualKinematic.a3_tilde[2]
                    + S_da3(0, r) * S_da3(0, s) + S_da3(1, r) * S_da3(1, s) + S_da3(2, r) * S_da3(2, s)
                    ) * inv_l_a3_3;

                double d = 3.0 * S_a3_da3[r] * S_a3_da3[s] * inv_l_a3_5;

                array_1d<double, 3> ddn = ZeroVector(3);
                ddn[0] = dda3[0] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(0, r) - S_a3_da3_l_a3_3[r] * S_da3(0, s) + (c + d) * rActualKinematic.a3_tilde[0];
                ddn[1] = dda3[1] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(1, r) - S_a3_da3_l_a3_3[r] * S_da3(1, s) + (c + d) * rActualKinematic.a3_tilde[1];
                ddn[2] = dda3[2] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(2, r) - S_a3_da3_l_a3_3[r] * S_da3(2, s) + (c + d) * rActualKinematic.a3_tilde[2];

                array_1d<double, 3> ddK_cu = ZeroVector(3);
                ddK_cu[0] = r_DDN_DDe(0, kr) * S_dn(dirr, s) + r_DDN_DDe(0, ks) * S_dn(dirs, r)
                    + H(0, 0) * ddn[0] + H(1, 0) * ddn[1] + H(2, 0) * ddn[2];
                ddK_cu[1] = r_DDN_DDe(1, kr) * S_dn(dirr, s) + r_DDN_DDe(1, ks) * S_dn(dirs, r)
                    + H(0, 1) * ddn[0] + H(1, 1) * ddn[1] + H(2, 1) * ddn[2];
                ddK_cu[2] = r_DDN_DDe(2, kr) * S_dn(dirr, s) + r_DDN_DDe(2, ks) * S_dn(dirs, r)
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

    inline void Shell5pElement::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rB,
        const Matrix& rD,
        const double IntegrationWeight
    )
    {
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(rB), Matrix(prod(rD, rB)));
    }

    inline void Shell5pElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& rSecondVariationsStrain,
        const Vector& rSD,
        const double IntegrationWeight)
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        for (int n = 0; n < mat_size; n++)
        {
            for (int m = 0; m <= n; m++)
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
    ///@name Dynamic Functions
    ///@{

    void Shell5pElement::GetValuesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 5;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const int index = i * 3; //TODO Geometry are 6 parameters but degrees of freedeom are 5 how to deal with this?

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void Shell5pElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * 3;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void Shell5pElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * 3;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void Shell5pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const SizeType index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    };

    void Shell5pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
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

    int Shell5pElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)

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

    ///@}

} // Namespace Kratos


