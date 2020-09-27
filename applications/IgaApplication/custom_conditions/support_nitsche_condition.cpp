//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//                   Tobias Tescheamacher
//                   Michael Breitenberger
//                   Riccardo Rossi
//

// System includes

// External includes
#include "custom_conditions/support_nitsche_condition.h"

// Project includes

namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void SupportNitscheCondition::Initialize()
    {
        KRATOS_TRY

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void SupportNitscheCondition::InitializeMaterial()
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();
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

    void SupportNitscheCondition::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ConfigurationType& rConfiguration
    )
    {
        const auto& r_geometry = GetGeometry();

        // pass/call this ShapeFunctionsLocalGradients[pnt]
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType number_of_nodes = r_geometry.size();
        Vector g1 = ZeroVector(dimension);
        Vector g2 = ZeroVector(dimension);

        Vector current_displacement = ZeroVector(dimension*number_of_nodes);

        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);

        for (SizeType i=0;i<number_of_nodes;++i){
            g1[0] += (r_geometry.GetPoint( i ).X0()+current_displacement[i*dimension]) * rShapeFunctionGradientValues(i, 0);
            g1[1] += (r_geometry.GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 0);
            g1[2] += (r_geometry.GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 0);

            g2[0] += (r_geometry.GetPoint( i ).X0()+current_displacement[i*dimension]) * rShapeFunctionGradientValues(i, 1);
            g2[1] += (r_geometry.GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 1);
            g2[2] += (r_geometry.GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 1);
        }

        rKinematicVariables.a1 = g1;
        rKinematicVariables.a2 = g2;

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

        //Compute the tangent and  the normal to the boundary vector
        array_1d<double, 2> local_tangents = GetProperties()[LOCAL_TANGENTS_MASTER];

        rKinematicVariables.t = local_tangents[0]*g1 + local_tangents[1]*g2;
        rKinematicVariables.t = rKinematicVariables.t/norm_2(rKinematicVariables.t);

        MathUtils<double>::CrossProduct(rKinematicVariables.n, rKinematicVariables.t, rKinematicVariables.a3);

        // transform the normal into the contavariant basis
        rKinematicVariables.n_contravariant[0] = rKinematicVariables.a1[0]*rKinematicVariables.n[0] + rKinematicVariables.a1[1]*rKinematicVariables.n[1] + rKinematicVariables.a1[2]*rKinematicVariables.n[2];
        rKinematicVariables.n_contravariant[1] = rKinematicVariables.a2[0]*rKinematicVariables.n[0] + rKinematicVariables.a2[1]*rKinematicVariables.n[1] + rKinematicVariables.a2[2]*rKinematicVariables.n[2];
    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void SupportNitscheCondition::CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& rT, Matrix& rT_hat, array_1d<array_1d<double, 3>,2>& rReferenceContraVariantBase
    )
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
        array_1d<double, 3> a_contravariant_1 = rKinematicVariables.a1*a_ab_contravariant[0] + rKinematicVariables.a2*a_ab_contravariant[2];
        array_1d<double, 3> a_contravariant_2 = rKinematicVariables.a1*a_ab_contravariant[2] + rKinematicVariables.a2*a_ab_contravariant[1];

        rReferenceContraVariantBase[0] = a_contravariant_1;
        rReferenceContraVariantBase[1] = a_contravariant_2;

        //Local cartesian coordinates
        double l_a1 = norm_2(rKinematicVariables.a1);
        array_1d<double, 3> e1 = rKinematicVariables.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        //Transformation matrix T from contravariant to local cartesian basis
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


        //Transformation matrix T from local cartesian basis to the covariant one
        if (rT_hat.size1() != 3 && rT_hat.size2() != 3)
            rT_hat.resize(3, 3);
        noalias(rT_hat) = ZeroMatrix(3, 3);

        rT_hat(0, 0) = pow(G(0, 0), 2);
        rT_hat(0, 1) = pow(G(1, 0), 2);
        rT_hat(0, 2) = 2*G(0, 0) * G(1, 0);

        rT_hat(1, 0) = pow(G(0, 1), 2);
        rT_hat(1, 1) = pow(G(1, 1), 2);
        rT_hat(1, 2) = 2*G(0, 1) * G(1, 1);

        rT_hat(2, 0) = G(0, 0) * G(0, 1);
        rT_hat(2, 1) = G(1, 0) * G(1, 1);
        rT_hat(2, 2) = (G(0, 0) * G(1, 1) + G(1, 0) * G(0, 1));


    }

    void SupportNitscheCondition::CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector[IntegrationPointIndex]);
        noalias(rThisConstitutiveVariablesMembrane.StrainVector) = prod(m_T_vector[IntegrationPointIndex], strain_vector);

        // Constitive Matrices DMembrane
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

        mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        rThisConstitutiveVariablesMembrane.ConstitutiveMatrix *= GetProperties()[THICKNESS];

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
    }

    //Prestress Transformation Matrix
    void SupportNitscheCondition::CalculateTransformationmatrixPrestress(
        const KinematicVariables& rActualKinematic,
        PrestresstransVariables& rPrestresstransVariables
    )
    {
        //define base vector in reference plane
        //ATTENTION: in some cases the vector must be modified (e.g. catenoid t3_unten=[0, 0, 1])

        array_1d<double, 3> t3_unten;
        t3_unten[0] = 0;
        t3_unten[1] = 1;
        t3_unten[2] = 0;

        array_1d<double, 3> t1_z;
        MathUtils<double>::CrossProduct(t1_z, t3_unten, rActualKinematic.a3);

        array_1d<double, 3> t2;
        MathUtils<double>::CrossProduct(t2, rActualKinematic.a3, t1_z);

        array_1d<double, 3> t1_n = t1_z/norm_2(t1_z);
        array_1d<double, 3> t2_n = t2/norm_2(t2);
        //array_1d<double, 3> t3_n = rActualKinematic.a3/norm_2(rActualKinematic.a3);

        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rActualKinematic.a_ab_covariant[0] * rActualKinematic.a_ab_covariant[1]
                - rActualKinematic.a_ab_covariant[2] * rActualKinematic.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] =  inv_det_g_ab * rActualKinematic.a_ab_covariant[1];
        a_ab_contravariant[1] =  inv_det_g_ab * rActualKinematic.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rActualKinematic.a_ab_covariant[2];

        //array_1d<double, 3> a_con_1 = rActualKinematic.a1*a_ab_contravariant[0] + rActualKinematic.a2*a_ab_contravariant[2];
        array_1d<double, 3> a_con_2 = rActualKinematic.a1*a_ab_contravariant[2] + rActualKinematic.a2*a_ab_contravariant[1]; 

        //local cartesian coordinates oriented along the 1st base vector in the ref. config.
        double la1 = norm_2(rActualKinematic.a1);
        array_1d<double, 3> e1 = rActualKinematic.a1 / la1;
        double la_con2 = norm_2(a_con_2);
        array_1d<double, 3> e2 = a_con_2 / la_con2;

        //Transformation matrix from the projected basis T to the local cartesian basis
        double eG11 = inner_prod(e1,t1_n);
        double eG12 = inner_prod(e1,t2_n);
        double eG21 = inner_prod(e2,t1_n);
        double eG22 = inner_prod(e2,t2_n);
    
        rPrestresstransVariables.Tpre = ZeroMatrix(3, 3);
        rPrestresstransVariables.Tpre(0,0) = eG11*eG11;
        rPrestresstransVariables.Tpre(0,1) = eG12*eG12;
        rPrestresstransVariables.Tpre(0,2) = 2.0*eG11*eG12;

        rPrestresstransVariables.Tpre(1,0) = eG21*eG21;
        rPrestresstransVariables.Tpre(1,1) = eG22*eG22;
        rPrestresstransVariables.Tpre(1,2) = 2.0*eG21*eG22;

        rPrestresstransVariables.Tpre(2,0) = eG11*eG21;
        rPrestresstransVariables.Tpre(2,1) = eG12*eG22;
        rPrestresstransVariables.Tpre(2,2) = eG11*eG22+eG12*eG21;       
    }  

    void SupportNitscheCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        const double stabilization_parameter = GetProperties()[NITSCHE_STABILIZATION_PARAMETER];

        const auto& r_geometry = GetGeometry();

        // Size definitions
        const SizeType number_of_nodes = r_geometry.size();

        const SizeType mat_size = 3 * number_of_nodes;

        // Memory allocation
        if (CalculateStiffnessMatrixFlag) {
            if (rLeftHandSideMatrix.size1() != mat_size) {
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
        }
        if (CalculateResidualVectorFlag) {
            if (rRightHandSideVector.size() != mat_size) {
                rRightHandSideVector.resize(mat_size, false);
            }
            rRightHandSideVector = ZeroVector(mat_size);
        }

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        const IntegrationMethod integration_method = r_geometry.GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geometry.ShapeFunctionsLocalGradients(integration_method);

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();
        
        // Prepare memory
        if (m_A_ab_covariant_vector.size() != r_number_of_integration_points)
            m_A_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_dA_vector.size() != r_number_of_integration_points)
            m_dA_vector.resize(r_number_of_integration_points);
        if (m_T_vector.size() != r_number_of_integration_points)
            m_T_vector.resize(r_number_of_integration_points);
        if (m_T_hat_vector.size() != r_number_of_integration_points)
            m_T_hat_vector.resize(r_number_of_integration_points);
        if (m_reference_contravariant_base.size() != r_number_of_integration_points)
            m_reference_contravariant_base.resize(r_number_of_integration_points);
        if (m_n_contravariant_vector.size() != r_number_of_integration_points)
            m_n_contravariant_vector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
        
            //Compute Kinematics Reference
            KinematicVariables kinematic_variables_reference(
                r_geometry.WorkingSpaceDimension());

            CalculateKinematics(
                point_number,
                kinematic_variables_reference,shape_functions_gradients_i, ConfigurationType::Reference);           
            
            m_A_ab_covariant_vector[point_number] = kinematic_variables_reference.a_ab_covariant;

            m_dA_vector[point_number] = kinematic_variables_reference.dA;

            m_n_contravariant_vector[point_number] = kinematic_variables_reference.n_contravariant;

            CalculateTransformation(kinematic_variables_reference, m_T_vector[point_number], m_T_hat_vector[point_number], m_reference_contravariant_base[point_number]);

            // Compute Kinematics Actual
            KinematicVariables kinematic_variables(
                r_geometry.WorkingSpaceDimension());
            
            CalculateKinematics(
                point_number,
                kinematic_variables,shape_functions_gradients_i, ConfigurationType::Current);
            
            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                r_geometry, GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            //Define Prestress
            array_1d<double, 3> prestress = GetProperties()[PRESTRESS];
            double thickness = GetProperties()[THICKNESS];

            PrestresstransVariables prestresstrans_variables(3);

            CalculateTransformationmatrixPrestress(
                kinematic_variables,
                prestresstrans_variables 
            );

            array_1d<double, 3> transformed_prestress = prod(prestresstrans_variables.Tpre, prestress);
            constitutive_variables_membrane.StressVector += transformed_prestress * thickness;

            // calculate traction vectors
            array_1d<double, 3> traction_vector;

            CalculateTraction(point_number, traction_vector, kinematic_variables, constitutive_variables_membrane);

            // calculate the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
            Matrix N_curvilinear = ZeroMatrix(3, 3*number_of_nodes);

            CalculateNCurvilinear(point_number, N_curvilinear, kinematic_variables, constitutive_variables_membrane);

            // calculate first variation of traction vectors
            Matrix d_traction = ZeroMatrix(3, 3*number_of_nodes);

            CalculateDTraction(point_number, d_traction, N_curvilinear, kinematic_variables, constitutive_variables_membrane);
            
            //Compute the NURBS basis functions
            Matrix N = r_geometry.ShapeFunctionsValues();
            Matrix r_N = ZeroMatrix(3, 3*number_of_nodes);

            for (IndexType r = 0; r < number_of_nodes; r++)
            {
                r_N(0, 3 * r) = N(0, r);
                r_N(1, 3 * r + 1) = N(0, r);
                r_N(2, 3 * r + 2) = N(0, r);
            }

            //Get the displacement vectors of the previous iteration step
            Vector current_displacement = ZeroVector(mat_size);

            GetValuesVector(current_displacement);

            array_1d<double, 3> disp_vector;

            disp_vector = prod(r_N, current_displacement);
            const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);
            disp_vector -= displacement;

            //Compute the necessary products needed for the second variations of the traction vectors
            Matrix Pi = ZeroMatrix(3, 3);

            CalculateDDTractionProduct(point_number, Pi, kinematic_variables, constitutive_variables_membrane);

            array_1d<double, 3> dd_traction_product_vector;

            dd_traction_product_vector = prod(trans(Pi), disp_vector);

            // calculate second variation of traction vectors
            Matrix dd_traction = ZeroMatrix(3 * number_of_nodes, 3 * number_of_nodes);

            CalculateDDTraction(point_number, dd_traction, kinematic_variables, N_curvilinear, disp_vector, dd_traction_product_vector);

            //Penalty part & RHS
            Matrix H = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < number_of_nodes; i++)
            {
                IndexType index = 3 * i;
                //if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    H(0, index) = N(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    H(1, index + 1) = N(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    H(2, index + 2) = N(point_number, i);
            }

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            //const double determinat_jacobian = r_geometry.DeterminantOfJacobian(point_number);
            const double determinat_jacobian = norm_2(kinematic_variables_reference.a2);

            const double gammaTilde = 1.0;

            // Assembly
             if (CalculateStiffnessMatrixFlag) {

                noalias(rLeftHandSideMatrix) += (prod(trans(d_traction), H) + prod(trans(H), d_traction))
                    * integration_weight * determinat_jacobian *  -gammaTilde;

                for (IndexType i = 0; i < 3 * number_of_nodes; i++)
                {
                    for (IndexType j = 0; j < 3 * number_of_nodes; j++)
                    {
                        rLeftHandSideMatrix(i, j) += dd_traction(i, j) 
                            * integration_weight * determinat_jacobian * -gammaTilde;
                    }
                }

                noalias(rLeftHandSideMatrix) += prod(trans(H), H)
                     * integration_weight * determinat_jacobian * stabilization_parameter;
            }

            if (CalculateResidualVectorFlag) {

                Vector u(mat_size);
                for (IndexType i = 0; i < number_of_nodes; i++)
                {
                    const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * i;
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }
                
                noalias(rRightHandSideVector) -= (prod(trans(H), traction_vector)) 
                    * integration_weight * determinat_jacobian * -gammaTilde;
                noalias(rRightHandSideVector) -= (prod(trans(d_traction), disp_vector))
                    * integration_weight * determinat_jacobian * -gammaTilde;   
                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
                    * integration_weight * determinat_jacobian * stabilization_parameter;
                
            }
        }

        KRATOS_CATCH("")
    }

    void SupportNitscheCondition::CalculateQMatrix(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // Size definitions
        const SizeType number_of_nodes = r_geometry.size();

        const SizeType mat_size = 3 * (number_of_nodes);

        // Memory allocation
        if (rLeftHandSideMatrix.size1() != mat_size) {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        if (rRightHandSideVector.size() != mat_size) {
                rRightHandSideVector.resize(mat_size, false);
            }
            rRightHandSideVector = ZeroVector(mat_size);
     
        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        const IntegrationMethod integration_method = r_geometry.GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geometry.ShapeFunctionsLocalGradients(integration_method);
 
        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();
        
        // Prepare memory
        if (m_A_ab_covariant_vector.size() != r_number_of_integration_points)
            m_A_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_dA_vector.size() != r_number_of_integration_points)
            m_dA_vector.resize(r_number_of_integration_points);
        if (m_T_vector.size() != r_number_of_integration_points)
            m_T_vector.resize(r_number_of_integration_points);
        if (m_T_hat_vector.size() != r_number_of_integration_points)
            m_T_hat_vector.resize(r_number_of_integration_points);
        if (m_reference_contravariant_base.size() != r_number_of_integration_points)
            m_reference_contravariant_base.resize(r_number_of_integration_points);
        if (m_n_contravariant_vector.size() != r_number_of_integration_points)
            m_n_contravariant_vector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
        
            //Compute Kinematics Reference
            KinematicVariables kinematic_variables_reference(
                r_geometry.WorkingSpaceDimension());

            CalculateKinematics(
                point_number,
                kinematic_variables_reference,shape_functions_gradients_i, ConfigurationType::Reference);           

            m_A_ab_covariant_vector[point_number] = kinematic_variables_reference.a_ab_covariant;

            m_dA_vector[point_number] = kinematic_variables_reference.dA;

            m_n_contravariant_vector[point_number] = kinematic_variables_reference.n_contravariant;

            CalculateTransformation(kinematic_variables_reference, m_T_vector[point_number], m_T_hat_vector[point_number], m_reference_contravariant_base[point_number]);

            // Compute Kinematics Actual
            KinematicVariables kinematic_variables(
                r_geometry.WorkingSpaceDimension());

            CalculateKinematics(
                point_number,
                kinematic_variables,shape_functions_gradients_i, ConfigurationType::Current);
            
            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                r_geometry, GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            //Define Prestress
            array_1d<double, 3> prestress = GetProperties()[PRESTRESS];
            double thickness = GetProperties()[THICKNESS];

            PrestresstransVariables prestresstrans_variables(3);

            CalculateTransformationmatrixPrestress(
                kinematic_variables,
                prestresstrans_variables 
            );

            array_1d<double, 3> transformed_prestress = prod(prestresstrans_variables.Tpre, prestress);
            constitutive_variables_membrane.StressVector += transformed_prestress * thickness;

            // calculate traction vectors
            array_1d<double, 3> traction_vector;

            CalculateTraction(point_number, traction_vector, kinematic_variables, constitutive_variables_membrane);

            // calculate the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
            Matrix N_curvilinear = ZeroMatrix(3, 3*number_of_nodes);

            CalculateNCurvilinear(point_number, N_curvilinear, kinematic_variables, constitutive_variables_membrane);

            // calculate first variation of traction vectors
            Matrix d_traction = ZeroMatrix(3, 3*number_of_nodes);

            CalculateDTraction(point_number, d_traction, N_curvilinear, kinematic_variables, constitutive_variables_membrane);

            //Compute the NURBS basis functions
            Matrix N = r_geometry.ShapeFunctionsValues();
            Matrix r_N = ZeroMatrix(3, 3*number_of_nodes);

            for (IndexType r = 0; r < number_of_nodes; r++)
            {
                r_N(0, 3 * r) = N(0, r);
                r_N(1, 3 * r + 1) = N(0, r);
                r_N(2, 3 * r + 2) = N(0, r);
            }

            //Compute the necessary products needed for the second variations of the traction vectors
            Matrix Pi = ZeroMatrix(3, 3);

            CalculateDDTractionProduct(point_number, Pi, kinematic_variables, constitutive_variables_membrane);

            array_1d<double, 3> dd_traction_product_vector;

            dd_traction_product_vector = prod(trans(Pi), traction_vector);

            // calculate second variation of traction vectors

            Matrix dd_traction = ZeroMatrix(3 * number_of_nodes, 3 * number_of_nodes);

            CalculateDDTraction(point_number, dd_traction, kinematic_variables, N_curvilinear, traction_vector, dd_traction_product_vector);
           
            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = r_geometry.DeterminantOfJacobian(point_number);

            const double gammaTilde = 1.0;

            // Assembly
            noalias(rLeftHandSideMatrix) += 2 * prod(trans(d_traction), d_traction)
                * integration_weight * determinat_jacobian * gammaTilde * gammaTilde;

            for (IndexType i = 0; i < 3 * number_of_nodes; i++)
            {
                for (IndexType j = 0; j < 3 * number_of_nodes; j++)
                {
                    rLeftHandSideMatrix(i, j) += dd_traction(i, j) * integration_weight * determinat_jacobian * gammaTilde * gammaTilde;
                }
            }
        }
        KRATOS_CATCH("")
    }

    void SupportNitscheCondition::CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe)
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

    void SupportNitscheCondition::CalculateTraction(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;
        array_1d<double, 2> n_contravariant_vector;

        stress_vector_covariant = prod(m_T_hat_vector[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

        // Compute the stress components
        Matrix Palphabeta = ZeroMatrix(2, 2);
        Palphabeta(0,0) = stress_vector_covariant[0];
        Palphabeta(1,1) = stress_vector_covariant[1];
        Palphabeta(0,1) = stress_vector_covariant[2];
        Palphabeta(1,0) = Palphabeta(0,1);

        array_1d<double, 2> P_n;
        P_n = prod(Palphabeta, n_contravariant_vector);
        
        // Compute the traction vectors
        rTraction[0] = rActualKinematic.a1[0]*P_n[0] + rActualKinematic.a2[0]*P_n[1];
        rTraction[1] = rActualKinematic.a1[1]*P_n[0] + rActualKinematic.a2[1]*P_n[1];
        rTraction[2] = rActualKinematic.a1[2]*P_n[0] + rActualKinematic.a2[2]*P_n[1];
    }

    void SupportNitscheCondition::CalculateDDTractionProduct(
        IndexType IntegrationPointIndex,
        Matrix& rPi,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        array_1d<double, 2> n_contravariant_vector;

        rPi = prod(m_T_hat_vector[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.ConstitutiveMatrix);
        rPi = prod(rPi, m_T_vector[IntegrationPointIndex]);

        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

        Matrix a_u = ZeroMatrix(3, 3); //covariant times normal vector

        for (IndexType r = 0; r < 3; r++)
        {
            a_u (r, 0) = rActualKinematic.a1[r] * n_contravariant_vector[0];
            a_u (r, 1) = rActualKinematic.a2[r] * n_contravariant_vector[1];
            a_u (r, 2) = rActualKinematic.a1[r] * n_contravariant_vector[1] + rActualKinematic.a2[r] * n_contravariant_vector[0];
        }

        rPi = prod(a_u, rPi);
    }

    void SupportNitscheCondition::CalculateNCurvilinear(
        IndexType IntegrationPointIndex,
        Matrix& rNCurvilinear,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        const auto& r_geometry = GetGeometry();
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        //Compute the first variation of the Green-Lagrange strains
        Matrix dE_cartesian = ZeroMatrix(3, mat_size);
        Matrix T_patch = ZeroMatrix(3, 3);

        array_1d<double, 2> n_contravariant_vector; 

        T_patch = m_T_vector[IntegrationPointIndex];
        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            array_1d<double, 3> dE_curvilinear;
            // strain
            dE_curvilinear[0] = r_DN_De(kr, 0)*rActualKinematic.a1(dirr);
            dE_curvilinear[1] = r_DN_De(kr, 1)*rActualKinematic.a2(dirr);
            dE_curvilinear[2] = 0.5*(r_DN_De(kr, 0)*rActualKinematic.a2(dirr) + rActualKinematic.a1(dirr)*r_DN_De(kr, 1));

            dE_cartesian(0, r) = T_patch(0, 0)*dE_curvilinear[0] + T_patch(0, 1)*dE_curvilinear[1] + T_patch(0, 2)*dE_curvilinear[2];
            dE_cartesian(1, r) = T_patch(1, 0)*dE_curvilinear[0] + T_patch(1, 1)*dE_curvilinear[1] + T_patch(1, 2)*dE_curvilinear[2];
            dE_cartesian(2, r) = T_patch(2, 0)*dE_curvilinear[0] + T_patch(2, 1)*dE_curvilinear[1] + T_patch(2, 2)*dE_curvilinear[2];
        }

        //Compute the first variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
        Matrix dN_cartesian = ZeroMatrix(3, mat_size);
        dN_cartesian = prod(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix,dE_cartesian);

        //Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
        rNCurvilinear = prod(m_T_hat_vector[IntegrationPointIndex], dN_cartesian);
    }
    

    void SupportNitscheCondition::CalculateDTraction(
        IndexType IntegrationPointIndex,
        Matrix& rDTraction,
        Matrix& rNCurvilinear,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        const auto& r_geometry = GetGeometry();
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        array_1d<double, 2> n_contravariant_vector; 

        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

        //Compute the first variation of the traction vectors
        Matrix n_a = ZeroMatrix(3, 3); //normal vector times covariant current

        n_a(0, 0) = n_contravariant_vector[0]*rActualKinematic.a1[0];
        n_a(1, 0) = n_contravariant_vector[0]*rActualKinematic.a1[1];
        n_a(2, 0) = n_contravariant_vector[0]*rActualKinematic.a1[2];
        n_a(0, 1) = n_contravariant_vector[1]*rActualKinematic.a2[0];
        n_a(1, 1) = n_contravariant_vector[1]*rActualKinematic.a2[1];
        n_a(2, 1) = n_contravariant_vector[1]*rActualKinematic.a2[2];
        n_a(0, 2) = n_contravariant_vector[1]*rActualKinematic.a1[0] + n_contravariant_vector[0]*rActualKinematic.a2[0];
        n_a(1, 2) = n_contravariant_vector[1]*rActualKinematic.a1[1] + n_contravariant_vector[0]*rActualKinematic.a2[1];
        n_a(2, 2) = n_contravariant_vector[1]*rActualKinematic.a1[2] + n_contravariant_vector[0]*rActualKinematic.a2[2];

        rDTraction = prod(n_a, rNCurvilinear); 

        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;

        stress_vector_covariant = prod(m_T_hat_vector[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);

        Matrix r_DN_Dxi = ZeroMatrix(3, mat_size);
        Matrix r_DN_Deta = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < number_of_control_points; r++)
        {
            r_DN_Dxi(0, 3 * r) = r_DN_De(r, 0);
            r_DN_Dxi(1, 3 * r + 1) = r_DN_De(r, 0);
            r_DN_Dxi(2, 3 * r + 2) = r_DN_De(r, 0);

            r_DN_Deta(0, 3 * r) = r_DN_De(r, 1);
            r_DN_Deta(1, 3 * r + 1) = r_DN_De(r, 1);
            r_DN_Deta(2, 3 * r + 2) = r_DN_De(r, 1);
        }

        rDTraction += r_DN_Dxi*(n_contravariant_vector[0]*stress_vector_covariant[0]+n_contravariant_vector[1]*stress_vector_covariant[2])+
                      r_DN_Deta*(n_contravariant_vector[1]*stress_vector_covariant[1]+n_contravariant_vector[0]*stress_vector_covariant[2]);
    }

    void SupportNitscheCondition::CalculateDDTraction(
        IndexType IntegrationPointIndex,
        Matrix& rDDTraction,
        const KinematicVariables& rActualKinematic,
        Matrix& rNCurvilinear, 
        array_1d<double, 3>& rDisp,
        array_1d<double, 3>& rDDTractionProduct)
    {

        const auto& r_geometry = GetGeometry();

        const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
 
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;

        Matrix r_DN_Dxi = ZeroMatrix(3, mat_size);
        Matrix r_DN_Deta = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < number_of_control_points; r++)
        {
            r_DN_Dxi(0, 3 * r) = r_DN_De(r, 0);
            r_DN_Dxi(1, 3 * r + 1) = r_DN_De(r, 0);
            r_DN_Dxi(2, 3 * r + 2) = r_DN_De(r, 0);

            r_DN_Deta(0, 3 * r) = r_DN_De(r, 1);
            r_DN_Deta(1, 3 * r + 1) = r_DN_De(r, 1);
            r_DN_Deta(2, 3 * r + 2) = r_DN_De(r, 1);
        }

        array_1d<double, 2> n_contravariant_vector; 

        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

        Vector dNCovariant1 = ZeroVector(mat_size);
        Vector dNCovariant2 = ZeroVector(mat_size);
        Vector dNCovariant3 = ZeroVector(mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            dNCovariant1(r) = rNCurvilinear(0, r);
            dNCovariant2(r) = rNCurvilinear(1, r);
            dNCovariant3(r) = rNCurvilinear(2, r);
        }

        Matrix dispdNCovariant1 = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(1,:)
        Matrix dispdNCovariant2 = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(2,:)
        Matrix dispdNCovariant3 = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(3,:)

        for (IndexType r = 0; r < mat_size; r++)
        {
            dispdNCovariant1(0, r) = rDisp(0)*dNCovariant1(r);
            dispdNCovariant1(1, r) = rDisp(1)*dNCovariant1(r);
            dispdNCovariant1(2, r) = rDisp(2)*dNCovariant1(r);

            dispdNCovariant2(0, r) = rDisp(0)*dNCovariant2(r);
            dispdNCovariant2(1, r) = rDisp(1)*dNCovariant2(r);
            dispdNCovariant2(2, r) = rDisp(2)*dNCovariant2(r);

            dispdNCovariant3(0, r) = rDisp(0)*dNCovariant3(r);
            dispdNCovariant3(1, r) = rDisp(1)*dNCovariant3(r);
            dispdNCovariant3(2, r) = rDisp(2)*dNCovariant3(r);
        }

        rDDTraction += prod(trans(r_DN_Dxi), r_DN_Dxi)*rDDTractionProduct(0);
        rDDTraction += prod(trans(r_DN_Deta), r_DN_Deta)*rDDTractionProduct(1);
        rDDTraction += 0.5*(prod(trans(r_DN_Dxi), r_DN_Deta) + prod(trans(r_DN_Deta), r_DN_Dxi))*rDDTractionProduct(2);

        rDDTraction += prod(trans(r_DN_Dxi), dispdNCovariant1)*n_contravariant_vector(0);
        rDDTraction += prod(trans(r_DN_Deta), dispdNCovariant2)*n_contravariant_vector(1);
        rDDTraction += (prod(trans(r_DN_Dxi), dispdNCovariant3)*n_contravariant_vector(1) + prod(trans(r_DN_Deta), dispdNCovariant3)*n_contravariant_vector(0));

        rDDTraction += prod(trans(dispdNCovariant1), r_DN_Dxi)*n_contravariant_vector(0);
        rDDTraction += prod(trans(dispdNCovariant2), r_DN_Deta)*n_contravariant_vector(1);
        rDDTraction += (prod(trans(dispdNCovariant3), r_DN_Dxi)*n_contravariant_vector(1) + prod(trans(dispdNCovariant3), r_DN_Deta)*n_contravariant_vector(0));
    }

    void SupportNitscheCondition::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const auto& r_geometry = GetGeometry();

        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void SupportNitscheCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();

        const SizeType number_of_nodes= r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void SupportNitscheCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const auto r_geometry = GetGeometry();

        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


