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
//                   

// System includes

// External includes

// Project includes
#include "custom_conditions/coupling_nitsche_condition.h"

namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void CouplingNitscheCondition::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ConfigurationType& rConfiguration,
        const PatchType& rPatch
    )
    {
        IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);

        // pass/call this ShapeFunctionsLocalGradients[pnt]
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType number_of_nodes = r_geometry.size();
        Vector g1 = ZeroVector(dimension);
        Vector g2 = ZeroVector(dimension);

        Vector current_displacement_total = ZeroVector(dimension*(GetGeometry().GetGeometryPart(0).size()+GetGeometry().GetGeometryPart(1).size()));
        Vector current_displacement = ZeroVector(dimension*number_of_nodes);

        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement_total);

        if (rPatch==PatchType::Master)
        {
            for (SizeType i=0;i<dimension*number_of_nodes;++i){
                current_displacement[i] = current_displacement_total[i];
            }
        }
        else
        {
            for (SizeType i=0;i<dimension*number_of_nodes;++i){
                current_displacement[i] = current_displacement_total[i+GetGeometry().GetGeometryPart(0).size()*3];
            }
        }

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
        array_1d<double, 3> local_tangent;
        GetGeometry().GetGeometryPart(GeometryPart).Calculate(LOCAL_TANGENT, local_tangent);

        rKinematicVariables.t = local_tangent[0]*g1 + local_tangent[1]*g2;
        MathUtils<double>::CrossProduct(rKinematicVariables.n, rKinematicVariables.t/norm_2(rKinematicVariables.t), rKinematicVariables.a3);

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
    void CouplingNitscheCondition::CalculateTransformation(
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

        //Transformation matrix T from local cartesian basis to covariant basis
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

    void CouplingNitscheCondition::CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure,
        const PatchType& rPatch
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        if (rPatch==PatchType::Master)
        {
            array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector_master[IntegrationPointIndex]);
            noalias(rThisConstitutiveVariablesMembrane.StrainVector) = prod(m_T_vector_master[IntegrationPointIndex], strain_vector);

            // Constitive Matrices DMembrane
            rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
            rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
            rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

            ConstitutiveLaw::Pointer constitutive_law_master = GetProperties().GetSubProperties().front()[CONSTITUTIVE_LAW];
            constitutive_law_master->InitializeMaterial(GetProperties(), GetGeometry().GetGeometryPart(0), row(GetGeometry().GetGeometryPart(0).ShapeFunctionsValues(), IntegrationPointIndex));

            constitutive_law_master->CalculateMaterialResponse(rValues, ThisStressMeasure);
            rThisConstitutiveVariablesMembrane.ConstitutiveMatrix *= GetProperties().GetSubProperties().front()[THICKNESS];
        }
        else
        {
            array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector_slave[IntegrationPointIndex]);
            noalias(rThisConstitutiveVariablesMembrane.StrainVector) = prod(m_T_vector_slave[IntegrationPointIndex], strain_vector);

            // Constitive Matrices DMembrane
            rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
            rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
            rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

            ConstitutiveLaw::Pointer constitutive_law_slave = GetProperties().GetSubProperties().back()[CONSTITUTIVE_LAW];
            constitutive_law_slave->InitializeMaterial(GetProperties(), GetGeometry().GetGeometryPart(1), row(GetGeometry().GetGeometryPart(1).ShapeFunctionsValues(), IntegrationPointIndex));

            constitutive_law_slave->CalculateMaterialResponse(rValues, ThisStressMeasure);
            rThisConstitutiveVariablesMembrane.ConstitutiveMatrix *= GetProperties().GetSubProperties().back()[THICKNESS];
        }

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
    }

    //Prestress Transformation Matrix
    void CouplingNitscheCondition::CalculateTransformationPrestress(
        Matrix& rTransformationPrestress,
        const KinematicVariables& rActualKinematic
    )
    {
        //define base vector in reference plane
        array_1d<double, 3> t1;
        array_1d<double, 3> t2;
        
        if (Has(LOCAL_PRESTRESS_AXIS_1))
        {
            t1 = GetValue(LOCAL_PRESTRESS_AXIS_1);
            MathUtils<double>::CrossProduct(t2, rActualKinematic.a3, t1);
        }
        else if (Has(LOCAL_PRESTRESS_AXIS_1) && Has(LOCAL_PRESTRESS_AXIS_2))
        {
            t1 = GetValue(LOCAL_PRESTRESS_AXIS_1);
            t2 = GetValue(LOCAL_PRESTRESS_AXIS_2); 
        }

        array_1d<double, 3> t1_n = t1/norm_2(t1);
        array_1d<double, 3> t2_n = t2/norm_2(t2);

        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rActualKinematic.a_ab_covariant[0] * rActualKinematic.a_ab_covariant[1]
                - rActualKinematic.a_ab_covariant[2] * rActualKinematic.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] =  inv_det_g_ab * rActualKinematic.a_ab_covariant[1];
        a_ab_contravariant[1] =  inv_det_g_ab * rActualKinematic.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rActualKinematic.a_ab_covariant[2];

        array_1d<double, 3> a_contravariant_2 = rActualKinematic.a1*a_ab_contravariant[2] + rActualKinematic.a2*a_ab_contravariant[1]; 

        //local cartesian coordinates oriented along the 1st base vector in the ref. config.
        double l_a1 = norm_2(rActualKinematic.a1);
        array_1d<double, 3> e1 = rActualKinematic.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        //Transformation matrix from the projected basis T to the local cartesian basis
        double eG11 = inner_prod(e1,t1_n);
        double eG12 = inner_prod(e1,t2_n);
        double eG21 = inner_prod(e2,t1_n);
        double eG22 = inner_prod(e2,t2_n);
    
        rTransformationPrestress(0,0) = eG11*eG11;
        rTransformationPrestress(0,1) = eG12*eG12;
        rTransformationPrestress(0,2) = 2.0*eG11*eG12;

        rTransformationPrestress(1,0) = eG21*eG21;
        rTransformationPrestress(1,1) = eG22*eG22;
        rTransformationPrestress(1,2) = 2.0*eG21*eG22;

        rTransformationPrestress(2,0) = eG11*eG21;
        rTransformationPrestress(2,1) = eG12*eG22;
        rTransformationPrestress(2,2) = eG11*eG22+eG12*eG21;      
    }

    void CouplingNitscheCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        const double stabilization_parameter = GetProperties()[NITSCHE_STABILIZATION_FACTOR];

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size = 3 * (number_of_nodes_master + number_of_nodes_slave);

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
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();

        // initial determinant of jacobian 
        Vector determinant_jacobian_vector_initial(integration_points.size());
        DeterminantOfJacobianInitial(r_geometry_master, determinant_jacobian_vector_initial);

        const IntegrationMethod integration_method_master = r_geometry_master.GetDefaultIntegrationMethod();
        const IntegrationMethod integration_method_slave = r_geometry_slave.GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_master = r_geometry_master.ShapeFunctionsLocalGradients(integration_method_master);
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_slave = r_geometry_slave.ShapeFunctionsLocalGradients(integration_method_slave);

        const SizeType r_number_of_integration_points_master = r_geometry_master.IntegrationPointsNumber();
        const SizeType r_number_of_integration_points_slave = r_geometry_slave.IntegrationPointsNumber();
        
        // Prepare memory
        if (m_A_ab_covariant_vector_master.size() != r_number_of_integration_points_master)
            m_A_ab_covariant_vector_master.resize(r_number_of_integration_points_master);
        if (m_A_ab_covariant_vector_slave.size() != r_number_of_integration_points_slave)
            m_A_ab_covariant_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_dA_vector_master.size() != r_number_of_integration_points_master)
            m_dA_vector_master.resize(r_number_of_integration_points_master);
        if (m_dA_vector_slave.size() != r_number_of_integration_points_slave)
            m_dA_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_T_vector_master.size() != r_number_of_integration_points_master)
            m_T_vector_master.resize(r_number_of_integration_points_master);
        if (m_T_vector_slave.size() != r_number_of_integration_points_slave)
            m_T_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_T_hat_vector_master.size() != r_number_of_integration_points_master)
            m_T_hat_vector_master.resize(r_number_of_integration_points_master);
        if (m_T_hat_vector_slave.size() != r_number_of_integration_points_slave)
            m_T_hat_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_reference_contravariant_base_master.size() != r_number_of_integration_points_master)
            m_reference_contravariant_base_master.resize(r_number_of_integration_points_master);
        if (m_reference_contravariant_base_slave.size() != r_number_of_integration_points_slave)
            m_reference_contravariant_base_slave.resize(r_number_of_integration_points_slave);
        if (m_n_contravariant_vector_master.size() != r_number_of_integration_points_master)
            m_n_contravariant_vector_master.resize(r_number_of_integration_points_master);
        if (m_n_contravariant_vector_slave.size() != r_number_of_integration_points_slave)
            m_n_contravariant_vector_slave.resize(r_number_of_integration_points_slave);

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const Matrix& shape_functions_gradients_i_master = r_shape_functions_gradients_master[point_number];
            const Matrix& shape_functions_gradients_i_slave = r_shape_functions_gradients_slave[point_number];
        
            //Compute Kinematics Reference
            KinematicVariables kinematic_variables_reference_master(
                r_geometry_master.WorkingSpaceDimension());
            
            KinematicVariables kinematic_variables_reference_slave(
                r_geometry_slave.WorkingSpaceDimension());

            CalculateKinematics(
                point_number,
                kinematic_variables_reference_master,shape_functions_gradients_i_master, ConfigurationType::Reference, PatchType::Master);           
            
            CalculateKinematics(
                point_number,
                kinematic_variables_reference_slave,shape_functions_gradients_i_slave, ConfigurationType::Reference, PatchType::Slave);

            m_A_ab_covariant_vector_master[point_number] = kinematic_variables_reference_master.a_ab_covariant;
            m_A_ab_covariant_vector_slave[point_number] = kinematic_variables_reference_slave.a_ab_covariant;

            m_dA_vector_master[point_number] = kinematic_variables_reference_master.dA;
            m_dA_vector_slave[point_number] = kinematic_variables_reference_slave.dA;

            m_n_contravariant_vector_master[point_number] = kinematic_variables_reference_master.n_contravariant;
            m_n_contravariant_vector_slave[point_number] = kinematic_variables_reference_slave.n_contravariant;

            CalculateTransformation(kinematic_variables_reference_master, m_T_vector_master[point_number], m_T_hat_vector_master[point_number], m_reference_contravariant_base_master[point_number]);
            CalculateTransformation(kinematic_variables_reference_slave, m_T_vector_slave[point_number], m_T_hat_vector_slave[point_number], m_reference_contravariant_base_slave[point_number]);

            // Compute Kinematics 
            KinematicVariables kinematic_variables_master(
                r_geometry_master.WorkingSpaceDimension());
            KinematicVariables kinematic_variables_slave(
                r_geometry_slave.WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables_master,shape_functions_gradients_i_master, ConfigurationType::Current, PatchType::Master);
            CalculateKinematics(
                point_number,
                kinematic_variables_slave,shape_functions_gradients_i_slave, ConfigurationType::Current, PatchType::Slave);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters_master(
                r_geometry_master, GetProperties().GetSubProperties().front(), rCurrentProcessInfo);
            ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
                r_geometry_slave, GetProperties().GetSubProperties().back(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane_master(3);
            ConstitutiveVariables constitutive_variables_membrane_slave(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_master,
                constitutive_variables_membrane_master,
                constitutive_law_parameters_master,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Master);
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_slave,
                constitutive_variables_membrane_slave,
                constitutive_law_parameters_slave,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Slave);

            //Prestress component
            array_1d<double, 3> prestress_master = GetProperties().GetSubProperties().front()[PRESTRESS]*GetProperties().GetSubProperties().front()[THICKNESS];
            array_1d<double, 3> prestress_slave = GetProperties().GetSubProperties().back()[PRESTRESS]*GetProperties().GetSubProperties().back()[THICKNESS];
            array_1d<double, 3> transformed_prestress_master, transformed_prestress_slave;

            Matrix T_pre_master = ZeroMatrix(3, 3);
            Matrix T_pre_slave = ZeroMatrix(3, 3);

            if (Has(LOCAL_PRESTRESS_AXIS_1)) //for anisotropic prestress case
            {
                CalculateTransformationPrestress(T_pre_master, kinematic_variables_master);
                CalculateTransformationPrestress(T_pre_slave, kinematic_variables_slave);
                transformed_prestress_master = prod(T_pre_master, prestress_master);
                transformed_prestress_slave = prod(T_pre_slave, prestress_slave);
            }
            else //for isotropic prestress case
            {
                transformed_prestress_master = prestress_master;
                transformed_prestress_slave = prestress_slave;
            }
            
            constitutive_variables_membrane_master.StressVector += transformed_prestress_master;
            constitutive_variables_membrane_slave.StressVector += transformed_prestress_slave;

            // calculate traction vectors
            array_1d<double, 3> traction_vector_master;
            array_1d<double, 3> traction_vector_slave;

            CalculateTraction(point_number, traction_vector_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateTraction(point_number, traction_vector_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            // calculate the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
            Matrix first_variations_stress_covariant_master = ZeroMatrix(3, 3 * number_of_nodes_master);
            Matrix first_variations_stress_covariant_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

            CalculateFirstVariationStressCovariant(point_number, first_variations_stress_covariant_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateFirstVariationStressCovariant(point_number, first_variations_stress_covariant_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            // calculate first variation of traction vectors
            Matrix first_variations_traction_master = ZeroMatrix(3, 3 * number_of_nodes_master);
            Matrix first_variations_traction_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

            CalculateFirstVariationTraction(point_number, first_variations_traction_master, first_variations_stress_covariant_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateFirstVariationTraction(point_number, first_variations_traction_slave, first_variations_stress_covariant_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            Matrix first_variations_traction = ZeroMatrix(3, mat_size);
            for (SizeType i=0;i<3 * number_of_nodes_master;++i){
                first_variations_traction(0, i) = first_variations_traction_master(0, i);
                first_variations_traction(1, i) = first_variations_traction_master(1, i);
                first_variations_traction(2, i) = first_variations_traction_master(2, i);
            }
            for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
                first_variations_traction(0, i + 3 * number_of_nodes_master) = -first_variations_traction_slave(0, i);
                first_variations_traction(1, i + 3 * number_of_nodes_master) = -first_variations_traction_slave(1, i);
                first_variations_traction(2, i + 3 * number_of_nodes_master) = -first_variations_traction_slave(2, i);
            }

            //Compute the NURBS basis functions
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            Matrix r_N_master = ZeroMatrix(3, 3 * number_of_nodes_master);
            Matrix r_N_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

            for (IndexType r = 0; r < number_of_nodes_master; r++)
            {
                r_N_master(0, 3 * r) = N_master(0, r);
                r_N_master(1, 3 * r + 1) = N_master(0, r);
                r_N_master(2, 3 * r + 2) = N_master(0, r);
            }
            for (IndexType r = 0; r < number_of_nodes_slave; r++)
            {
                r_N_slave(0, 3 * r) = N_slave(0, r);
                r_N_slave(1, 3 * r + 1) = N_slave(0, r);
                r_N_slave(2, 3 * r + 2) = N_slave(0, r);
            }

            //Get the displacement vectors of the previous iteration step
            Vector current_displacement_total = ZeroVector(mat_size);
            Vector current_displacement_master = ZeroVector(3 * number_of_nodes_master);
            Vector current_displacement_slave = ZeroVector(3 * number_of_nodes_slave);

            GetValuesVector(current_displacement_total);

            for (SizeType i=0;i<3 * number_of_nodes_master;++i){
                current_displacement_master[i] = current_displacement_total[i];
            }
            for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
                current_displacement_slave[i] = current_displacement_total[i + 3 * number_of_nodes_master];
            }

            array_1d<double, 3> displacement_vector_master;
            array_1d<double, 3> displacement_vector_slave;

            displacement_vector_master = prod(r_N_master, current_displacement_master);
            displacement_vector_slave = prod(r_N_slave, current_displacement_slave);

            //Compute the necessary products needed for the second variations of the traction vectors
            Matrix Pi_master = ZeroMatrix(3, 3);
            Matrix Pi_slave = ZeroMatrix(3, 3);

            CalculateSecondVariationTractionProduct(point_number, Pi_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateSecondVariationTractionProduct(point_number, Pi_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            array_1d<double, 3> second_variations_traction_product_vector_master = prod(trans(Pi_master), displacement_vector_master);
            array_1d<double, 3> second_variations_traction_product_vector_slave = prod(trans(Pi_slave), displacement_vector_slave);
            array_1d<double, 3> second_variations_traction_product_vector_master_slave = prod(trans(Pi_slave), displacement_vector_master);
            array_1d<double, 3> second_variations_traction_product_vector_slave_master = prod(trans(Pi_master), displacement_vector_slave);

            // calculate second variation of traction vectors
            Matrix second_variations_traction_master = ZeroMatrix(3 * number_of_nodes_master, 3 * number_of_nodes_master);
            Matrix second_variations_traction_slave = ZeroMatrix(3 * number_of_nodes_slave, 3 * number_of_nodes_slave);

            CalculateSecondVariationTraction(point_number, second_variations_traction_master, kinematic_variables_master, first_variations_stress_covariant_master, displacement_vector_master, displacement_vector_slave, 
                                             second_variations_traction_product_vector_master, second_variations_traction_product_vector_slave_master, PatchType::Master);
            CalculateSecondVariationTraction(point_number, second_variations_traction_slave, kinematic_variables_slave, first_variations_stress_covariant_slave, displacement_vector_master, displacement_vector_slave, 
                                             second_variations_traction_product_vector_slave, second_variations_traction_product_vector_master_slave, PatchType::Slave);

            //Penalty part & RHS
            Matrix H = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                IndexType index = 3 * i;
                    H(0, index) = N_master(point_number, i);
                    H(1, index + 1) = N_master(point_number, i);
                    H(2, index + 2) = N_master(point_number, i);
            }
            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                IndexType index = 3 * (i + number_of_nodes_master);
                    H(0, index) = -N_slave(point_number, i);
                    H(1, index + 1) = -N_slave(point_number, i);
                    H(2, index + 2) = -N_slave(point_number, i);
            }

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinant_jacobian = determinant_jacobian_vector_initial[point_number];
            const double gammaTilde = 0.5;

            // Assembly
             if (CalculateStiffnessMatrixFlag) {

                noalias(rLeftHandSideMatrix) += (prod(trans(first_variations_traction), H) + prod(trans(H), first_variations_traction))
                    * integration_weight * determinant_jacobian * -gammaTilde;

                for (IndexType i = 0; i < 3 * number_of_nodes_master; i++)
                {
                    for (IndexType j = 0; j < 3 * number_of_nodes_master; j++)
                    {
                        rLeftHandSideMatrix(i, j) += second_variations_traction_master(i, j) * integration_weight * determinant_jacobian * -gammaTilde;
                    }
                }

                for (IndexType i = 0; i < 3 * number_of_nodes_slave; i++)
                {
                    for (IndexType j = 0; j < 3 * number_of_nodes_slave; j++)
                    {
                        rLeftHandSideMatrix(i + 3 * number_of_nodes_master, j + 3 * number_of_nodes_master) += second_variations_traction_slave(i, j) * integration_weight * determinant_jacobian * -gammaTilde;
                    }
                }

                noalias(rLeftHandSideMatrix) += prod(trans(H), H)
                    * integration_weight * determinant_jacobian * stabilization_parameter;
             }

            if (CalculateResidualVectorFlag) {

                Vector u(mat_size);
                for (IndexType i = 0; i < number_of_nodes_master; i++)
                {
                    const array_1d<double, 3> disp = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * i;
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }
                for (IndexType i = 0; i < number_of_nodes_slave; i++)
                {
                    const array_1d<double, 3> disp = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * (i + number_of_nodes_master);
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }
                
                noalias(rRightHandSideVector) -= (prod(trans(H), traction_vector_master) - prod(trans(H), traction_vector_slave))
                    * integration_weight * determinant_jacobian * -gammaTilde;
                noalias(rRightHandSideVector) -= (prod(trans(first_variations_traction), displacement_vector_master) - prod(trans(first_variations_traction), displacement_vector_slave))
                    * integration_weight * determinant_jacobian * -gammaTilde;   
                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
                    * integration_weight * determinant_jacobian * stabilization_parameter;
            }
        }
        KRATOS_CATCH("")
    }

    void CouplingNitscheCondition::DeterminantOfJacobianInitial(
        const GeometryType& rGeometry,
        Vector& rDeterminantOfJacobian)
    {
        const IndexType nb_integration_points = rGeometry.IntegrationPointsNumber();
        if (rDeterminantOfJacobian.size() != nb_integration_points) {
            rDeterminantOfJacobian.resize(nb_integration_points, false);
        }

        const SizeType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const SizeType local_space_dimension = rGeometry.LocalSpaceDimension();
        const SizeType nb_nodes = rGeometry.PointsNumber();

        Matrix J = ZeroMatrix(working_space_dimension, local_space_dimension);
        for (IndexType pnt = 0; pnt < nb_integration_points; pnt++)
        {
            const Matrix& r_DN_De = rGeometry.ShapeFunctionsLocalGradients()[pnt];
            J.clear();
            for (IndexType i = 0; i < nb_nodes; ++i) {
                const array_1d<double, 3>& r_coordinates = rGeometry[i].GetInitialPosition();
                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    const double value = r_coordinates[k];
                    for (IndexType m = 0; m < local_space_dimension; ++m) {
                        J(k, m) += value * r_DN_De(i, m);
                    }
                }
            }

            //Compute the tangent and  the normal to the boundary vector
            array_1d<double, 3> local_tangent;
            GetGeometry().GetGeometryPart(0).Calculate(LOCAL_TANGENT, local_tangent);

            array_1d<double, 3> a_1 = column(J, 0);
            array_1d<double, 3> a_2 = column(J, 1);

            rDeterminantOfJacobian[pnt] = norm_2(a_1 * local_tangent[0] + a_2 * local_tangent[1]);
        }
    }

    void CouplingNitscheCondition::CalculateNitscheStabilizationMatrix(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size = 3 * (number_of_nodes_master + number_of_nodes_slave);

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
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();

        // initial determinant of jacobian 
        Vector determinant_jacobian_vector_initial(integration_points.size());
        DeterminantOfJacobianInitial(r_geometry_master, determinant_jacobian_vector_initial);

        const IntegrationMethod integration_method_master = r_geometry_master.GetDefaultIntegrationMethod();
        const IntegrationMethod integration_method_slave = r_geometry_slave.GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_master = r_geometry_master.ShapeFunctionsLocalGradients(integration_method_master);
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_slave = r_geometry_slave.ShapeFunctionsLocalGradients(integration_method_slave);

        const SizeType r_number_of_integration_points_master = r_geometry_master.IntegrationPointsNumber();
        const SizeType r_number_of_integration_points_slave = r_geometry_slave.IntegrationPointsNumber();
        
        // Prepare memory
        if (m_A_ab_covariant_vector_master.size() != r_number_of_integration_points_master)
            m_A_ab_covariant_vector_master.resize(r_number_of_integration_points_master);
        if (m_A_ab_covariant_vector_slave.size() != r_number_of_integration_points_slave)
            m_A_ab_covariant_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_dA_vector_master.size() != r_number_of_integration_points_master)
            m_dA_vector_master.resize(r_number_of_integration_points_master);
        if (m_dA_vector_slave.size() != r_number_of_integration_points_slave)
            m_dA_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_T_vector_master.size() != r_number_of_integration_points_master)
            m_T_vector_master.resize(r_number_of_integration_points_master);
        if (m_T_vector_slave.size() != r_number_of_integration_points_slave)
            m_T_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_T_hat_vector_master.size() != r_number_of_integration_points_master)
            m_T_hat_vector_master.resize(r_number_of_integration_points_master);
        if (m_T_hat_vector_slave.size() != r_number_of_integration_points_slave)
            m_T_hat_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_reference_contravariant_base_master.size() != r_number_of_integration_points_master)
            m_reference_contravariant_base_master.resize(r_number_of_integration_points_master);
        if (m_reference_contravariant_base_slave.size() != r_number_of_integration_points_slave)
            m_reference_contravariant_base_slave.resize(r_number_of_integration_points_slave);
        if (m_n_contravariant_vector_master.size() != r_number_of_integration_points_master)
            m_n_contravariant_vector_master.resize(r_number_of_integration_points_master);
        if (m_n_contravariant_vector_slave.size() != r_number_of_integration_points_slave)
            m_n_contravariant_vector_slave.resize(r_number_of_integration_points_slave);

        //check wheter the size of the element larger than tolerance or not
        array_1d<double, 3> characteristic_length_master;
        array_1d<double, 3> characteristic_length_slave;
        r_geometry_master.GetGeometryParent(0).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->Calculate(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_master);
        r_geometry_slave.GetGeometryParent(0).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->Calculate(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_slave);

        const double tol_basic = 0.01;
        double characteristic_length = std::max(norm_2(characteristic_length_master), norm_2(characteristic_length_slave));
        double characteristic_area = characteristic_length*characteristic_length;
        double tol_surface_normal = tol_basic*characteristic_area;

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const Matrix& shape_functions_gradients_i_master = r_shape_functions_gradients_master[point_number];
            const Matrix& shape_functions_gradients_i_slave = r_shape_functions_gradients_slave[point_number];
        
            //Compute Kinematics Reference
            KinematicVariables kinematic_variables_reference_master(
                r_geometry_master.WorkingSpaceDimension());
            
            KinematicVariables kinematic_variables_reference_slave(
                r_geometry_slave.WorkingSpaceDimension());

            CalculateKinematics(
                point_number,
                kinematic_variables_reference_master,shape_functions_gradients_i_master, ConfigurationType::Reference, PatchType::Master);           
            
            CalculateKinematics(
                point_number,
                kinematic_variables_reference_slave,shape_functions_gradients_i_slave, ConfigurationType::Reference, PatchType::Slave);

            m_A_ab_covariant_vector_master[point_number] = kinematic_variables_reference_master.a_ab_covariant;
            m_A_ab_covariant_vector_slave[point_number] = kinematic_variables_reference_slave.a_ab_covariant;

            m_dA_vector_master[point_number] = kinematic_variables_reference_master.dA;
            m_dA_vector_slave[point_number] = kinematic_variables_reference_slave.dA;

            m_n_contravariant_vector_master[point_number] = kinematic_variables_reference_master.n_contravariant;
            m_n_contravariant_vector_slave[point_number] = kinematic_variables_reference_slave.n_contravariant;

            CalculateTransformation(kinematic_variables_reference_master, m_T_vector_master[point_number], m_T_hat_vector_master[point_number], m_reference_contravariant_base_master[point_number]);
            CalculateTransformation(kinematic_variables_reference_slave, m_T_vector_slave[point_number], m_T_hat_vector_slave[point_number], m_reference_contravariant_base_slave[point_number]);

            // Compute Kinematics Actual
            KinematicVariables kinematic_variables_master(
                r_geometry_master.WorkingSpaceDimension());
            KinematicVariables kinematic_variables_slave(
                r_geometry_slave.WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables_master,shape_functions_gradients_i_master, ConfigurationType::Current, PatchType::Master);
            CalculateKinematics(
                point_number,
                kinematic_variables_slave,shape_functions_gradients_i_slave, ConfigurationType::Current, PatchType::Slave);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters_master(
                r_geometry_master, GetProperties().GetSubProperties().front(), rCurrentProcessInfo);
            ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
                r_geometry_slave, GetProperties().GetSubProperties().back(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane_master(3);
            ConstitutiveVariables constitutive_variables_membrane_slave(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_master,
                constitutive_variables_membrane_master,
                constitutive_law_parameters_master,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Master);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_slave,
                constitutive_variables_membrane_slave,
                constitutive_law_parameters_slave,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Slave);

            //Prestress component
            array_1d<double, 3> prestress_master = GetProperties().GetSubProperties().front()[PRESTRESS]*GetProperties().GetSubProperties().front()[THICKNESS];
            array_1d<double, 3> prestress_slave = GetProperties().GetSubProperties().back()[PRESTRESS]*GetProperties().GetSubProperties().back()[THICKNESS];
            array_1d<double, 3> transformed_prestress_master, transformed_prestress_slave;

            Matrix T_pre_master = ZeroMatrix(3, 3);
            Matrix T_pre_slave = ZeroMatrix(3, 3);

            if (Has(LOCAL_PRESTRESS_AXIS_1)) //for anisotropic prestress case
            {
                CalculateTransformationPrestress(T_pre_master, kinematic_variables_master);
                CalculateTransformationPrestress(T_pre_slave, kinematic_variables_slave);
                transformed_prestress_master = prod(T_pre_master, prestress_master);
                transformed_prestress_slave = prod(T_pre_slave, prestress_slave);
            }
            else //for isotropic prestress case
            {
                transformed_prestress_master = prestress_master;
                transformed_prestress_slave = prestress_slave;
            }
            
            constitutive_variables_membrane_master.StressVector += transformed_prestress_master;
            constitutive_variables_membrane_slave.StressVector += transformed_prestress_slave;

            // calculate traction vectors
            array_1d<double, 3> traction_vector_master;
            array_1d<double, 3> traction_vector_slave;

            CalculateTraction(point_number, traction_vector_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateTraction(point_number, traction_vector_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            // calculate the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
            Matrix first_variations_stress_covariant_master = ZeroMatrix(3, 3*number_of_nodes_master);
            Matrix first_variations_stress_covariant_slave = ZeroMatrix(3, 3*number_of_nodes_slave);

            CalculateFirstVariationStressCovariant(point_number, first_variations_stress_covariant_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateFirstVariationStressCovariant(point_number, first_variations_stress_covariant_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            // calculate first variation of traction vectors
            Matrix first_variations_traction_master = ZeroMatrix(3, 3*number_of_nodes_master);
            Matrix first_variations_traction_slave = ZeroMatrix(3, 3*number_of_nodes_slave);

            CalculateFirstVariationTraction(point_number, first_variations_traction_master, first_variations_stress_covariant_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateFirstVariationTraction(point_number, first_variations_traction_slave, first_variations_stress_covariant_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            Matrix first_variations_traction = ZeroMatrix(3, mat_size);
            for (SizeType i=0;i<3 * number_of_nodes_master;++i){
                first_variations_traction(0, i) = first_variations_traction_master(0, i);
                first_variations_traction(1, i) = first_variations_traction_master(1, i);
                first_variations_traction(2, i) = first_variations_traction_master(2, i);
            }
            for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
                first_variations_traction(0, i + 3 * number_of_nodes_master) = -first_variations_traction_slave(0, i);
                first_variations_traction(1, i + 3 * number_of_nodes_master) = -first_variations_traction_slave(1, i);
                first_variations_traction(2, i + 3 * number_of_nodes_master) = -first_variations_traction_slave(2, i);
            }

            //Compute the NURBS basis functions
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            Matrix r_N_master = ZeroMatrix(3, 3*number_of_nodes_master);
            Matrix r_N_slave = ZeroMatrix(3, 3*number_of_nodes_slave);

            for (IndexType r = 0; r < number_of_nodes_master; r++)
            {
                r_N_master(0, 3 * r) = N_master(0, r);
                r_N_master(1, 3 * r + 1) = N_master(0, r);
                r_N_master(2, 3 * r + 2) = N_master(0, r);
            }

            for (IndexType r = 0; r < number_of_nodes_slave; r++)
            {
                r_N_slave(0, 3 * r) = N_slave(0, r);
                r_N_slave(1, 3 * r + 1) = N_slave(0, r);
                r_N_slave(2, 3 * r + 2) = N_slave(0, r);
            }

            //Compute the necessary products needed for the second variations of the traction vectors
            Matrix Pi_master = ZeroMatrix(3, 3);
            Matrix Pi_slave = ZeroMatrix(3, 3);

            CalculateSecondVariationTractionProduct(point_number, Pi_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateSecondVariationTractionProduct(point_number, Pi_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            array_1d<double, 3> second_variations_traction_product_vector_master = prod(trans(Pi_master), traction_vector_master);
            array_1d<double, 3> second_variations_traction_product_vector_slave = prod(trans(Pi_slave), traction_vector_slave);
            array_1d<double, 3> second_variations_traction_product_vector_master_slave = prod(trans(Pi_slave), traction_vector_master);
            array_1d<double, 3> second_variations_traction_product_vector_slave_master = prod(trans(Pi_master), traction_vector_slave);

            // calculate second variation of traction vectors
            Matrix second_variations_traction_master = ZeroMatrix(3 * number_of_nodes_master, 3 * number_of_nodes_master);
            Matrix second_variations_traction_slave = ZeroMatrix(3 * number_of_nodes_slave, 3 * number_of_nodes_slave);

            if(norm_2(kinematic_variables_reference_master.a3_tilde) > tol_surface_normal && norm_2(kinematic_variables_reference_slave.a3_tilde) > tol_surface_normal)
            {
                CalculateSecondVariationTraction(point_number, second_variations_traction_master, kinematic_variables_master, first_variations_stress_covariant_master, traction_vector_master, traction_vector_slave, 
                                                 second_variations_traction_product_vector_master, second_variations_traction_product_vector_slave_master, PatchType::Master);
                CalculateSecondVariationTraction(point_number, second_variations_traction_slave, kinematic_variables_slave, first_variations_stress_covariant_slave, traction_vector_master, traction_vector_slave, 
                                                 second_variations_traction_product_vector_slave, second_variations_traction_product_vector_master_slave, PatchType::Slave);
            }
            
            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinant_jacobian = determinant_jacobian_vector_initial[point_number];
            const double gammaTilde = 0.5;

            // Assembly
            noalias(rLeftHandSideMatrix) += 2 * prod(trans(first_variations_traction), first_variations_traction)
                * integration_weight * determinant_jacobian * gammaTilde * gammaTilde;
        

            for (IndexType i = 0; i < 3 * number_of_nodes_master; i++)
            {
                for (IndexType j = 0; j < 3 * number_of_nodes_master; j++)
                {
                    rLeftHandSideMatrix(i, j) += second_variations_traction_master(i, j) * integration_weight * determinant_jacobian * gammaTilde * gammaTilde;
                }
            }

            for (IndexType i = 0; i < 3 * number_of_nodes_slave; i++)
            {
                for (IndexType j = 0; j < 3 * number_of_nodes_slave; j++)
                {
                    rLeftHandSideMatrix(i + 3 * number_of_nodes_master, j + 3 * number_of_nodes_master) += second_variations_traction_slave(i, j) * integration_weight * determinant_jacobian * gammaTilde * gammaTilde;
                }
            }

            if(norm_2(kinematic_variables_reference_master.a3_tilde) < tol_surface_normal || norm_2(kinematic_variables_reference_slave.a3_tilde) < tol_surface_normal)
            {
                for (IndexType i = 0; i < 3 * number_of_nodes_master; i++)
                {
                    for (IndexType j = 0; j < 3 * number_of_nodes_master; j++)
                    {
                        rLeftHandSideMatrix(i, j) *= 0.0;
                    }
                }

                for (IndexType i = 0; i < 3 * number_of_nodes_slave; i++)
                {
                    for (IndexType j = 0; j < 3 * number_of_nodes_slave; j++)
                    {
                        rLeftHandSideMatrix(i + 3 * number_of_nodes_master, j + 3 * number_of_nodes_master) *= 0.0;
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    void CouplingNitscheCondition::CalculateTraction(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch)
    {
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;
        array_1d<double, 2> n_contravariant_vector;

        if (rPatch==PatchType::Master)
        {
            stress_vector_covariant = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
            stress_vector_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
            n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        }

        // Compute the stress components
        Matrix Palphabeta = ZeroMatrix(2, 2);
        Palphabeta(0,0) = stress_vector_covariant[0];
        Palphabeta(1,1) = stress_vector_covariant[1];
        Palphabeta(0,1) = stress_vector_covariant[2];
        Palphabeta(1,0) = Palphabeta(0,1);
        
        // Compute the traction vectors
        rTraction[0] = rActualKinematic.a1[0]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                     + rActualKinematic.a2[0]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
        rTraction[1] = rActualKinematic.a1[1]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                     + rActualKinematic.a2[1]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
        rTraction[2] = rActualKinematic.a1[2]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                     + rActualKinematic.a2[2]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
    }

    void CouplingNitscheCondition::CalculateFirstVariationStressCovariant(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch)
    {
        IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        //Compute the first variation of the Green-Lagrange strains
        Matrix dE_cartesian = ZeroMatrix(3, mat_size);
        Matrix T_patch = ZeroMatrix(3, 3);

        if (rPatch==PatchType::Master)
        {
            T_patch = m_T_vector_master[IntegrationPointIndex];
        }
        else
        {
            T_patch = m_T_vector_slave[IntegrationPointIndex];
        }

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
        Matrix first_variations_stress_cartesian = ZeroMatrix(3, mat_size);
        first_variations_stress_cartesian = prod(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix,dE_cartesian);

        //Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
        if (rPatch==PatchType::Master)
        {
            rFirstVariationStressCovariant = prod(m_T_hat_vector_master[IntegrationPointIndex], first_variations_stress_cartesian);
        }
        else
        {
            rFirstVariationStressCovariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], first_variations_stress_cartesian);
        }
    }
    
    void CouplingNitscheCondition::CalculateFirstVariationTraction(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationTraction,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch)
    {
        IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        //get the normal vector
        array_1d<double, 2> n_contravariant_vector; 
        if (rPatch==PatchType::Master)
        {
            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
            n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        }

        //Compute the first variation of the traction vectors:
        //1. normal vector * derivative stress covariant
        
        // normal vector * covariant base vector
        Matrix n_a = ZeroMatrix(3, 3); 

        for (IndexType r = 0; r < 3; r++)
        {
            n_a (r, 0) = rActualKinematic.a1[r] * n_contravariant_vector[0];
            n_a (r, 1) = rActualKinematic.a2[r] * n_contravariant_vector[1];
            n_a (r, 2) = rActualKinematic.a1[r] * n_contravariant_vector[1] + rActualKinematic.a2[r] * n_contravariant_vector[0];
        }

        rFirstVariationTraction = prod(n_a, rFirstVariationStressCovariant);

        //2. derivative normal vector * stress covariant
        
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;
        if (rPatch==PatchType::Master)
        {
            stress_vector_covariant = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
        }
        else
        {
            stress_vector_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
        }

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

        rFirstVariationTraction += r_DN_Dxi*(n_contravariant_vector[0]*stress_vector_covariant[0] + n_contravariant_vector[1]*stress_vector_covariant[2])+
                                   r_DN_Deta*(n_contravariant_vector[1]*stress_vector_covariant[1] + n_contravariant_vector[0]*stress_vector_covariant[2]);
    }

    void CouplingNitscheCondition::CalculateSecondVariationTractionProduct(
        IndexType IntegrationPointIndex,
        Matrix& rPi,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch)
    {
        array_1d<double, 2> n_contravariant_vector;

        if (rPatch==PatchType::Master)
        {
            rPi = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.ConstitutiveMatrix);
            rPi = prod(rPi, m_T_vector_master[IntegrationPointIndex]);

            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
            rPi = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.ConstitutiveMatrix);
            rPi = prod(rPi, m_T_vector_slave[IntegrationPointIndex]);

            n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        }

        // normal vector * covariant base vector
        Matrix n_a = ZeroMatrix(3, 3);

        for (IndexType r = 0; r < 3; r++)
        {
            n_a (r, 0) = rActualKinematic.a1[r] * n_contravariant_vector[0];
            n_a (r, 1) = rActualKinematic.a2[r] * n_contravariant_vector[1];
            n_a (r, 2) = rActualKinematic.a1[r] * n_contravariant_vector[1] + rActualKinematic.a2[r] * n_contravariant_vector[0];
        }

        rPi = prod(n_a, rPi);
    }

    void CouplingNitscheCondition::CalculateSecondVariationTraction(
        IndexType IntegrationPointIndex,
        Matrix& rSecondVariationTraction,
        const KinematicVariables& rActualKinematic,
        Matrix& rFirstVariationStressCovariant, 
        array_1d<double, 3>& rDisplacementMaster,
        array_1d<double, 3>& rDisplacementSlave,
        array_1d<double, 3>& rSecondVariationTractionProduct,
        array_1d<double, 3>& rSecondVariationTractionProductMasterSlave,
        const PatchType& rPatch)
    {
        IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);

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

        if (rPatch==PatchType::Master)
        {
            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
            n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        }

        Vector first_variations_stress_covariant_1 = ZeroVector(mat_size);
        Vector first_variations_stress_covariant_2 = ZeroVector(mat_size);
        Vector first_variations_stress_covariant_3 = ZeroVector(mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            first_variations_stress_covariant_1(r) = rFirstVariationStressCovariant(0, r);
            first_variations_stress_covariant_2(r) = rFirstVariationStressCovariant(1, r);
            first_variations_stress_covariant_3(r) = rFirstVariationStressCovariant(2, r);
        }

        // displacement * first variation stress covariant
        Matrix displacement_dNCovariant1_master = ZeroMatrix(3, mat_size); 
        Matrix displacement_dNCovariant2_master = ZeroMatrix(3, mat_size); 
        Matrix displacement_dNCovariant3_master = ZeroMatrix(3, mat_size); 

        Matrix displacement_dNCovariant1_slave = ZeroMatrix(3, mat_size); 
        Matrix displacement_dNCovariant2_slave = ZeroMatrix(3, mat_size); 
        Matrix displacement_dNCovariant3_slave = ZeroMatrix(3, mat_size); 

        for (IndexType r = 0; r < mat_size; r++)
        {
            displacement_dNCovariant1_master(0, r) = rDisplacementMaster(0)*first_variations_stress_covariant_1(r);
            displacement_dNCovariant1_master(1, r) = rDisplacementMaster(1)*first_variations_stress_covariant_1(r);
            displacement_dNCovariant1_master(2, r) = rDisplacementMaster(2)*first_variations_stress_covariant_1(r);

            displacement_dNCovariant2_master(0, r) = rDisplacementMaster(0)*first_variations_stress_covariant_2(r);
            displacement_dNCovariant2_master(1, r) = rDisplacementMaster(1)*first_variations_stress_covariant_2(r);
            displacement_dNCovariant2_master(2, r) = rDisplacementMaster(2)*first_variations_stress_covariant_2(r);

            displacement_dNCovariant3_master(0, r) = rDisplacementMaster(0)*first_variations_stress_covariant_3(r);
            displacement_dNCovariant3_master(1, r) = rDisplacementMaster(1)*first_variations_stress_covariant_3(r);
            displacement_dNCovariant3_master(2, r) = rDisplacementMaster(2)*first_variations_stress_covariant_3(r);

            displacement_dNCovariant1_slave(0, r) = rDisplacementSlave(0)*first_variations_stress_covariant_1(r);
            displacement_dNCovariant1_slave(1, r) = rDisplacementSlave(1)*first_variations_stress_covariant_1(r);
            displacement_dNCovariant1_slave(2, r) = rDisplacementSlave(2)*first_variations_stress_covariant_1(r);

            displacement_dNCovariant2_slave(0, r) = rDisplacementSlave(0)*first_variations_stress_covariant_2(r);
            displacement_dNCovariant2_slave(1, r) = rDisplacementSlave(1)*first_variations_stress_covariant_2(r);
            displacement_dNCovariant2_slave(2, r) = rDisplacementSlave(2)*first_variations_stress_covariant_2(r);

            displacement_dNCovariant3_slave(0, r) = rDisplacementSlave(0)*first_variations_stress_covariant_3(r);
            displacement_dNCovariant3_slave(1, r) = rDisplacementSlave(1)*first_variations_stress_covariant_3(r);
            displacement_dNCovariant3_slave(2, r) = rDisplacementSlave(2)*first_variations_stress_covariant_3(r);
        }

        if (rPatch==PatchType::Slave)
        {
            displacement_dNCovariant1_master *= -1.0;
            displacement_dNCovariant2_master *= -1.0;
            displacement_dNCovariant3_master *= -1.0;

            displacement_dNCovariant1_slave *= -1.0;
            displacement_dNCovariant2_slave *= -1.0;
            displacement_dNCovariant3_slave *= -1.0;
        }

        rSecondVariationTraction += prod(trans(r_DN_Dxi), r_DN_Dxi)*rSecondVariationTractionProduct(0);
        rSecondVariationTraction += prod(trans(r_DN_Deta), r_DN_Deta)*rSecondVariationTractionProduct(1);
        rSecondVariationTraction += 0.5*(prod(trans(r_DN_Dxi), r_DN_Deta) + prod(trans(r_DN_Deta), r_DN_Dxi))*rSecondVariationTractionProduct(2);

        rSecondVariationTraction += prod(trans(r_DN_Dxi), displacement_dNCovariant1_master)*n_contravariant_vector(0);
        rSecondVariationTraction += prod(trans(r_DN_Deta), displacement_dNCovariant2_master)*n_contravariant_vector(1);
        rSecondVariationTraction += (prod(trans(r_DN_Dxi), displacement_dNCovariant3_master)*n_contravariant_vector(1) + prod(trans(r_DN_Deta), displacement_dNCovariant3_master)*n_contravariant_vector(0));

        rSecondVariationTraction += prod(trans(displacement_dNCovariant1_master), r_DN_Dxi)*n_contravariant_vector(0);
        rSecondVariationTraction += prod(trans(displacement_dNCovariant2_master), r_DN_Deta)*n_contravariant_vector(1);
        rSecondVariationTraction += (prod(trans(displacement_dNCovariant3_master), r_DN_Dxi)*n_contravariant_vector(1) + prod(trans(displacement_dNCovariant3_master), r_DN_Deta)*n_contravariant_vector(0));

        rSecondVariationTraction -= prod(trans(r_DN_Dxi), r_DN_Dxi)*rSecondVariationTractionProductMasterSlave(0);
        rSecondVariationTraction -= prod(trans(r_DN_Deta), r_DN_Deta)*rSecondVariationTractionProductMasterSlave(1);
        rSecondVariationTraction -= 0.5*(prod(trans(r_DN_Dxi), r_DN_Deta) + prod(trans(r_DN_Deta), r_DN_Dxi))*rSecondVariationTractionProductMasterSlave(2);

        rSecondVariationTraction -= prod(trans(r_DN_Dxi), displacement_dNCovariant1_slave)*n_contravariant_vector(0);
        rSecondVariationTraction -= prod(trans(r_DN_Deta), displacement_dNCovariant2_slave)*n_contravariant_vector(1);
        rSecondVariationTraction -= (prod(trans(r_DN_Dxi), displacement_dNCovariant3_slave)*n_contravariant_vector(1) + prod(trans(r_DN_Deta), displacement_dNCovariant3_slave)*n_contravariant_vector(0));

        rSecondVariationTraction -= prod(trans(displacement_dNCovariant1_slave), r_DN_Dxi)*n_contravariant_vector(0);
        rSecondVariationTraction -= prod(trans(displacement_dNCovariant2_slave), r_DN_Deta)*n_contravariant_vector(1);
        rSecondVariationTraction -= (prod(trans(displacement_dNCovariant3_slave), r_DN_Dxi)*n_contravariant_vector(1) + prod(trans(displacement_dNCovariant3_slave), r_DN_Deta)*n_contravariant_vector(0));
    }

    void CouplingNitscheCondition::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_control_points_master = r_geometry_master.size();
        const SizeType number_of_control_points_slave = r_geometry_slave.size();
        const SizeType mat_size = (number_of_control_points_master + number_of_control_points_slave) * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points_master; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }

        for (IndexType i = 0; i < number_of_control_points_slave; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = 3 * (i + number_of_control_points_master) ;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void CouplingNitscheCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        if (rResult.size() != 3 * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(3 * (number_of_nodes_master + number_of_nodes_slave), false);

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry_master[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = 3 * (i + number_of_nodes_master);
            const auto& r_node = r_geometry_slave[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void CouplingNitscheCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * (number_of_nodes_master + number_of_nodes_slave));

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos



