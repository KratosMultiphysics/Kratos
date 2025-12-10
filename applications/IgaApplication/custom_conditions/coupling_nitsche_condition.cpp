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

        if (Is(IgaFlags::FIX_ROTATION_X))
        {
            Matrix H = ZeroMatrix(3, 3);
            CalculateHessian(H, r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod()), rPatch);

            rKinematicVariables.b_ab_covariant[0] = H(0, 0) * rKinematicVariables.a3[0] + H(1, 0) * rKinematicVariables.a3[1] + H(2, 0) * rKinematicVariables.a3[2];
            rKinematicVariables.b_ab_covariant[1] = H(0, 1) * rKinematicVariables.a3[0] + H(1, 1) * rKinematicVariables.a3[1] + H(2, 1) * rKinematicVariables.a3[2];
            rKinematicVariables.b_ab_covariant[2] = H(0, 2) * rKinematicVariables.a3[0] + H(1, 2) * rKinematicVariables.a3[1] + H(2, 2) * rKinematicVariables.a3[2];
        }
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
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
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

            if (Is(IgaFlags::FIX_ROTATION_X))
            {
                array_1d<double, 3> curvature_vector = rActualKinematic.b_ab_covariant - m_B_ab_covariant_vector_master[IntegrationPointIndex];
                noalias(rThisConstitutiveVariablesCurvature.StrainVector) = prod(m_T_vector_master[IntegrationPointIndex], curvature_vector);
                noalias(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix) = rThisConstitutiveVariablesMembrane.ConstitutiveMatrix * (pow(GetProperties().GetSubProperties().front()[THICKNESS], 2) / 12);
            }
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

            if (Is(IgaFlags::FIX_ROTATION_X))
            {
                array_1d<double, 3> curvature_vector = rActualKinematic.b_ab_covariant - m_B_ab_covariant_vector_slave[IntegrationPointIndex];
                noalias(rThisConstitutiveVariablesCurvature.StrainVector) = prod(m_T_vector_slave[IntegrationPointIndex], curvature_vector);
                noalias(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix) = rThisConstitutiveVariablesMembrane.ConstitutiveMatrix * (pow(GetProperties().GetSubProperties().back()[THICKNESS], 2) / 12);
            }
        }

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
        
        if (Is(IgaFlags::FIX_ROTATION_X))
        {
            noalias(rThisConstitutiveVariablesCurvature.StressVector) = prod(
                trans(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix), rThisConstitutiveVariablesCurvature.StrainVector);
        }
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
        double stabilization_rotation_parameter;

        if (Is(IgaFlags::FIX_ROTATION_X))
        {
            stabilization_rotation_parameter = GetProperties()[NITSCHE_STABILIZATION_ROTATION_FACTOR];
        }

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
        if (Is(IgaFlags::FIX_ROTATION_X))
        {
            if (m_B_ab_covariant_vector_master.size() != r_number_of_integration_points_master)
                m_B_ab_covariant_vector_master.resize(r_number_of_integration_points_master);
            if (m_B_ab_covariant_vector_slave.size() != r_number_of_integration_points_slave)
                m_B_ab_covariant_vector_slave.resize(r_number_of_integration_points_slave);
        }

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

            if (Is(IgaFlags::FIX_ROTATION_X))
            {
                m_B_ab_covariant_vector_master[point_number] = kinematic_variables_reference_master.b_ab_covariant;
                m_B_ab_covariant_vector_slave[point_number] = kinematic_variables_reference_slave.b_ab_covariant;
            }

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
            ConstitutiveVariables constitutive_variables_curvature_master(3);
            ConstitutiveVariables constitutive_variables_curvature_slave(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_master,
                constitutive_variables_membrane_master,
                constitutive_variables_curvature_master,
                constitutive_law_parameters_master,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Master);
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_slave,
                constitutive_variables_membrane_slave,
                constitutive_variables_curvature_slave,
                constitutive_law_parameters_slave,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Slave);

            //Prestress component
            array_1d<double, 3> prestress_master = ZeroVector(3);
            array_1d<double, 3> prestress_slave = ZeroVector(3);

            if(this->GetProperties().Has(PRESTRESS))
            {
                prestress_master = GetProperties().GetSubProperties().front()[PRESTRESS]*GetProperties().GetSubProperties().front()[THICKNESS];
                prestress_slave = GetProperties().GetSubProperties().back()[PRESTRESS]*GetProperties().GetSubProperties().back()[THICKNESS];
            }
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
            double determinant_jacobian = determinant_jacobian_vector_initial[point_number];
            const double gammaTilde = 0.5;

            // Rotation coupling
            if (Is(IgaFlags::FIX_ROTATION_X))
            {
                //Get the rotation vectors of the previous iteration step
                array_1d<double, 3> rotation_vector_master;
                array_1d<double, 3> rotation_vector_slave;

                Vector phi_r = ZeroVector(mat_size);
                Matrix phi_rs = ZeroMatrix(mat_size, mat_size);
                array_1d<double, 2> diff_phi;

                CalculateRotationalShapeFunctions(point_number, phi_r, phi_rs, diff_phi);

                //Calculate moment vectors
                array_1d<double, 3> moment_vector_master;
                array_1d<double, 3> moment_vector_slave;

                CalculateMoment(point_number, moment_vector_master, kinematic_variables_master, constitutive_variables_curvature_master, PatchType::Master);
                CalculateMoment(point_number, moment_vector_slave, kinematic_variables_slave, constitutive_variables_curvature_slave, PatchType::Slave);

                //Calculate moment T2
                array_1d<double, 3> T2_master = kinematic_variables_master.t / norm_2(kinematic_variables_master.t);
                array_1d<double, 3> T2_slave = kinematic_variables_slave.t / norm_2(kinematic_variables_slave.t);

                array_1d<double, 3> T1_master = kinematic_variables_master.n / norm_2(kinematic_variables_master.n);
                array_1d<double, 3> T1_slave = kinematic_variables_slave.n / norm_2(kinematic_variables_slave.n);

                bool opposite_direction_of_trims = true; 
                if (inner_prod(T2_master, T2_slave) > 0) // tangents have the same direction
                {
                    opposite_direction_of_trims = false;
                }

                double moment_T2;
                double moment_T2_master = inner_prod(moment_vector_master, T2_master);
                double moment_T2_slave = inner_prod(moment_vector_slave, T2_slave);

                if (opposite_direction_of_trims)
                {
                    moment_T2 = moment_T2_master + moment_T2_slave;
                }
                else
                {
                    moment_T2 = moment_T2_master - moment_T2_slave;
                }

                //Calculate the first variations of the 2nd Piola-Kichhoff moments at the covariant bases
                Matrix first_variations_moment_covariant_master = ZeroMatrix(3, 3 * number_of_nodes_master);
                Matrix first_variations_moment_covariant_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

                CalculateFirstVariationMomentCovariant(point_number, first_variations_moment_covariant_master, kinematic_variables_master, constitutive_variables_curvature_master, PatchType::Master);
                CalculateFirstVariationMomentCovariant(point_number, first_variations_moment_covariant_slave, kinematic_variables_slave, constitutive_variables_curvature_slave, PatchType::Slave);

                //Calculate first variation of moment vectors
                Matrix first_variations_moment_master = ZeroMatrix(3, 3 * number_of_nodes_master);
                Matrix first_variations_moment_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

                CalculateFirstVariationMoment(point_number, first_variations_moment_master, first_variations_moment_covariant_master, kinematic_variables_master, constitutive_variables_curvature_master, PatchType::Master);
                CalculateFirstVariationMoment(point_number, first_variations_moment_slave, first_variations_moment_covariant_slave, kinematic_variables_slave, constitutive_variables_curvature_slave, PatchType::Slave);

                // Matrix first_variations_moment = ZeroMatrix(3, mat_size);
                // for (SizeType i=0;i<3 * number_of_nodes_master;++i){
                //     first_variations_moment(0, i) = first_variations_moment_master(0, i);
                //     first_variations_moment(1, i) = first_variations_moment_master(1, i);
                //     first_variations_moment(2, i) = first_variations_moment_master(2, i);
                // }
                // for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
                //     first_variations_moment(0, i + 3 * number_of_nodes_master) = -first_variations_moment_slave(0, i);
                //     first_variations_moment(1, i + 3 * number_of_nodes_master) = -first_variations_moment_slave(1, i);
                //     first_variations_moment(2, i + 3 * number_of_nodes_master) = -first_variations_moment_slave(2, i);
                // }

                //Calculate first variation of moment vectors T2
                Vector first_variations_moment_master_T2 = ZeroVector(number_of_nodes_master * 3);
                Vector first_variations_moment_slave_T2 = ZeroVector(number_of_nodes_slave * 3);
                Vector first_variations_moment_T2 = ZeroVector(mat_size);

                for (SizeType i=0;i<3 * number_of_nodes_master;++i){
                    first_variations_moment_master_T2(i) =  T1_master[0]*first_variations_moment_master(0,i) + T1_master[1]*first_variations_moment_master(1,i) +
                                                        T1_master[2]*first_variations_moment_master(2,i);
                }
                for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
                    first_variations_moment_slave_T2(i) = T1_slave[0]*first_variations_moment_slave(0,i) + T1_slave[1]*first_variations_moment_slave(1,i) +
                                                        T1_slave[2]*first_variations_moment_slave(2,i);
                } 

                for (IndexType i = 0; i < 3 * number_of_nodes_master; i++)
                {
                    first_variations_moment_T2(i) = first_variations_moment_master_T2(i);
                }

                for (IndexType i = 0; i < 3 * number_of_nodes_slave; i++)
                {
                    if (opposite_direction_of_trims)
                    {
                        first_variations_moment_T2(i + 3 * number_of_nodes_master) = first_variations_moment_slave_T2(i);
                    }
                    else
                    {
                        first_variations_moment_T2(i + 3 * number_of_nodes_master) = -first_variations_moment_slave_T2(i);
                    }
                }

                //Compute the necessary products needed for the second variations of the traction vectors
                Matrix Pi_master = ZeroMatrix(3, 3);
                Matrix Pi_slave = ZeroMatrix(3, 3);

                CalculateSecondVariationMomentProduct(point_number, Pi_master, kinematic_variables_master, constitutive_variables_curvature_master, PatchType::Master);
                CalculateSecondVariationMomentProduct(point_number, Pi_slave, kinematic_variables_slave, constitutive_variables_curvature_slave, PatchType::Slave);

                // array_1d<double, 3> second_variations_moment_product_vector_master = prod(trans(Pi_master), rotation_vector_master);
                // array_1d<double, 3> second_variations_moment_product_vector_slave = prod(trans(Pi_slave), rotation_vector_slave);
                // array_1d<double, 3> second_variations_moment_product_vector_master_slave = prod(trans(Pi_slave), rotation_vector_master);
                // array_1d<double, 3> second_variations_moment_product_vector_slave_master = prod(trans(Pi_master), rotation_vector_slave);

                array_1d<double, 3> second_variations_moment_product_vector_master = prod(trans(Pi_master), T2_master);
                array_1d<double, 3> second_variations_moment_product_vector_slave = prod(trans(Pi_slave), T2_slave);
                array_1d<double, 3> second_variations_moment_product_vector_master_slave = prod(trans(Pi_slave), T2_master);
                array_1d<double, 3> second_variations_moment_product_vector_slave_master = prod(trans(Pi_master), T2_slave);

                //Calculate second variation of moment vectors
                Matrix second_variations_moment_master = ZeroMatrix(3 * number_of_nodes_master, 3 * number_of_nodes_master);
                Matrix second_variations_moment_slave = ZeroMatrix(3 * number_of_nodes_slave, 3 * number_of_nodes_slave);

                // CalculateSecondVariationMoment(point_number, second_variations_moment_master, kinematic_variables_master, first_variations_moment_covariant_master, rotation_vector_master, rotation_vector_slave, 
                //                                second_variations_moment_product_vector_master, second_variations_moment_product_vector_slave_master, PatchType::Master);
                // CalculateSecondVariationMoment(point_number, second_variations_moment_slave, kinematic_variables_slave, first_variations_moment_covariant_slave, rotation_vector_master, rotation_vector_slave, 
                //                                second_variations_moment_product_vector_slave, second_variations_moment_product_vector_master_slave, PatchType::Slave);

                CalculateSecondVariationMomentT2(point_number, second_variations_moment_master, kinematic_variables_master, first_variations_moment_covariant_master, T2_master, T2_slave, 
                                                second_variations_moment_product_vector_master, second_variations_moment_product_vector_slave_master, PatchType::Master);
                CalculateSecondVariationMomentT2(point_number, second_variations_moment_slave, kinematic_variables_slave, first_variations_moment_covariant_slave, T2_master, T2_slave, 
                                                second_variations_moment_product_vector_slave, second_variations_moment_product_vector_master_slave, PatchType::Slave);

                second_variations_moment_master *= diff_phi(0);
                second_variations_moment_slave *= diff_phi(0);

                if (!opposite_direction_of_trims)
                    second_variations_moment_slave *= -1;

                if (CalculateStiffnessMatrixFlag) {
                    for (IndexType i = 0; i < mat_size; ++i)
                    {
                        for (IndexType j = 0; j < mat_size; ++j)
                        {
                            rLeftHandSideMatrix(i, j) += (first_variations_moment_T2(i) * phi_r(j) + phi_r(i) * first_variations_moment_T2(j)) * integration_weight * determinant_jacobian * -gammaTilde;
                            rLeftHandSideMatrix(i, j) += moment_T2 * phi_rs(i, j) * integration_weight * determinant_jacobian * -gammaTilde;
                            rLeftHandSideMatrix(i, j) += (phi_r(i) * phi_r(j) + diff_phi(0) * phi_rs(i, j)) * integration_weight * determinant_jacobian * stabilization_rotation_parameter;
                        }
                    }
                    for (IndexType i = 0; i < 3 * number_of_nodes_master; i++)
                    {
                        for (IndexType j = 0; j < 3 * number_of_nodes_master; j++)
                        {
                            rLeftHandSideMatrix(i, j) += second_variations_moment_master(i, j) * integration_weight * determinant_jacobian * -gammaTilde;
                        }
                    }
                    for (IndexType i = 0; i < 3 * number_of_nodes_slave; i++)
                    {
                        for (IndexType j = 0; j < 3 * number_of_nodes_slave; j++)
                        {
                            rLeftHandSideMatrix(i + 3 * number_of_nodes_master, j + 3 * number_of_nodes_master) += second_variations_moment_slave(i, j) * integration_weight * determinant_jacobian * -gammaTilde;
                        }
                    }
                } 

                if (CalculateResidualVectorFlag) {
                    for (IndexType i = 0; i < mat_size; ++i)
                    {
                        rRightHandSideVector[i] -= (diff_phi(0) * first_variations_moment_T2(i)) * integration_weight * determinant_jacobian * -gammaTilde;
                        rRightHandSideVector[i] -= (moment_T2 * phi_r(i)) * integration_weight * determinant_jacobian * -gammaTilde;
                        rRightHandSideVector[i] -= (diff_phi(0) * phi_r(i)) * integration_weight * determinant_jacobian * stabilization_rotation_parameter;
                    }
                }
            }

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
        if (Is(IgaFlags::FIX_ROTATION_X))
        {
            if (m_B_ab_covariant_vector_master.size() != r_number_of_integration_points_master)
                m_B_ab_covariant_vector_master.resize(r_number_of_integration_points_master);
            if (m_B_ab_covariant_vector_slave.size() != r_number_of_integration_points_slave)
                m_B_ab_covariant_vector_slave.resize(r_number_of_integration_points_slave);
        }

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

            if (Is(IgaFlags::FIX_ROTATION_X))
            {
                m_B_ab_covariant_vector_master[point_number] = kinematic_variables_reference_master.b_ab_covariant;
                m_B_ab_covariant_vector_slave[point_number] = kinematic_variables_reference_slave.b_ab_covariant;
            }

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
            ConstitutiveVariables constitutive_variables_curvature_master(3);
            ConstitutiveVariables constitutive_variables_curvature_slave(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_master,
                constitutive_variables_membrane_master,
                constitutive_variables_curvature_master,
                constitutive_law_parameters_master,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Master);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_slave,
                constitutive_variables_membrane_slave,
                constitutive_variables_curvature_slave,
                constitutive_law_parameters_slave,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Slave);

            //Prestress component
            array_1d<double, 3> prestress_master = ZeroVector(3);
            array_1d<double, 3> prestress_slave = ZeroVector(3);

            if(this->GetProperties().Has(PRESTRESS))
            {
                prestress_master = GetProperties().GetSubProperties().front()[PRESTRESS]*GetProperties().GetSubProperties().front()[THICKNESS];
                prestress_slave = GetProperties().GetSubProperties().back()[PRESTRESS]*GetProperties().GetSubProperties().back()[THICKNESS];
            }
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
            double determinant_jacobian = determinant_jacobian_vector_initial[point_number];
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
        }
        KRATOS_CATCH("")
    }

    void CouplingNitscheCondition::CalculateNitscheStabilizationRotationMatrix(
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
        if (Is(IgaFlags::FIX_ROTATION_X))
        {
            if (m_B_ab_covariant_vector_master.size() != r_number_of_integration_points_master)
                m_B_ab_covariant_vector_master.resize(r_number_of_integration_points_master);
            if (m_B_ab_covariant_vector_slave.size() != r_number_of_integration_points_slave)
                m_B_ab_covariant_vector_slave.resize(r_number_of_integration_points_slave);
        }

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

            if (Is(IgaFlags::FIX_ROTATION_X))
            {
                m_B_ab_covariant_vector_master[point_number] = kinematic_variables_reference_master.b_ab_covariant;
                m_B_ab_covariant_vector_slave[point_number] = kinematic_variables_reference_slave.b_ab_covariant;
            }

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
            ConstitutiveVariables constitutive_variables_curvature_master(3);
            ConstitutiveVariables constitutive_variables_curvature_slave(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_master,
                constitutive_variables_membrane_master,
                constitutive_variables_curvature_master,
                constitutive_law_parameters_master,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Master);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_slave,
                constitutive_variables_membrane_slave,
                constitutive_variables_curvature_slave,
                constitutive_law_parameters_slave,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Slave);

            //Prestress component
            array_1d<double, 3> prestress_master = ZeroVector(3);
            array_1d<double, 3> prestress_slave = ZeroVector(3);
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

            //Calculate moment vectors
            array_1d<double, 3> moment_vector_master;
            array_1d<double, 3> moment_vector_slave;

            CalculateMoment(point_number, moment_vector_master, kinematic_variables_master, constitutive_variables_curvature_master, PatchType::Master);
            CalculateMoment(point_number, moment_vector_slave, kinematic_variables_slave, constitutive_variables_curvature_slave, PatchType::Slave);

            //Calculate moment T2
            array_1d<double, 3> T2_master = kinematic_variables_master.t / norm_2(kinematic_variables_master.t);
            array_1d<double, 3> T2_slave = kinematic_variables_slave.t / norm_2(kinematic_variables_slave.t);

            array_1d<double, 3> T1_master = kinematic_variables_master.n / norm_2(kinematic_variables_master.n);
            array_1d<double, 3> T1_slave = kinematic_variables_slave.n / norm_2(kinematic_variables_slave.n);

            bool opposite_direction_of_trims = true; 
            if (inner_prod(T2_master, T2_slave) > 0) // tangents have the same direction
            {
                opposite_direction_of_trims = false;
            }

            double moment_T2;
            double moment_T2_master = inner_prod(moment_vector_master, T2_master);
            double moment_T2_slave = inner_prod(moment_vector_slave, T2_slave);

            if (opposite_direction_of_trims)
            {
                moment_T2 = moment_T2_master + moment_T2_slave;
            }
            else
            {
                moment_T2 = moment_T2_master - moment_T2_slave;
            }

            //Calculate the first variations of the 2nd Piola-Kichhoff moments at the covariant bases
            Matrix first_variations_moment_covariant_master = ZeroMatrix(3, 3 * number_of_nodes_master);
            Matrix first_variations_moment_covariant_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

            CalculateFirstVariationMomentCovariant(point_number, first_variations_moment_covariant_master, kinematic_variables_master, constitutive_variables_curvature_master, PatchType::Master);
            CalculateFirstVariationMomentCovariant(point_number, first_variations_moment_covariant_slave, kinematic_variables_slave, constitutive_variables_curvature_slave, PatchType::Slave);

            //Calculate first variation of moment vectors
            Matrix first_variations_moment_master = ZeroMatrix(3, 3 * number_of_nodes_master);
            Matrix first_variations_moment_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

            CalculateFirstVariationMoment(point_number, first_variations_moment_master, first_variations_moment_covariant_master, kinematic_variables_master, constitutive_variables_curvature_master, PatchType::Master);
            CalculateFirstVariationMoment(point_number, first_variations_moment_slave, first_variations_moment_covariant_slave, kinematic_variables_slave, constitutive_variables_curvature_slave, PatchType::Slave);

            // Matrix first_variations_moment = ZeroMatrix(3, mat_size);
            // for (SizeType i=0;i<3 * number_of_nodes_master;++i){
            //     first_variations_moment(0, i) = first_variations_moment_master(0, i);
            //     first_variations_moment(1, i) = first_variations_moment_master(1, i);
            //     first_variations_moment(2, i) = first_variations_moment_master(2, i);
            // }
            // for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
            //     first_variations_moment(0, i + 3 * number_of_nodes_master) = -first_variations_moment_slave(0, i);
            //     first_variations_moment(1, i + 3 * number_of_nodes_master) = -first_variations_moment_slave(1, i);
            //     first_variations_moment(2, i + 3 * number_of_nodes_master) = -first_variations_moment_slave(2, i);
            // }

            //Calculate first variation of moment vectors T2
            Vector first_variations_moment_master_T2 = ZeroVector(number_of_nodes_master * 3);
            Vector first_variations_moment_slave_T2 = ZeroVector(number_of_nodes_slave * 3);
            Vector first_variations_moment_T2 = ZeroVector(mat_size);

            for (SizeType i=0;i<3 * number_of_nodes_master;++i){
                first_variations_moment_master_T2(i) = T1_master[0]*first_variations_moment_master(0,i) + T1_master[1]*first_variations_moment_master(1,i) +
                                                       T1_master[2]*first_variations_moment_master(2,i);
            }
            for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
                first_variations_moment_slave_T2(i) = T1_slave[0]*first_variations_moment_slave(0,i) + T1_slave[1]*first_variations_moment_slave(1,i) +
                                                      T1_slave[2]*first_variations_moment_slave(2,i);
            } 

            for (IndexType i = 0; i < 3 * number_of_nodes_master; i++)
            {
                first_variations_moment_T2(i) = first_variations_moment_master_T2(i);
            }

            for (IndexType i = 0; i < 3 * number_of_nodes_slave; i++)
            {
                if (opposite_direction_of_trims)
                {
                    first_variations_moment_T2(i + 3 * number_of_nodes_master) = first_variations_moment_slave_T2(i);
                }
                else
                {
                    first_variations_moment_T2(i + 3 * number_of_nodes_master) = -first_variations_moment_slave_T2(i);
                }
            }

            //Compute the necessary products needed for the second variations of the traction vectors
            Matrix Pi_master = ZeroMatrix(3, 3);
            Matrix Pi_slave = ZeroMatrix(3, 3);

            CalculateSecondVariationMomentProduct(point_number, Pi_master, kinematic_variables_master, constitutive_variables_curvature_master, PatchType::Master);
            CalculateSecondVariationMomentProduct(point_number, Pi_slave, kinematic_variables_slave, constitutive_variables_curvature_slave, PatchType::Slave);

            // array_1d<double, 3> second_variations_moment_product_vector_master = prod(trans(Pi_master), rotation_vector_master);
            // array_1d<double, 3> second_variations_moment_product_vector_slave = prod(trans(Pi_slave), rotation_vector_slave);
            // array_1d<double, 3> second_variations_moment_product_vector_master_slave = prod(trans(Pi_slave), rotation_vector_master);
            // array_1d<double, 3> second_variations_moment_product_vector_slave_master = prod(trans(Pi_master), rotation_vector_slave);

            array_1d<double, 3> second_variations_moment_product_vector_master = prod(trans(Pi_master), T2_master);
            array_1d<double, 3> second_variations_moment_product_vector_slave = prod(trans(Pi_slave), T2_slave);
            array_1d<double, 3> second_variations_moment_product_vector_master_slave = prod(trans(Pi_slave), T2_master);
            array_1d<double, 3> second_variations_moment_product_vector_slave_master = prod(trans(Pi_master), T2_slave);

            //Calculate second variation of moment vectors
            Matrix second_variations_moment_master = ZeroMatrix(3 * number_of_nodes_master, 3 * number_of_nodes_master);
            Matrix second_variations_moment_slave = ZeroMatrix(3 * number_of_nodes_slave, 3 * number_of_nodes_slave);

            // CalculateSecondVariationMoment(point_number, second_variations_moment_master, kinematic_variables_master, first_variations_moment_covariant_master, rotation_vector_master, rotation_vector_slave, 
            //                                second_variations_moment_product_vector_master, second_variations_moment_product_vector_slave_master, PatchType::Master);
            // CalculateSecondVariationMoment(point_number, second_variations_moment_slave, kinematic_variables_slave, first_variations_moment_covariant_slave, rotation_vector_master, rotation_vector_slave, 
            //                                second_variations_moment_product_vector_slave, second_variations_moment_product_vector_master_slave, PatchType::Slave);

            if(norm_2(kinematic_variables_reference_master.a3_tilde) > tol_surface_normal && norm_2(kinematic_variables_reference_slave.a3_tilde) > tol_surface_normal)
            {
                CalculateSecondVariationMomentT2(point_number, second_variations_moment_master, kinematic_variables_master, first_variations_moment_covariant_master, T2_master, T2_slave, 
                                                    second_variations_moment_product_vector_master, second_variations_moment_product_vector_slave_master, PatchType::Master);
                CalculateSecondVariationMomentT2(point_number, second_variations_moment_slave, kinematic_variables_slave, first_variations_moment_covariant_slave, T2_master, T2_slave, 
                                                    second_variations_moment_product_vector_slave, second_variations_moment_product_vector_master_slave, PatchType::Slave);
                
                second_variations_moment_master *= moment_T2;
                second_variations_moment_slave *= moment_T2;

                if (!opposite_direction_of_trims)
                    second_variations_moment_slave *= -1;
            }
            
            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinant_jacobian = determinant_jacobian_vector_initial[point_number];
            const double gammaTilde = 0.5;

            // Assembly
            for (IndexType i = 0; i < mat_size; ++i)
            {
                for (IndexType j = 0; j < mat_size; ++j)
                {
                    rLeftHandSideMatrix(i, j) += (first_variations_moment_T2(i) * first_variations_moment_T2(j) * 2) * integration_weight * determinant_jacobian  * gammaTilde * gammaTilde;
                }
            }

            for (IndexType i = 0; i < 3 * number_of_nodes_master; i++)
            {
                for (IndexType j = 0; j < 3 * number_of_nodes_master; j++)
                {
                    rLeftHandSideMatrix(i, j) += second_variations_moment_master(i, j) * integration_weight * determinant_jacobian * gammaTilde * gammaTilde;
                }
            }

            for (IndexType i = 0; i < 3 * number_of_nodes_slave; i++)
            {
                for (IndexType j = 0; j < 3 * number_of_nodes_slave; j++)
                {
                    rLeftHandSideMatrix(i + 3 * number_of_nodes_master, j + 3 * number_of_nodes_master) += second_variations_moment_slave(i, j) * integration_weight * determinant_jacobian * gammaTilde * gammaTilde;
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
    
    ///@}
    ///@name Results on Gauss Points
    ///@{

    void CouplingNitscheCondition::CalculateOnIntegrationPoints(
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

        if(rVariable==NITSCHE_STABILIZATION_FACTOR)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            {
                rOutput[point_number] = GetProperties()[NITSCHE_STABILIZATION_FACTOR];
            }
        }
    }

    void CouplingNitscheCondition::CalculateMoment(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rMoment,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature, 
        const PatchType& rPatch)
    {
        // Transform the 2nd Piola-kirchhoff moment in the covariant systems
        array_1d<double, 3> moment_vector_covariant;
        array_1d<double, 2> n_contravariant_vector;

        if (rPatch==PatchType::Master)
        {
            moment_vector_covariant = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesCurvature.StressVector);
            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
            moment_vector_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesCurvature.StressVector);
            n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        }

        // Compute the stress components
        Matrix Palphabeta = ZeroMatrix(2, 2);
        Palphabeta(0,0) = moment_vector_covariant[0];
        Palphabeta(1,1) = moment_vector_covariant[1];
        Palphabeta(0,1) = moment_vector_covariant[2];
        Palphabeta(1,0) = Palphabeta(0,1);
        
        // Compute the traction vectors
        rMoment[0] = rActualKinematic.a1[0]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                   + rActualKinematic.a2[0]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
        rMoment[1] = rActualKinematic.a1[1]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                   + rActualKinematic.a2[1]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
        rMoment[2] = rActualKinematic.a1[2]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                   + rActualKinematic.a2[2]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
    }

    void CouplingNitscheCondition::CalculateFirstVariationMomentCovariant(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationMomentCovariant,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature, 
        const PatchType& rPatch)
    {
        IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod());

        Matrix da3 = ZeroMatrix(3, 3);
        Matrix dn = ZeroMatrix(3, 3);
        Matrix b = ZeroMatrix(3, mat_size);

        double inv_dA = 1 / rActualKinematic.dA;
        double inv_dA3 = 1 / std::pow(rActualKinematic.dA, 3);

        Matrix H = ZeroMatrix(3, 3);
        CalculateHessian(H, r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex), rPatch);

        //Compute the first variation of the Green-Lagrange curvature
        Matrix dK_curvilinear = ZeroMatrix(3, mat_size);
        Matrix dK_cartesian = ZeroMatrix(3, mat_size);
        Matrix T_patch = ZeroMatrix(3, 3);

        if (rPatch==PatchType::Master)
        {
            T_patch = m_T_vector_master[IntegrationPointIndex];
        }
        else
        {
            T_patch = m_T_vector_slave[IntegrationPointIndex];
        }

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
            dK_curvilinear(0, index) = 0 - (r_DDN_DDe(i, 0) * rActualKinematic.a3[0] + H(0, 0) * dn(0, 0) + H(1, 0) * dn(0, 1) + H(2, 0) * dn(0, 2));
            dK_curvilinear(0, index + 1) = 0 - (r_DDN_DDe(i, 0) * rActualKinematic.a3[1] + H(0, 0) * dn(1, 0) + H(1, 0) * dn(1, 1) + H(2, 0) * dn(1, 2));
            dK_curvilinear(0, index + 2) = 0 - (r_DDN_DDe(i, 0) * rActualKinematic.a3[2] + H(0, 0) * dn(2, 0) + H(1, 0) * dn(2, 1) + H(2, 0) * dn(2, 2));

            //second line
            dK_curvilinear(1, index) = 0 - (r_DDN_DDe(i, 2) * rActualKinematic.a3[0] + H(0, 1) * dn(0, 0) + H(1, 1) * dn(0, 1) + H(2, 1) * dn(0, 2));
            dK_curvilinear(1, index + 1) = 0 - (r_DDN_DDe(i, 2) * rActualKinematic.a3[1] + H(0, 1) * dn(1, 0) + H(1, 1) * dn(1, 1) + H(2, 1) * dn(1, 2));
            dK_curvilinear(1, index + 2) = 0 - (r_DDN_DDe(i, 2) * rActualKinematic.a3[2] + H(0, 1) * dn(2, 0) + H(1, 1) * dn(2, 1) + H(2, 1) * dn(2, 2));

            //third line
            dK_curvilinear(2, index) = 0 - (r_DDN_DDe(i, 1) * rActualKinematic.a3[0] + H(0, 2) * dn(0, 0) + H(1, 2) * dn(0, 1) + H(2, 2) * dn(0, 2));
            dK_curvilinear(2, index + 1) = 0 - (r_DDN_DDe(i, 1) * rActualKinematic.a3[1] + H(0, 2) * dn(1, 0) + H(1, 2) * dn(1, 1) + H(2, 2) * dn(1, 2));
            dK_curvilinear(2, index + 2) = 0 - (r_DDN_DDe(i, 1) * rActualKinematic.a3[2] + H(0, 2) * dn(2, 0) + H(1, 2) * dn(2, 1) + H(2, 2) * dn(2, 2));
        }

        noalias(dK_cartesian) = -prod(T_patch, dK_curvilinear);

        //Compute the first variations of the 2nd Piola-Kichhoff moments in the local Cartesian bases
        Matrix first_variations_moment_cartesian = ZeroMatrix(3, mat_size);
        first_variations_moment_cartesian = prod(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix,dK_cartesian);

        //Transform the first variations of the 2nd Piola-Kichhoff moments at the covariant bases
        if (rPatch==PatchType::Master)
        {
            rFirstVariationMomentCovariant = prod(m_T_hat_vector_master[IntegrationPointIndex], first_variations_moment_cartesian);
        }
        else
        {
            rFirstVariationMomentCovariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], first_variations_moment_cartesian);
        }
    }

    void CouplingNitscheCondition::CalculateFirstVariationMoment(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationMoment,
        Matrix& rFirstVariationMomentCovariant,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature, 
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

        rFirstVariationMoment = prod(n_a, rFirstVariationMomentCovariant);

        //2. derivative normal vector * moment covariant
        
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> moment_vector_covariant;
        if (rPatch==PatchType::Master)
        {
            moment_vector_covariant = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesCurvature.StressVector);
        }
        else
        {
            moment_vector_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesCurvature.StressVector);
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

        rFirstVariationMoment += r_DN_Dxi*(n_contravariant_vector[0]*moment_vector_covariant[0] + n_contravariant_vector[1]*moment_vector_covariant[2])+
                                 r_DN_Deta*(n_contravariant_vector[1]*moment_vector_covariant[1] + n_contravariant_vector[0]*moment_vector_covariant[2]);
    }

    void CouplingNitscheCondition::CalculateSecondVariationMomentProduct(
        IndexType IntegrationPointIndex,
        Matrix& rPi,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature, 
        const PatchType& rPatch)
    {
        array_1d<double, 2> n_contravariant_vector;

        if (rPatch==PatchType::Master)
        {
            rPi = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesCurvature.ConstitutiveMatrix);
            rPi = prod(rPi, m_T_vector_master[IntegrationPointIndex]);

            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
            rPi = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesCurvature.ConstitutiveMatrix);
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

    void CouplingNitscheCondition::CalculateSecondVariationMoment(
        IndexType IntegrationPointIndex,
        Matrix& rSecondVariationMoment,
        const KinematicVariables& rActualKinematic,
        Matrix& rFirstVariationMomentCovariant, 
        array_1d<double, 3>& rRotationMaster,
        array_1d<double, 3>& rRotationSlave,
        array_1d<double, 3>& rSecondVariationMomentProduct,
        array_1d<double, 3>& rSecondVariationMomentProductMasterSlave,
        const PatchType& rPatch)
    {
        IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);

        const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod());
 
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

        // 1. second variation of moment * base vector
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
        CalculateHessian(H, r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod()), rPatch);

        // first variation of curvature w.r.t. dof
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

        // second variation of curvature w.r.t. dofs
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

                rSecondVariationMoment(r, s) += ddK_cu[0]*rSecondVariationMomentProduct(0);
                rSecondVariationMoment(r, s) += ddK_cu[1]*rSecondVariationMomentProduct(1);
                rSecondVariationMoment(r, s) += ddK_cu[2]*rSecondVariationMomentProduct(2);

                rSecondVariationMoment(r, s) -= ddK_cu[0]*rSecondVariationMomentProductMasterSlave(0);
                rSecondVariationMoment(r, s) -= ddK_cu[1]*rSecondVariationMomentProductMasterSlave(1);
                rSecondVariationMoment(r, s) -= ddK_cu[2]*rSecondVariationMomentProductMasterSlave(2);

                rSecondVariationMoment(s, r) = rSecondVariationMoment(r, s);
            }
        }
        
        // 2. first variation of moment * first variation of base vector
        array_1d<double, 2> n_contravariant_vector; 

        if (rPatch==PatchType::Master)
        {
            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
            n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        }

        Vector first_variations_moment_covariant_1 = ZeroVector(mat_size);
        Vector first_variations_moment_covariant_2 = ZeroVector(mat_size);
        Vector first_variations_moment_covariant_3 = ZeroVector(mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            first_variations_moment_covariant_1(r) = rFirstVariationMomentCovariant(0, r);
            first_variations_moment_covariant_2(r) = rFirstVariationMomentCovariant(1, r);
            first_variations_moment_covariant_3(r) = rFirstVariationMomentCovariant(2, r);
        }

        // rotation * first variation stress covariant
        Matrix rotation_dNCovariant1_master = ZeroMatrix(3, mat_size); 
        Matrix rotation_dNCovariant2_master = ZeroMatrix(3, mat_size); 
        Matrix rotation_dNCovariant3_master = ZeroMatrix(3, mat_size); 

        Matrix rotation_dNCovariant1_slave = ZeroMatrix(3, mat_size); 
        Matrix rotation_dNCovariant2_slave = ZeroMatrix(3, mat_size); 
        Matrix rotation_dNCovariant3_slave = ZeroMatrix(3, mat_size); 

        for (IndexType r = 0; r < mat_size; r++)
        {
            rotation_dNCovariant1_master(0, r) = rRotationMaster(0)*first_variations_moment_covariant_1(r);
            rotation_dNCovariant1_master(1, r) = rRotationMaster(1)*first_variations_moment_covariant_1(r);
            rotation_dNCovariant1_master(2, r) = rRotationMaster(2)*first_variations_moment_covariant_1(r);

            rotation_dNCovariant2_master(0, r) = rRotationMaster(0)*first_variations_moment_covariant_2(r);
            rotation_dNCovariant2_master(1, r) = rRotationMaster(1)*first_variations_moment_covariant_2(r);
            rotation_dNCovariant2_master(2, r) = rRotationMaster(2)*first_variations_moment_covariant_2(r);

            rotation_dNCovariant3_master(0, r) = rRotationMaster(0)*first_variations_moment_covariant_3(r);
            rotation_dNCovariant3_master(1, r) = rRotationMaster(1)*first_variations_moment_covariant_3(r);
            rotation_dNCovariant3_master(2, r) = rRotationMaster(2)*first_variations_moment_covariant_3(r);

            rotation_dNCovariant1_slave(0, r) = rRotationSlave(0)*first_variations_moment_covariant_1(r);
            rotation_dNCovariant1_slave(1, r) = rRotationSlave(1)*first_variations_moment_covariant_1(r);
            rotation_dNCovariant1_slave(2, r) = rRotationSlave(2)*first_variations_moment_covariant_1(r);

            rotation_dNCovariant2_slave(0, r) = rRotationSlave(0)*first_variations_moment_covariant_2(r);
            rotation_dNCovariant2_slave(1, r) = rRotationSlave(1)*first_variations_moment_covariant_2(r);
            rotation_dNCovariant2_slave(2, r) = rRotationSlave(2)*first_variations_moment_covariant_2(r);

            rotation_dNCovariant3_slave(0, r) = rRotationSlave(0)*first_variations_moment_covariant_3(r);
            rotation_dNCovariant3_slave(1, r) = rRotationSlave(1)*first_variations_moment_covariant_3(r);
            rotation_dNCovariant3_slave(2, r) = rRotationSlave(2)*first_variations_moment_covariant_3(r);
        }

        if (rPatch==PatchType::Slave)
        {
            rotation_dNCovariant1_master *= -1.0;
            rotation_dNCovariant2_master *= -1.0;
            rotation_dNCovariant3_master *= -1.0;

            rotation_dNCovariant1_slave *= -1.0;
            rotation_dNCovariant2_slave *= -1.0;
            rotation_dNCovariant3_slave *= -1.0;
        }

        rSecondVariationMoment += prod(trans(r_DN_Dxi), rotation_dNCovariant1_master)*n_contravariant_vector(0);
        rSecondVariationMoment += prod(trans(r_DN_Deta), rotation_dNCovariant2_master)*n_contravariant_vector(1);
        rSecondVariationMoment += (prod(trans(r_DN_Dxi), rotation_dNCovariant3_master)*n_contravariant_vector(1) + prod(trans(r_DN_Deta), rotation_dNCovariant3_master)*n_contravariant_vector(0));

        rSecondVariationMoment += prod(trans(rotation_dNCovariant1_master), r_DN_Dxi)*n_contravariant_vector(0);
        rSecondVariationMoment += prod(trans(rotation_dNCovariant2_master), r_DN_Deta)*n_contravariant_vector(1);
        rSecondVariationMoment += (prod(trans(rotation_dNCovariant3_master), r_DN_Dxi)*n_contravariant_vector(1) + prod(trans(rotation_dNCovariant3_master), r_DN_Deta)*n_contravariant_vector(0));

        rSecondVariationMoment -= prod(trans(r_DN_Dxi), rotation_dNCovariant1_slave)*n_contravariant_vector(0);
        rSecondVariationMoment -= prod(trans(r_DN_Deta), rotation_dNCovariant2_slave)*n_contravariant_vector(1);
        rSecondVariationMoment -= (prod(trans(r_DN_Dxi), rotation_dNCovariant3_slave)*n_contravariant_vector(1) + prod(trans(r_DN_Deta), rotation_dNCovariant3_slave)*n_contravariant_vector(0));

        rSecondVariationMoment -= prod(trans(rotation_dNCovariant1_slave), r_DN_Dxi)*n_contravariant_vector(0);
        rSecondVariationMoment -= prod(trans(rotation_dNCovariant2_slave), r_DN_Deta)*n_contravariant_vector(1);
        rSecondVariationMoment -= (prod(trans(rotation_dNCovariant3_slave), r_DN_Dxi)*n_contravariant_vector(1) + prod(trans(rotation_dNCovariant3_slave), r_DN_Deta)*n_contravariant_vector(0));
    }

    void CouplingNitscheCondition::CalculateSecondVariationMomentT2(
        IndexType IntegrationPointIndex,
        Matrix& rSecondVariationMoment,
        const KinematicVariables& rActualKinematic,
        Matrix& rFirstVariationMomentCovariant, 
        array_1d<double, 3>& T2Master,
        array_1d<double, 3>& T2Slave,
        array_1d<double, 3>& rSecondVariationMomentProduct,
        array_1d<double, 3>& rSecondVariationMomentProductMasterSlave,
        const PatchType& rPatch)
    {
        IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);

        const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod());
 
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

        //Second variation of moment * base vector
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
        CalculateHessian(H, r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod()), rPatch);

        //First variation of curvature w.r.t. dof
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

        //Second variation of curvature w.r.t. dofs
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

                rSecondVariationMoment(r, s) += ddK_cu[0]*rSecondVariationMomentProduct(0);
                rSecondVariationMoment(r, s) += ddK_cu[1]*rSecondVariationMomentProduct(1);
                rSecondVariationMoment(r, s) += ddK_cu[2]*rSecondVariationMomentProduct(2);

                rSecondVariationMoment(r, s) -= ddK_cu[0]*rSecondVariationMomentProductMasterSlave(0);
                rSecondVariationMoment(r, s) -= ddK_cu[1]*rSecondVariationMomentProductMasterSlave(1);
                rSecondVariationMoment(r, s) -= ddK_cu[2]*rSecondVariationMomentProductMasterSlave(2);

                rSecondVariationMoment(s, r) = rSecondVariationMoment(r, s);
            }
        }
        
        //First variation of moment * first variation of base vector
        array_1d<double, 2> n_contravariant_vector; 

        if (rPatch==PatchType::Master)
        {
            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
            n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        }

        Vector first_variations_moment_covariant_1 = ZeroVector(mat_size);
        Vector first_variations_moment_covariant_2 = ZeroVector(mat_size);
        Vector first_variations_moment_covariant_3 = ZeroVector(mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            first_variations_moment_covariant_1(r) = rFirstVariationMomentCovariant(0, r);
            first_variations_moment_covariant_2(r) = rFirstVariationMomentCovariant(1, r);
            first_variations_moment_covariant_3(r) = rFirstVariationMomentCovariant(2, r);
        }

        //Rotation * first variation stress covariant
        Matrix rotation_dNCovariant1_master = ZeroMatrix(3, mat_size); 
        Matrix rotation_dNCovariant2_master = ZeroMatrix(3, mat_size); 
        Matrix rotation_dNCovariant3_master = ZeroMatrix(3, mat_size); 

        Matrix rotation_dNCovariant1_slave = ZeroMatrix(3, mat_size); 
        Matrix rotation_dNCovariant2_slave = ZeroMatrix(3, mat_size); 
        Matrix rotation_dNCovariant3_slave = ZeroMatrix(3, mat_size); 

        for (IndexType r = 0; r < mat_size; r++)
        {
            rotation_dNCovariant1_master(0, r) = T2Master(0)*first_variations_moment_covariant_1(r);
            rotation_dNCovariant1_master(1, r) = T2Master(1)*first_variations_moment_covariant_1(r);
            rotation_dNCovariant1_master(2, r) = T2Master(2)*first_variations_moment_covariant_1(r);

            rotation_dNCovariant2_master(0, r) = T2Master(0)*first_variations_moment_covariant_2(r);
            rotation_dNCovariant2_master(1, r) = T2Master(1)*first_variations_moment_covariant_2(r);
            rotation_dNCovariant2_master(2, r) = T2Master(2)*first_variations_moment_covariant_2(r);

            rotation_dNCovariant3_master(0, r) = T2Master(0)*first_variations_moment_covariant_3(r);
            rotation_dNCovariant3_master(1, r) = T2Master(1)*first_variations_moment_covariant_3(r);
            rotation_dNCovariant3_master(2, r) = T2Master(2)*first_variations_moment_covariant_3(r);

            rotation_dNCovariant1_slave(0, r) = T2Slave(0)*first_variations_moment_covariant_1(r);
            rotation_dNCovariant1_slave(1, r) = T2Slave(1)*first_variations_moment_covariant_1(r);
            rotation_dNCovariant1_slave(2, r) = T2Slave(2)*first_variations_moment_covariant_1(r);

            rotation_dNCovariant2_slave(0, r) = T2Slave(0)*first_variations_moment_covariant_2(r);
            rotation_dNCovariant2_slave(1, r) = T2Slave(1)*first_variations_moment_covariant_2(r);
            rotation_dNCovariant2_slave(2, r) = T2Slave(2)*first_variations_moment_covariant_2(r);

            rotation_dNCovariant3_slave(0, r) = T2Slave(0)*first_variations_moment_covariant_3(r);
            rotation_dNCovariant3_slave(1, r) = T2Slave(1)*first_variations_moment_covariant_3(r);
            rotation_dNCovariant3_slave(2, r) = T2Slave(2)*first_variations_moment_covariant_3(r);
        }

        if (rPatch==PatchType::Slave)
        {
            rotation_dNCovariant1_master *= -1.0;
            rotation_dNCovariant2_master *= -1.0;
            rotation_dNCovariant3_master *= -1.0;

            rotation_dNCovariant1_slave *= -1.0;
            rotation_dNCovariant2_slave *= -1.0;
            rotation_dNCovariant3_slave *= -1.0;
        }

        rSecondVariationMoment += prod(trans(r_DN_Dxi), rotation_dNCovariant1_master)*n_contravariant_vector(0);
        rSecondVariationMoment += prod(trans(r_DN_Deta), rotation_dNCovariant2_master)*n_contravariant_vector(1);
        rSecondVariationMoment += (prod(trans(r_DN_Dxi), rotation_dNCovariant3_master)*n_contravariant_vector(1) + prod(trans(r_DN_Deta), rotation_dNCovariant3_master)*n_contravariant_vector(0));

        rSecondVariationMoment += prod(trans(rotation_dNCovariant1_master), r_DN_Dxi)*n_contravariant_vector(0);
        rSecondVariationMoment += prod(trans(rotation_dNCovariant2_master), r_DN_Deta)*n_contravariant_vector(1);
        rSecondVariationMoment += (prod(trans(rotation_dNCovariant3_master), r_DN_Dxi)*n_contravariant_vector(1) + prod(trans(rotation_dNCovariant3_master), r_DN_Deta)*n_contravariant_vector(0));

        rSecondVariationMoment -= prod(trans(r_DN_Dxi), rotation_dNCovariant1_slave)*n_contravariant_vector(0);
        rSecondVariationMoment -= prod(trans(r_DN_Deta), rotation_dNCovariant2_slave)*n_contravariant_vector(1);
        rSecondVariationMoment -= (prod(trans(r_DN_Dxi), rotation_dNCovariant3_slave)*n_contravariant_vector(1) + prod(trans(r_DN_Deta), rotation_dNCovariant3_slave)*n_contravariant_vector(0));

        rSecondVariationMoment -= prod(trans(rotation_dNCovariant1_slave), r_DN_Dxi)*n_contravariant_vector(0);
        rSecondVariationMoment -= prod(trans(rotation_dNCovariant2_slave), r_DN_Deta)*n_contravariant_vector(1);
        rSecondVariationMoment -= (prod(trans(rotation_dNCovariant3_slave), r_DN_Dxi)*n_contravariant_vector(1) + prod(trans(rotation_dNCovariant3_slave), r_DN_Deta)*n_contravariant_vector(0));
    }

    void CouplingNitscheCondition::CalculateRotationalShapeFunctions(
        IndexType IntegrationPointIndex,
        Vector &phi_r, 
        Matrix &phi_rs, 
        array_1d<double, 2> &diff_phi)
    {
        // compute rotation (master)
        array_1d<double, 3> local_tangent_master;
        GetGeometry().GetGeometryPart(0).Calculate(LOCAL_TANGENT, local_tangent_master);

        const IntegrationMethod integration_method_master = GetGeometry().GetGeometryPart(0).GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_master = GetGeometry().GetGeometryPart(0).ShapeFunctionsLocalGradients(integration_method_master);
        const Matrix& shape_functions_gradients_master = r_shape_functions_gradients_master(IntegrationPointIndex);

        const SizeType number_of_nodes_master = GetGeometry().GetGeometryPart(0).size();

        Vector phi_r_master = ZeroVector(number_of_nodes_master * 3);
        Matrix phi_rs_master = ZeroMatrix(number_of_nodes_master * 3, number_of_nodes_master * 3);
        array_1d<double, 2> phi_master;
        array_1d<double, 3> trim_tangents_master;

        CalculateRotation(IntegrationPointIndex, shape_functions_gradients_master, phi_r_master, phi_rs_master, phi_master, trim_tangents_master, local_tangent_master, true);

        // compute rotation (slave)
        array_1d<double, 3> local_tangent_slave;
        GetGeometry().GetGeometryPart(1).Calculate(LOCAL_TANGENT, local_tangent_slave);

        const IntegrationMethod integration_method_slave = GetGeometry().GetGeometryPart(1).GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_slave = GetGeometry().GetGeometryPart(1).ShapeFunctionsLocalGradients(integration_method_slave);
        const Matrix& shape_functions_gradients_slave = r_shape_functions_gradients_slave(IntegrationPointIndex);

        const SizeType number_of_nodes_slave = GetGeometry().GetGeometryPart(1).size();

        Vector phi_r_slave = ZeroVector(number_of_nodes_slave * 3);
        Matrix phi_rs_slave = ZeroMatrix(number_of_nodes_slave * 3, number_of_nodes_slave * 3);
        array_1d<double, 2> phi_slave;
        array_1d<double, 3> trim_tangents_slave;

        CalculateRotation(IntegrationPointIndex, shape_functions_gradients_slave, phi_r_slave, phi_rs_slave, phi_slave, trim_tangents_slave, local_tangent_slave, false);

        // compute phi_r, phi_rs and diff_phi
        bool opposite_direction_of_trims = true;
        if (inner_prod(trim_tangents_master, trim_tangents_slave) > 0) // tangents have the same direction (assumption for coordinates system changes)
        {
            opposite_direction_of_trims = false;
        }

        if (opposite_direction_of_trims)
        {
            diff_phi = phi_slave + phi_master;
        }
        else
        {
            diff_phi = -(phi_slave - phi_master);
        }
        
        for (IndexType i = 0; i < phi_r_master.size(); i++)
        {
            phi_r(i) = phi_r_master(i);
        }

        SizeType index = phi_r_master.size();
        for (IndexType i = 0; i < phi_r_slave.size(); i++)
        {
            if (opposite_direction_of_trims)
            {
                phi_r(i + index) = phi_r_slave(i);
            }
            else
            {
                phi_r(i + index) = -phi_r_slave(i);
            }
        }

        for (IndexType i = 0; i < phi_rs_master.size1(); i++)
        {
            for (IndexType j = 0; j < phi_rs_master.size2(); j++)
            {
                phi_rs(i, j) = phi_rs_master(i, j);
            }
        }

        SizeType index_1 = phi_rs_master.size1();
        SizeType index_2 = phi_rs_master.size2();
        for (IndexType i = 0; i < phi_rs_slave.size1(); i++)
        {
            for (IndexType j = 0; j < phi_rs_slave.size2(); j++)
            {
                if (opposite_direction_of_trims)
                    phi_rs(i + index_1, j + index_2) = phi_rs_slave(i, j);
                else
                    phi_rs(i + index_1, j + index_2) = -phi_rs_slave(i, j);
            }
        }
    } 

    void CouplingNitscheCondition::CalculateRotation(
        IndexType IntegrationPointIndex,
        const Matrix &rShapeFunctionGradientValues,
        Vector &phi_r,
        Matrix &phi_rs,
        array_1d<double, 2> &phi,
        array_1d<double, 3> &trim_tangent,
        const Vector &local_tangent,
        const bool master)
    {
        KRATOS_TRY

        const SizeType number_of_points = rShapeFunctionGradientValues.size1();
        
        // compute the initialize base vectors of master or slave 
        Vector g10 = ZeroVector(3);
        Vector g20 = ZeroVector(3);
        Vector g30 = ZeroVector(3);
        if (master)
        {
            for (SizeType i = 0; i < GetGeometry().GetGeometryPart(0).size(); ++i){
                g10[0] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 0);
                g10[1] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 0);
                g10[2] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 0);

                g20[0] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 1);
                g20[1] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 1);
                g20[2] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 1);

                MathUtils<double>::CrossProduct(g30, g10, g20);
                g30 = g30 / norm_2(g30);
            }
        }
        else
        {
            for (SizeType i = 0; i < GetGeometry().GetGeometryPart(1).size(); ++i){
                g10[0] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 0);
                g10[1] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 0);
                g10[2] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 0);

                g20[0] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 1);
                g20[1] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 1);
                g20[2] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 1);

                MathUtils<double>::CrossProduct(g30, g10, g20);
                g30 = g30 / norm_2(g30);
            }
        }

        // compute the actual base vectors of master or slave
        array_1d<double, 3> g1, g2, g3;
        Matrix J;

        if (master)
        {
            GetGeometry().GetGeometryPart(0).Jacobian(J, IntegrationPointIndex);
        }
        else 
        {
            GetGeometry().GetGeometryPart(1).Jacobian(J, IntegrationPointIndex);
        }

        g1 = column(J, 0);
        g2 = column(J, 1);

        MathUtils<double>::CrossProduct(g3, g1, g2);
        g3 = g3 / norm_2(g3);

        // compute the tangent (T2) and the normal (T1) to the boundary vector
        array_1d<double, 3> T1, T2;
        T2 = local_tangent[0] * g10 + local_tangent[1] * g20;
        trim_tangent = T2;
        MathUtils<double>::CrossProduct(T1, T2, g30);
        T2 = T2 / norm_2(T2);
        T1 = T1 / norm_2(T1);

        // compute the a3 displacement
        array_1d<double, 3> w = g3 - g30;
        array_1d<double, 3> sinus_omega_vector;
        MathUtils<double>::CrossProduct(sinus_omega_vector, g30, w);

        array_1d<double, 2> sinus_omega;
        sinus_omega(0) = inner_prod(sinus_omega_vector, T2);
        sinus_omega(1) = inner_prod(sinus_omega_vector, T1);

        array_1d<double, 3> omega;
        if (sinus_omega(0) > 1.0)
            sinus_omega(0) = 0.999999;
        if (sinus_omega(1) > 1.0)
            sinus_omega(1) = 0.999999;
        omega(0) = asin(sinus_omega(0));
        omega(1) = asin(sinus_omega(1));

        phi(0) = omega(0);
        phi(1) = omega(1);

        // compute variation of the a3 
        array_1d<double, 3> t3 = g3;
        array_1d<double, 3> tilde_t3; 
        MathUtils<double>::CrossProduct(tilde_t3, g1, g2);
        double length_t3 = norm_2(tilde_t3);

        std::vector<array_1d<double, 3>> t3_r(number_of_points * 3);
        std::vector<array_1d<double, 3>> tilde_3_r(number_of_points * 3);
        Vector line_t3_r = ZeroVector(number_of_points * 3);
        std::vector<array_1d<double, 3>> sinus_omega_r(number_of_points * 3);

        for (IndexType n = 0; n < number_of_points; n++)
        {
            for (IndexType i = 0; i < 3; i++)
            {
                int nb_dof = n * 3 + i;

                //variations of the basis vectors
                array_1d<double, 3> a1_r = ZeroVector(3);
                array_1d<double, 3> a2_r = ZeroVector(3);

                a1_r(i) = rShapeFunctionGradientValues(n, 0);
                a2_r(i) = rShapeFunctionGradientValues(n, 1);
                
                array_1d<double, 3> a1_r__g2, g1__a2_r = ZeroVector(3);
                MathUtils<double>::CrossProduct(a1_r__g2, a1_r, g2);
                MathUtils<double>::CrossProduct(g1__a2_r, g1, a2_r);

                //variation of the non normalized local vector
                tilde_3_r[nb_dof] = a1_r__g2 + g1__a2_r;
                line_t3_r[nb_dof] = inner_prod(t3, tilde_3_r[nb_dof]);
                t3_r[nb_dof] = tilde_3_r[nb_dof] / length_t3 - line_t3_r[nb_dof] * t3 / length_t3;

                MathUtils<double>::CrossProduct(sinus_omega_r[nb_dof], g30, t3_r[nb_dof]);
                phi_r(nb_dof) = 1.0 / sqrt(1.0 - pow(sinus_omega(0), 2))*inner_prod(sinus_omega_r[nb_dof], T2);
            }
        }

        for (IndexType n = 0; n < number_of_points; n++)
        {
            for (IndexType i = 0; i < 3; i++)
            {
                int nb_dof_n = n * 3 + i;
                
                //variations of the basis vectors
                array_1d<double, 3> a1_r_n = ZeroVector(3);
                array_1d<double, 3> a2_r_n = ZeroVector(3);

                a1_r_n(i) = rShapeFunctionGradientValues(n, 0);
                a2_r_n(i) = rShapeFunctionGradientValues(n, 1);

                for (IndexType m = 0; m < number_of_points; m++)
                {
                    for (IndexType j = 0; j < 3; j++)
                    {
                        int nb_dof_m = m * 3 + j;

                        //variations of the basis vectors
                        array_1d<double, 3> a1_r_m = ZeroVector(3);
                        array_1d<double, 3> a2_r_m = ZeroVector(3);

                        a1_r_m(j) = rShapeFunctionGradientValues(m, 0);
                        a2_r_m(j) = rShapeFunctionGradientValues(m, 1);

                        //variation of the non normalized local vector
                        array_1d<double, 3> a1_r_n__a2_r_m, a1_r_m__a2_r_n = ZeroVector(3);
                        MathUtils<double>::CrossProduct(a1_r_n__a2_r_m, a1_r_n, a2_r_m);
                        MathUtils<double>::CrossProduct(a1_r_m__a2_r_n, a1_r_m, a2_r_n);

                        array_1d<double, 3> tilde_t3_rs = a1_r_n__a2_r_m + a1_r_m__a2_r_n;
                        double line_t3_rs = inner_prod(t3_r[nb_dof_m], tilde_3_r[nb_dof_n]) + inner_prod(t3, tilde_t3_rs);

                        array_1d<double, 3> t3_rs = (tilde_t3_rs*length_t3 - line_t3_r[nb_dof_m] * tilde_3_r[nb_dof_n]) / pow(length_t3, 2)
                            - line_t3_rs * t3 / length_t3 - line_t3_r[nb_dof_n] * (t3_r[nb_dof_m] * length_t3 - line_t3_r[nb_dof_m] * t3) / pow(length_t3, 2);

                        array_1d<double, 3> sinus_omega_rs = ZeroVector(3);
                        MathUtils<double>::CrossProduct(sinus_omega_rs, g30, t3_rs);

                        phi_rs(n * 3 + i, m * 3 + j) = inner_prod(sinus_omega_rs, T2) / sqrt(1.0 - pow(sinus_omega(0), 2))
                            + inner_prod(sinus_omega_r[nb_dof_m], T2)*inner_prod(sinus_omega_r[nb_dof_n], T2)*sinus_omega(0) / pow(1.0
                                - pow(sinus_omega(0), 2), 1.5);
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    void CouplingNitscheCondition::CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe,
        const PatchType& rPatch) const
    {
        IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);

        const SizeType number_of_points = r_geometry.size();
        const SizeType working_space_dimension = r_geometry.WorkingSpaceDimension();
        Hessian.resize(working_space_dimension, working_space_dimension);
        Hessian = ZeroMatrix(working_space_dimension, working_space_dimension);

        for (IndexType k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3> coords = r_geometry[k].Coordinates();

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



