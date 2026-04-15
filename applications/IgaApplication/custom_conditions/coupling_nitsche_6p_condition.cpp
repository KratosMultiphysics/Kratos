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
#include "custom_conditions/coupling_nitsche_6p_condition.h"

namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void CouplingNitsche6pCondition::CalculateKinematics(
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

        Vector current_displacement_total = ZeroVector(2*dimension*(GetGeometry().GetGeometryPart(0).size()+GetGeometry().GetGeometryPart(1).size()));
        Vector current_displacement = ZeroVector(2*dimension*number_of_nodes);

        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement_total);

        if (rPatch==PatchType::Master)
        {
            for (SizeType i=0;i<2*dimension*number_of_nodes;++i){
                current_displacement[i] = current_displacement_total[i];
            }
        }
        else
        {
            for (SizeType i=0;i<2*dimension*number_of_nodes;++i){
                current_displacement[i] = current_displacement_total[i+GetGeometry().GetGeometryPart(0).size()*6];
            }
        }

        for (SizeType i=0;i<number_of_nodes;++i){
            g1[0] += (r_geometry.GetPoint( i ).X0()+current_displacement[i*2*dimension]) * rShapeFunctionGradientValues(i, 0);
            g1[1] += (r_geometry.GetPoint( i ).Y0()+current_displacement[(i*2*dimension)+1]) * rShapeFunctionGradientValues(i, 0);
            g1[2] += (r_geometry.GetPoint( i ).Z0()+current_displacement[(i*2*dimension)+2]) * rShapeFunctionGradientValues(i, 0);

            g2[0] += (r_geometry.GetPoint( i ).X0()+current_displacement[i*2*dimension]) * rShapeFunctionGradientValues(i, 1);
            g2[1] += (r_geometry.GetPoint( i ).Y0()+current_displacement[(i*2*dimension)+1]) * rShapeFunctionGradientValues(i, 1);
            g2[2] += (r_geometry.GetPoint( i ).Z0()+current_displacement[(i*2*dimension)+2]) * rShapeFunctionGradientValues(i, 1);
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
        rKinematicVariables.n_contravariant[2] = rKinematicVariables.a3[0] * rKinematicVariables.n[0] + rKinematicVariables.a3[1] * rKinematicVariables.n[1] +rKinematicVariables.a3[2] * rKinematicVariables.n[2];
    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void CouplingNitsche6pCondition::CalculateTransformation(
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
        

        //Transformation matrix T from contravariant to local cartesian basis
        // e * a_contravariant
     
        double l_a1 = norm_2(rKinematicVariables.a1);
        array_1d<double, 3> e1 = rKinematicVariables.a1 / l_a1;
        array_1d<double, 3> e3 =  rKinematicVariables.a3;
        array_1d<double, 3> e2;
        MathUtils<double>::CrossProduct(e2, e3, e1);
        
        //Transformation matrix 
        if (rT.size1() != 6 && rT.size2() != 6)                                                                 
            rT.resize(6, 6);
        noalias(rT) = ZeroMatrix(6, 6);

        if (rT_hat.size1() != 6 && rT_hat.size2() != 6)                                                                 
            rT_hat.resize(6, 6);
        noalias(rT_hat) = ZeroMatrix(6, 6);

        for (std::size_t i = 0; i < 3; ++i)
        {
            std::size_t j = (i + 1) % 3;

            rT_hat(i, 0) = e1[i] * e1[i];
            rT_hat(i, 1) = e2[i] * e2[i];
            rT_hat(i, 2) = e3[i] * e3[i];
            rT_hat(i, 3) = 2 * e1[i] * e2[i];
            rT_hat(i, 4) = 2 * e2[i] * e3[i];
            rT_hat(i, 5) = 2 * e1[i] * e3[i];

            rT_hat(i + 3, 0) = e1[i] * e1[j];
            rT_hat(i + 3, 1) = e2[i] * e2[j];
            rT_hat(i + 3, 2) = e3[i] * e3[j];
            rT_hat(i + 3, 3) = (e1[i] * e2[j]) + (e2[i] * e1[j]);
            rT_hat(i + 3, 4) = (e2[i] * e3[j]) + (e3[i] * e2[j]);
            rT_hat(i + 3, 5) = (e1[i] * e3[j]) + (e3[i] * e1[j]);
        }
    }

    void CouplingNitsche6pCondition::CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariables,
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
            array_1d<double, 6> strain_vector(6, 0.0);
        //To do: compute strain vector
        noalias(rThisConstitutiveVariables.StrainVector) = prod(m_T_vector_master[IntegrationPointIndex], strain_vector);

        // Constitive Matrices D
        rValues.SetStrainVector(rThisConstitutiveVariables.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.ConstitutiveMatrix); //this is an ouput parameter

        const double nu = this->GetProperties()[POISSON_RATIO];
        const double Emodul = this->GetProperties()[YOUNG_MODULUS];
        double lambda = Emodul / (1.0 - nu * nu);
        double Gmodul = Emodul / (2.0 * (1.0 + nu));
        
        rThisConstitutiveVariables.ConstitutiveMatrix(0, 0) = lambda;
        rThisConstitutiveVariables.ConstitutiveMatrix(0, 1) = lambda * nu;
        rThisConstitutiveVariables.ConstitutiveMatrix(1, 0) = lambda * nu;
        rThisConstitutiveVariables.ConstitutiveMatrix(1, 1) = lambda;
        rThisConstitutiveVariables.ConstitutiveMatrix(3, 3) = lambda * (1 - nu) / 2;
        rThisConstitutiveVariables.ConstitutiveMatrix(4, 4) = Gmodul * 5.0 / 6.0;
        rThisConstitutiveVariables.ConstitutiveMatrix(5, 5) = Gmodul * 5.0 / 6.0;

        //Local Cartesian Stresses
        noalias(rThisConstitutiveVariables.StressVector) = prod(
            trans(rThisConstitutiveVariables.ConstitutiveMatrix), rThisConstitutiveVariables.StrainVector);
        }
        else
        {
             array_1d<double, 6> strain_vector(6, 0.0);
        //To do: compute strain vector
        noalias(rThisConstitutiveVariables.StrainVector) = prod(m_T_vector_slave[IntegrationPointIndex], strain_vector);

        // Constitive Matrices D
        rValues.SetStrainVector(rThisConstitutiveVariables.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.ConstitutiveMatrix); //this is an ouput parameter

        const double nu = this->GetProperties()[POISSON_RATIO];
        const double Emodul = this->GetProperties()[YOUNG_MODULUS];
        double lambda = Emodul / (1.0 - nu * nu);
        double Gmodul = Emodul / (2.0 * (1.0 + nu));
        
        rThisConstitutiveVariables.ConstitutiveMatrix(0, 0) = lambda;
        rThisConstitutiveVariables.ConstitutiveMatrix(0, 1) = lambda * nu;
        rThisConstitutiveVariables.ConstitutiveMatrix(1, 0) = lambda * nu;
        rThisConstitutiveVariables.ConstitutiveMatrix(1, 1) = lambda;
        rThisConstitutiveVariables.ConstitutiveMatrix(3, 3) = lambda * (1 - nu) / 2;
        rThisConstitutiveVariables.ConstitutiveMatrix(4, 4) = Gmodul * 5.0 / 6.0;
        rThisConstitutiveVariables.ConstitutiveMatrix(5, 5) = Gmodul * 5.0 / 6.0;

        //Local Cartesian Stresses
        noalias(rThisConstitutiveVariables.StressVector) = prod(
            trans(rThisConstitutiveVariables.ConstitutiveMatrix), rThisConstitutiveVariables.StrainVector);
        }

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariables.StressVector) = prod(
            trans(rThisConstitutiveVariables.ConstitutiveMatrix), rThisConstitutiveVariables.StrainVector);
    }

    //Prestress Transformation Matrix
    void CouplingNitsche6pCondition::CalculateTransformationPrestress(
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

    void CouplingNitsche6pCondition::CalculateAll(
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

        const SizeType mat_size = 6 * (number_of_nodes_master + number_of_nodes_slave);

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

            ConstitutiveVariables constitutive_variables_master(6);
            ConstitutiveVariables constitutive_variables_slave(6);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_master,
                constitutive_variables_master,
                constitutive_law_parameters_master,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Master);
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables_slave,
                constitutive_variables_master,
                constitutive_law_parameters_slave,
                ConstitutiveLaw::StressMeasure_PK2,
                PatchType::Slave);



             for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {
                // Initialization
                Matrix B_master= ZeroMatrix(6, number_of_nodes_master*6);
                Matrix dn_master = ZeroMatrix(3, 3);
                Matrix DN_De_Jn_master = ZeroMatrix(number_of_nodes_master,3);
                Matrix J_inv_master = ZeroMatrix(3, 3);
                double area_master = 0.0;
                Matrix B_slave = ZeroMatrix(6, number_of_nodes_slave*6);
                Matrix dn_slave = ZeroMatrix(3, 3);
                Matrix DN_De_Jn_slave = ZeroMatrix(number_of_nodes_slave,3);
                Matrix J_inv_slave = ZeroMatrix(3, 3);
                double area_slave = 0.0;
                mZeta = mGaussIntegrationThickness.zeta(Gauss_index);

                CalculateJn(point_number, kinematic_variables_master, mZeta, DN_De_Jn_master, J_inv_master, dn_master, area_master,PatchType::Master);
                CalculateJn(point_number, kinematic_variables_slave, mZeta, DN_De_Jn_slave, J_inv_slave, dn_slave, area_slave,PatchType::Slave);
                CalculateB(point_number, B_master, mZeta, DN_De_Jn_master, J_inv_master, dn_master, kinematic_variables_master,PatchType::Master);
                CalculateB(point_number, B_slave, mZeta, DN_De_Jn_slave, J_inv_slave, dn_slave, kinematic_variables_slave,PatchType::Slave);

                

               
            }

             // calculate traction vectors
            array_1d<double, 6> traction_vector_master;
            array_1d<double, 6> traction_vector_slave;

            
            Matrix r_UnitN_master = ZeroMatrix(3, 6);
            Matrix r_UnitN_slave = ZeroMatrix(3, 6);
            BuildNormalOperator6p(point_number, r_UnitN_master, PatchType::Master);
            BuildNormalOperator6p(point_number, r_UnitN_slave, PatchType::Slave);




            // CalculateTraction(point_number, traction_vector_master, kinematic_variables_master, constitutive_variables_master, PatchType::Master);
            // CalculateTraction(point_number, traction_vector_slave, kinematic_variables_slave, constitutive_variables_slave, PatchType::Slave);

            // calculate the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
            // Matrix first_variations_stress_covariant_master = ZeroMatrix(6, 6*number_of_nodes_master);
            // Matrix first_variations_stress_covariant_slave = ZeroMatrix(6, 6*number_of_nodes_slave);

            // CalculateFirstVariationStressCovariant(point_number, first_variations_stress_covariant_master, kinematic_variables_master, constitutive_variables_master, PatchType::Master);
            // CalculateFirstVariationStressCovariant(point_number, first_variations_stress_covariant_slave, kinematic_variables_slave, constitutive_variables_slave, PatchType::Slave);

            // calculate first variation of traction vectors
            // Matrix first_variations_traction_master = ZeroMatrix(6, 6*number_of_nodes_master);
            // Matrix first_variations_traction_slave = ZeroMatrix(6, 6*number_of_nodes_slave);

            // CalculateFirstVariationTraction(point_number, first_variations_traction_master, first_variations_stress_covariant_master, kinematic_variables_master, constitutive_variables_master, PatchType::Master);
            // CalculateFirstVariationTraction(point_number, first_variations_traction_slave, first_variations_stress_covariant_slave, kinematic_variables_slave, constitutive_variables_slave, PatchType::Slave);

            // Matrix first_variations_traction = ZeroMatrix(6, mat_size);
            // for (SizeType i=0;i<6 * number_of_nodes_master;++i){
            //     first_variations_traction(0, i) = first_variations_traction_master(0, i);
            //     first_variations_traction(1, i) = first_variations_traction_master(1, i);
            //     first_variations_traction(2, i) = first_variations_traction_master(2, i);
            //     first_variations_traction(3, i) = first_variations_traction_master(3, i);
            //     first_variations_traction(4, i) = first_variations_traction_master(4, i);
            //     first_variations_traction(5, i) = first_variations_traction_master(5, i);
            // }
            // for (SizeType i=0;i<6 * number_of_nodes_slave;++i){
            //     first_variations_traction(0, i + 6 * number_of_nodes_master) = -first_variations_traction_slave(0, i);
            //     first_variations_traction(1, i + 6 * number_of_nodes_master) = -first_variations_traction_slave(1, i);
            //     first_variations_traction(2, i + 6 * number_of_nodes_master) = -first_variations_traction_slave(2, i);
            //     first_variations_traction(3, i + 6 * number_of_nodes_master) = -first_variations_traction_slave(3, i);
            //     first_variations_traction(4, i + 6 * number_of_nodes_master) = -first_variations_traction_slave(4, i);
            //     first_variations_traction(5, i + 6 * number_of_nodes_master) = -first_variations_traction_slave(5, i);
            // }



            
            //Compute the NURBS basis functions
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            Matrix r_N_master = ZeroMatrix(6, 6 * number_of_nodes_master);
            Matrix r_N_slave = ZeroMatrix(6, 6* number_of_nodes_slave);

            for (IndexType r = 0; r < number_of_nodes_master; r++)
            {
                r_N_master(0, 6 * r) = N_master(0, r);
                r_N_master(1, 6 * r + 1) = N_master(0, r);
                r_N_master(2, 6 * r + 2) = N_master(0, r);
                r_N_master(3, 6 * r) = N_master(0, r);
                r_N_master(4, 6 * r + 1) = N_master(0, r);
                r_N_master(5, 6 * r + 2) = N_master(0, r);
            }
            for (IndexType r = 0; r < number_of_nodes_slave; r++)
            {
                r_N_slave(0, 6 * r) = N_slave(0, r);
                r_N_slave(1, 6 * r + 1) = N_slave(0, r);
                r_N_slave(2, 6 * r + 2) = N_slave(0, r);
                r_N_slave(3, 6 * r) = N_slave(0, r);
                r_N_slave(4, 6 * r + 1) = N_slave(0, r);
                r_N_slave(5, 6 * r + 2) = N_slave(0, r);
            }

            //Get the displacement vectors of the previous iteration step
            Vector current_displacement_total = ZeroVector(mat_size);
            Vector current_displacement_master = ZeroVector(6 * number_of_nodes_master);
            Vector current_displacement_slave = ZeroVector(6 * number_of_nodes_slave);

            GetValuesVector(current_displacement_total);

            for (SizeType i=0;i<6 * number_of_nodes_master;++i){
                current_displacement_master[i] = current_displacement_total[i];
            }
            for (SizeType i=0;i<6 * number_of_nodes_slave;++i){
                current_displacement_slave[i] = current_displacement_total[i + 6 * number_of_nodes_master];
            }

            array_1d<double, 6> displacement_vector_master;
            array_1d<double, 6> displacement_vector_slave;

            displacement_vector_master = prod(r_N_master, current_displacement_master);
            displacement_vector_slave = prod(r_N_slave, current_displacement_slave);


            

         
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
             

            // if (CalculateResidualVectorFlag) {

            //     Vector u(mat_size);
            //     for (IndexType i = 0; i < number_of_nodes_master; i++)
            //     {
            //         const array_1d<double, 3> disp = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT);
            //         IndexType index = 3 * i;
            //         u[index]     = disp[0];
            //         u[index + 1] = disp[1];
            //         u[index + 2] = disp[2];
            //     }
            //     for (IndexType i = 0; i < number_of_nodes_slave; i++)
            //     {
            //         const array_1d<double, 3> disp = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT);
            //         IndexType index = 3 * (i + number_of_nodes_master);
            //         u[index]     = disp[0];
            //         u[index + 1] = disp[1];
            //         u[index + 2] = disp[2];
            //     }
                
            //     noalias(rRightHandSideVector) -= (prod(trans(H), traction_vector_master) - prod(trans(H), traction_vector_slave))
            //         * integration_weight * determinant_jacobian * -gammaTilde;
            //     noalias(rRightHandSideVector) -= (prod(trans(first_variations_traction), displacement_vector_master) - prod(trans(first_variations_traction), displacement_vector_slave))
            //         * integration_weight * determinant_jacobian * -gammaTilde;   
            //     noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
            //         * integration_weight * determinant_jacobian * stabilization_parameter;
            // }
        }
        KRATOS_CATCH("")
    }

    void CouplingNitsche6pCondition::DeterminantOfJacobianInitial(
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

    void CouplingNitsche6pCondition::CalculateNitscheStabilizationMatrix(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const double beta = GetProperties()[NITSCHE_STABILIZATION_FACTOR];

    const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
    const auto& r_geometry_slave  = GetGeometry().GetGeometryPart(1);

    const SizeType number_of_nodes_master = r_geometry_master.size();
    const SizeType number_of_nodes_slave  = r_geometry_slave.size();

    const SizeType mat_size = 6 * (number_of_nodes_master + number_of_nodes_slave);
    const SizeType master_size = 6 * number_of_nodes_master;
    const SizeType slave_size  = 6 * number_of_nodes_slave;

    if (rLeftHandSideMatrix.size1() != mat_size) {
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    if (rRightHandSideVector.size() != mat_size) {
        rRightHandSideVector.resize(mat_size, false);
    }
    rRightHandSideVector = ZeroVector(mat_size);

    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry_master.IntegrationPoints();

    Vector determinant_jacobian_vector_initial(integration_points.size());
    DeterminantOfJacobianInitial(r_geometry_master, determinant_jacobian_vector_initial);

    const IntegrationMethod integration_method_master = r_geometry_master.GetDefaultIntegrationMethod();
    const IntegrationMethod integration_method_slave  = r_geometry_slave.GetDefaultIntegrationMethod();

    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_master =
        r_geometry_master.ShapeFunctionsLocalGradients(integration_method_master);
    const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_slave =
        r_geometry_slave.ShapeFunctionsLocalGradients(integration_method_slave);

    const SizeType n_gp_master = r_geometry_master.IntegrationPointsNumber();
    const SizeType n_gp_slave  = r_geometry_slave.IntegrationPointsNumber();

    if (m_A_ab_covariant_vector_master.size() != n_gp_master)
        m_A_ab_covariant_vector_master.resize(n_gp_master);
    if (m_A_ab_covariant_vector_slave.size() != n_gp_slave)
        m_A_ab_covariant_vector_slave.resize(n_gp_slave);
    if (m_dA_vector_master.size() != n_gp_master)
        m_dA_vector_master.resize(n_gp_master);
    if (m_dA_vector_slave.size() != n_gp_slave)
        m_dA_vector_slave.resize(n_gp_slave);
    if (m_T_vector_master.size() != n_gp_master)
        m_T_vector_master.resize(n_gp_master);
    if (m_T_vector_slave.size() != n_gp_slave)
        m_T_vector_slave.resize(n_gp_slave);
    if (m_T_hat_vector_master.size() != n_gp_master)
        m_T_hat_vector_master.resize(n_gp_master);
    if (m_T_hat_vector_slave.size() != n_gp_slave)
        m_T_hat_vector_slave.resize(n_gp_slave);
    if (m_reference_contravariant_base_master.size() != n_gp_master)
        m_reference_contravariant_base_master.resize(n_gp_master);
    if (m_reference_contravariant_base_slave.size() != n_gp_slave)
        m_reference_contravariant_base_slave.resize(n_gp_slave);
    if (m_n_contravariant_vector_master.size() != n_gp_master)
        m_n_contravariant_vector_master.resize(n_gp_master);
    if (m_n_contravariant_vector_slave.size() != n_gp_slave)
        m_n_contravariant_vector_slave.resize(n_gp_slave);

    array_1d<double, 3> characteristic_length_master;
    array_1d<double, 3> characteristic_length_slave;
    r_geometry_master.GetGeometryParent(0).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)
        ->Calculate(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_master);
    r_geometry_slave.GetGeometryParent(0).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)
        ->Calculate(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_slave);

    const double tol_basic = 0.01;
    const double characteristic_length =
        std::max(norm_2(characteristic_length_master), norm_2(characteristic_length_slave));
    const double tol_surface_normal = tol_basic * characteristic_length * characteristic_length;

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        const Matrix& shape_functions_gradients_i_master = r_shape_functions_gradients_master[point_number];
        const Matrix& shape_functions_gradients_i_slave  = r_shape_functions_gradients_slave[point_number];

        KinematicVariables kv_ref_master(r_geometry_master.WorkingSpaceDimension());
        KinematicVariables kv_ref_slave(r_geometry_slave.WorkingSpaceDimension());

        CalculateKinematics(point_number, kv_ref_master, shape_functions_gradients_i_master,
                            ConfigurationType::Reference, PatchType::Master);
        CalculateKinematics(point_number, kv_ref_slave, shape_functions_gradients_i_slave,
                            ConfigurationType::Reference, PatchType::Slave);

        m_A_ab_covariant_vector_master[point_number] = kv_ref_master.a_ab_covariant;
        m_A_ab_covariant_vector_slave[point_number]  = kv_ref_slave.a_ab_covariant;
        m_dA_vector_master[point_number] = kv_ref_master.dA;
        m_dA_vector_slave[point_number]  = kv_ref_slave.dA;
        m_n_contravariant_vector_master[point_number] = kv_ref_master.n_contravariant;
        m_n_contravariant_vector_slave[point_number]  = kv_ref_slave.n_contravariant;

        CalculateTransformation(kv_ref_master,
                                m_T_vector_master[point_number],
                                m_T_hat_vector_master[point_number],
                                m_reference_contravariant_base_master[point_number]);

        CalculateTransformation(kv_ref_slave,
                                m_T_vector_slave[point_number],
                                m_T_hat_vector_slave[point_number],
                                m_reference_contravariant_base_slave[point_number]);

        KinematicVariables kv_master(r_geometry_master.WorkingSpaceDimension());
        KinematicVariables kv_slave(r_geometry_slave.WorkingSpaceDimension());

        CalculateKinematics(point_number, kv_master, shape_functions_gradients_i_master,
                            ConfigurationType::Current, PatchType::Master);
        CalculateKinematics(point_number, kv_slave, shape_functions_gradients_i_slave,
                            ConfigurationType::Current, PatchType::Slave);

        ConstitutiveLaw::Parameters cl_params_master(
            r_geometry_master, GetProperties().GetSubProperties().front(), rCurrentProcessInfo);
        ConstitutiveLaw::Parameters cl_params_slave(
            r_geometry_slave, GetProperties().GetSubProperties().back(), rCurrentProcessInfo);

        ConstitutiveVariables cv_master(6);
        ConstitutiveVariables cv_slave(6);

        CalculateConstitutiveVariables(point_number, kv_master, cv_master,
                                       cl_params_master, ConstitutiveLaw::StressMeasure_PK2,
                                       PatchType::Master);

        CalculateConstitutiveVariables(point_number, kv_slave, cv_slave,
                                       cl_params_slave, ConstitutiveLaw::StressMeasure_PK2,
                                       PatchType::Slave);

        Matrix N_master = r_geometry_master.ShapeFunctionsValues();
        Matrix N_slave  = r_geometry_slave.ShapeFunctionsValues();

        Matrix R_master = ZeroMatrix(6, master_size);
        Matrix R_slave  = ZeroMatrix(6, slave_size);

        for (IndexType r = 0; r < number_of_nodes_master; ++r) {
            for (IndexType a = 0; a < 6; ++a) {
                R_master(a, 6 * r + a) = N_master(point_number, r);
            }
        }

        for (IndexType r = 0; r < number_of_nodes_slave; ++r) {
            for (IndexType a = 0; a < 6; ++a) {
                R_slave(a, 6 * r + a) = N_slave(point_number, r);
            }
        }

        Matrix K11 = ZeroMatrix(master_size, master_size);
        Matrix K22 = ZeroMatrix(slave_size, slave_size);
        Matrix K12 = ZeroMatrix(master_size, slave_size);
        Matrix K21 = ZeroMatrix(slave_size, master_size);

        for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; ++Gauss_index)
        {
            Matrix B_master = ZeroMatrix(6, master_size);
            Matrix dn_master = ZeroMatrix(3, 3);
            Matrix DN_De_Jn_master = ZeroMatrix(number_of_nodes_master, 3);
            Matrix J_inv_master = ZeroMatrix(3, 3);
            double area_master = 0.0;

            Matrix B_slave = ZeroMatrix(6, slave_size);
            Matrix dn_slave = ZeroMatrix(3, 3);
            Matrix DN_De_Jn_slave = ZeroMatrix(number_of_nodes_slave, 3);
            Matrix J_inv_slave = ZeroMatrix(3, 3);
            double area_slave = 0.0;

            mZeta = mGaussIntegrationThickness.zeta(Gauss_index);

            CalculateJn(point_number, kv_master, mZeta, DN_De_Jn_master, J_inv_master, dn_master, area_master, PatchType::Master);
            CalculateJn(point_number, kv_slave,  mZeta, DN_De_Jn_slave,  J_inv_slave,  dn_slave,  area_slave,  PatchType::Slave);

            CalculateB(point_number, B_master, mZeta, DN_De_Jn_master, J_inv_master, dn_master, kv_master, PatchType::Master);
            CalculateB(point_number, B_slave,  mZeta, DN_De_Jn_slave,  J_inv_slave,  dn_slave,  kv_slave,  PatchType::Slave);

            Matrix C_master = cv_master.ConstitutiveMatrix;
            Matrix C_slave  = cv_slave.ConstitutiveMatrix;

            Matrix n_master = ZeroMatrix(6, 6);
            Matrix n_slave  = ZeroMatrix(6, 6);

            BuildNormalOperator6p(point_number, n_master, PatchType::Master);
            BuildNormalOperator6p(point_number, n_slave, PatchType::Slave);

            Matrix CB_master = ZeroMatrix(6, master_size);
            noalias(CB_master) = prod(C_master, B_master);

            Matrix G_master = ZeroMatrix(6, master_size);
            noalias(G_master) = prod(n_master, CB_master);
            Matrix CB_slave = ZeroMatrix(6, slave_size);
            noalias(CB_slave) = prod(C_slave, B_slave);

            Matrix G_slave = ZeroMatrix(6, slave_size);
            noalias(G_slave) = prod(n_slave, CB_slave);

            const double w =
                integration_points[point_number].Weight();

            noalias(K11) += (
                -0.5 * prod(trans(R_master), G_master)
                -0.5 * prod(trans(G_master), R_master)
                + beta * prod(trans(R_master), R_master)
            ) * w;

            noalias(K22) += (
                -0.5 * prod(trans(R_slave), G_slave)
                -0.5 * prod(trans(G_slave), R_slave)
                + beta * prod(trans(R_slave), R_slave)
            ) * w;

            noalias(K12) += (
                -beta * prod(trans(R_master), R_slave)
                -0.5 * prod(trans(R_master), G_slave)
                +0.5 * prod(trans(G_master), R_slave)
            ) * w;

            noalias(K21) += (
                -beta * prod(trans(R_slave), R_master)
                -0.5 * prod(trans(R_slave), G_master)
                +0.5 * prod(trans(G_slave), R_master)
            ) * w;
        }

        if (norm_2(kv_ref_master.a3_tilde) > tol_surface_normal &&
            norm_2(kv_ref_slave.a3_tilde) > tol_surface_normal)
        {
            for (IndexType i = 0; i < master_size; ++i) {
                for (IndexType j = 0; j < master_size; ++j) {
                    rLeftHandSideMatrix(i, j) += K11(i, j);
                }
            }

            for (IndexType i = 0; i < master_size; ++i) {
                for (IndexType j = 0; j < slave_size; ++j) {
                    rLeftHandSideMatrix(i, j + master_size) += K12(i, j);
                }
            }

            for (IndexType i = 0; i < slave_size; ++i) {
                for (IndexType j = 0; j < master_size; ++j) {
                    rLeftHandSideMatrix(i + master_size, j) += K21(i, j);
                }
            }

            for (IndexType i = 0; i < slave_size; ++i) {
                for (IndexType j = 0; j < slave_size; ++j) {
                    rLeftHandSideMatrix(i + master_size, j + master_size) += K22(i, j);
                }
            }
        }
    }

    KRATOS_CATCH("")
}
    
    void CouplingNitsche6pCondition::CalculateJn(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        double zeta,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        double& area, 
        const PatchType& rPatch) const
    {

        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        double thickness = this->GetProperties().GetValue(THICKNESS);  

        Matrix J;
        GetGeometry().Jacobian(J, IntegrationPointIndex);
        const std::size_t number_of_control_points = GetGeometry().size();

        rKinematicVariables.a1 = column(J, 0);
        rKinematicVariables.a2 = column(J, 1);
        MathUtils<double>::CrossProduct(rKinematicVariables.a3_tilde, rKinematicVariables.a1, rKinematicVariables.a2);
        rKinematicVariables.dA = norm_2(rKinematicVariables.a3_tilde);
        noalias(rKinematicVariables.a3) = rKinematicVariables.a3_tilde / rKinematicVariables.dA;


        Matrix da3 = ZeroMatrix(3, 3);
        double inv_dA = 1 / rKinematicVariables.dA;
        double inv_dA3 = 1 / std::pow(rKinematicVariables.dA, 3);

        //compute da3
        array_1d<double, 3> da1_d1;
        array_1d<double, 3> da1_d2; //da1_d2 = da2_d1
        array_1d<double, 3> da2_d2;

        noalias(da1_d1) = ZeroVector(3);
        noalias(da1_d2) = ZeroVector(3);
        noalias(da2_d2) = ZeroVector(3);

        for (std::size_t i=0;i<number_of_control_points;++i){
            da1_d1[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 0);
            da1_d1[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 0);
            da1_d1[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 0);

            da1_d2[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 1);
            da1_d2[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 1);
            da1_d2[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 1);

            da2_d2[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 2);
            da2_d2[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 2);
            da2_d2[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 2);
        }

        array_1d<double, 3> da3_tilde_d1;
        array_1d<double, 3> da3_tilde_d1_1;
        array_1d<double, 3> da3_tilde_d1_2;
        array_1d<double, 3> da3_tilde_d2;
        array_1d<double, 3> da3_tilde_d2_1;
        array_1d<double, 3> da3_tilde_d2_2;

        MathUtils<double>::CrossProduct(da3_tilde_d1_1, da1_d1, rKinematicVariables.a2);
        MathUtils<double>::CrossProduct(da3_tilde_d1_2, rKinematicVariables.a1, da1_d2);
        da3_tilde_d1 = da3_tilde_d1_1 + da3_tilde_d1_2;

        MathUtils<double>::CrossProduct(da3_tilde_d2_1, da1_d2, rKinematicVariables.a2);
        MathUtils<double>::CrossProduct(da3_tilde_d2_2, rKinematicVariables.a1, da2_d2);
        da3_tilde_d2 = da3_tilde_d2_1 + da3_tilde_d2_2;

        for (IndexType j = 0; j < 3; j++)
        {
            dn(0, j) = da3_tilde_d1[j] * inv_dA - rKinematicVariables.a3_tilde[j] * inner_prod(rKinematicVariables.a3_tilde, da3_tilde_d1) * inv_dA3;
            dn(1, j) = da3_tilde_d2[j] * inv_dA - rKinematicVariables.a3_tilde[j] * inner_prod(rKinematicVariables.a3_tilde, da3_tilde_d2) * inv_dA3;
            dn(2, j) = 0.0;
        }

        Matrix Jn = ZeroMatrix(3,3);
        for (int i = 0; i < 3; i++) {
            Jn(0, i) = rKinematicVariables.a1[i] + (thickness/2) * zeta * dn(0, i); 
            Jn(1, i) = rKinematicVariables.a2[i] + (thickness/2) * zeta * dn(1, i) ; 
            Jn(2, i) = rKinematicVariables.a3[i] * (thickness/2) ;
        }

        // Jn(0,0) = 4.0; //Jn must be computed in the reference configuration! Warning: This is just a hard coded test

        Matrix Jn_inv = ZeroMatrix(3,3);
        double det_Jn = 0.0;
        MathUtils<double>::InvertMatrix(Jn, Jn_inv, det_Jn);

        Matrix new_DN_De = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De(i, 0) = r_DN_De(i, 0);  // Copy first column
            new_DN_De(i, 1) = r_DN_De(i, 1);  // Copy second column
            new_DN_De(i, 2) = 0.0;            // Set third column to zero
        }

        DN_De_Jn = trans(prod(Jn_inv, trans(new_DN_De)));  // DN_De_Jn = Jn_inv * r_DN_De that has 3 columns
        J_inv = Jn_inv;
        area = det_Jn;
    }
    
     void CouplingNitsche6pCondition::CalculateB(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        double zeta,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic,
        const PatchType& rPatch) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 6;
        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        double thickness = this->GetProperties().GetValue(THICKNESS);   

        //Bending part
        Matrix DN_De_Jn_bending = ZeroMatrix(number_of_control_points,3);
        Matrix new_DN_De_bending = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De_bending(i, 0) = r_DN_De(i, 0) * (thickness/2) * zeta;  // Copy first column
            new_DN_De_bending(i, 1) = r_DN_De(i, 1) * (thickness/2) * zeta;  // Copy second column
            new_DN_De_bending(i, 2) = 0.0;            // Set third column to zero
        }
        DN_De_Jn_bending = trans(prod(J_inv, trans(new_DN_De_bending)));
                                                
        // y hat 
        const double y1= rActualKinematic.a3[0];
        const double y2= rActualKinematic.a3[1];
        const double y3= rActualKinematic.a3[2];

        Matrix da3 = ZeroMatrix(3, 3);
        Matrix Dn = ZeroMatrix(3, 3);
        Matrix b = ZeroMatrix(3, mat_size);

        //compute da3
        array_1d<double, 3> da1_d1;
        array_1d<double, 3> da1_d2; //da1_d2 = da2_d1
        array_1d<double, 3> da2_d2;

        array_1d<double, 3> da3_tilde_d1;
        array_1d<double, 3> da3_tilde_d1_1;
        array_1d<double, 3> da3_tilde_d1_2;
        array_1d<double, 3> da3_tilde_d2;
        array_1d<double, 3> da3_tilde_d2_1;
        array_1d<double, 3> da3_tilde_d2_2;

        if (rB.size1() != 6 || rB.size2() != mat_size)
            rB.resize(6, mat_size);
        noalias(rB) = ZeroMatrix(6, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            // zeta  Derivatives 
            const double dzetadx= J_inv(0,2);
            const double dzetady= J_inv(1,2);
            const double dzetadz= J_inv(2,2); 
          
            Dn = prod(J_inv, dn);

            // y hat Derivative w.r.t x
            const double dy1x= Dn(0,0);
            const double dy2x= Dn(0,1);
            const double dy3x= Dn(0,2);

            // y hat Derivative w.r.t y
            const double dy1y= Dn(1,0);
            const double dy2y= Dn(1,1);
            const double dy3y= Dn(1,2);

            // y hat Derivative w.r.t z
            const double dy1z= Dn(2,0);
            const double dy2z= Dn(2,1);
            const double dy3z= Dn(2,2);

            // Compute rB
            IndexType index = i * 6;
            
            rB(0, index)     = DN_De_Jn(i, 0);
            rB(0, index + 1) = 0;
            rB(0, index + 2) = 0;

            rB(0, index + 3) = 0;
            rB(0, index + 4) = (DN_De_Jn_bending(i, 0) * y3) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy3x + y3 * dzetadx));
            rB(0, index + 5) = - ((DN_De_Jn_bending(i, 0) * y2) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy2x + dzetadx * y2)));

            rB(1, index)     = 0;
            rB(1, index + 1) = DN_De_Jn(i, 1);
            rB(1, index + 2) = 0;

            rB(1, index + 3) = - ((DN_De_Jn_bending(i, 1) * y3) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy3y + dzetady * y3))); 
            rB(1, index + 4) = 0;
            rB(1, index + 5) = (DN_De_Jn_bending(i, 1) * y1) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy1y + dzetady * y1)); 

            rB(2, index)     = 0;
            rB(2, index + 1) = 0;
            rB(2, index + 2) = DN_De_Jn(i, 2);

            rB(2, index + 3) = (DN_De_Jn_bending(i, 2)   * y2) + (r_N(IntegrationPointIndex, i)  * (thickness/2) * (zeta * dy2z + dzetadz * y2));  
            rB(2, index + 4) = - ((DN_De_Jn_bending(i, 2) * y1)  + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta  * dy1z + dzetadz *y1))); 
            rB(2, index + 5) = 0;
            
            rB(3, index)     = DN_De_Jn(i, 1);                    
            rB(3, index + 1) = DN_De_Jn(i, 0);  
            rB(3, index + 2) = 0;

            rB(3, index + 3) = - ((DN_De_Jn_bending(i, 0) * y3) +(r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy3x + dzetadx * y3)));
            rB(3, index + 4) = (DN_De_Jn_bending(i, 1) * y3) +(r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy3y + dzetady * y3));
            rB(3, index + 5) = ((DN_De_Jn_bending(i, 0) * y1) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy1x + dzetadx * y1))) - ((DN_De_Jn_bending(i, 1) * y2)+ (r_N(IntegrationPointIndex, i) * (thickness/2) *  (zeta * dy2y + dzetady * y2))); 

            rB(4, index)     = 0;
            rB(4, index + 1) = DN_De_Jn(i, 2);
            rB(4, index + 2) = DN_De_Jn(i, 1);

            rB(4, index + 3) = ((DN_De_Jn_bending(i, 1) * y2) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy2y + dzetady * y2)))  - ((DN_De_Jn_bending(i, 2) * y3)+ (r_N(IntegrationPointIndex, i)  * (thickness/2) * (zeta *dy3z + dzetadz * y3))); 
            rB(4, index + 4) = - ((DN_De_Jn_bending(i, 1) * y1) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy1y + dzetady * y1))); 
            rB(4, index + 5) = (DN_De_Jn_bending(i, 2) *  y1 ) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy1z + dzetadz * y1)); 

            rB(5, index)   = DN_De_Jn (i,2);                    
            rB(5, index + 1) = 0;  
            rB(5, index + 2) = DN_De_Jn(i, 0);

            rB(5, index + 3) = ((DN_De_Jn_bending(i, 0)  * y2) + (r_N(IntegrationPointIndex, i) * (thickness/2) * ( zeta * dy2x +  dzetadx * y2)) ); 
            rB(5, index + 4) = ((DN_De_Jn_bending (i,2)  * y3) + (r_N(IntegrationPointIndex, i)   * (thickness/2) * ( zeta * dy3z + dzetadz * y3 ))) - ((DN_De_Jn_bending(i, 0) * y1) + (r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta  * dy1x + dzetadx * y1) )); 
            rB(5, index + 5) = - ((DN_De_Jn_bending (i,2)  * y2 )+(r_N(IntegrationPointIndex, i) * (thickness/2) * (zeta * dy2z + dzetadz * y2)));  
        }
    }
    
    
    void CouplingNitsche6pCondition::BuildNormalOperator6p(
    IndexType IntegrationPointIndex,
    Matrix& r_n, 
    const PatchType& rPatch)
    {
    if (r_n.size1() != 3 || r_n.size2() != 6)
        r_n.resize(3, 6, false);

    array_1d<double, 3> n_contravariant_vector;

        if (rPatch==PatchType::Master)
        {
            
            n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        }
        else
        {
           
            n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        }


    r_n(0,0) = n_contravariant_vector[0];
    r_n(0,3) = n_contravariant_vector[1];
    r_n(0,4) = n_contravariant_vector[2];

    r_n(1,1) = n_contravariant_vector[0];
    r_n(1,3) = n_contravariant_vector[1];
    r_n(1,5) = n_contravariant_vector[2];

    r_n(2,4) = n_contravariant_vector[0];
    r_n(2,5) = n_contravariant_vector[1];
  
    }


    // void CouplingNitsche6pCondition::CalculateTraction(
    //     IndexType IntegrationPointIndex,
    //     array_1d<double, 6>& rTraction,
    //     const KinematicVariables& rActualKinematic,
    //     ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
    //     const PatchType& rPatch)
    // {
    //     // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
    //     array_1d<double, 6> stress_vector_covariant;
    //     array_1d<double, 3> n_contravariant_vector;

    //     if (rPatch==PatchType::Master)
    //     {
    //         stress_vector_covariant = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
    //         n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
    //     }
    //     else
    //     {
    //         stress_vector_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
    //         n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
    //     }

    //     // Compute the stress components
    //     Matrix Palphabeta = ZeroMatrix(2, 2);
    //     Palphabeta(0,0) = stress_vector_covariant[0];
    //     Palphabeta(1,1) = stress_vector_covariant[1];
    //     Palphabeta(0,1) = stress_vector_covariant[2];
    //     Palphabeta(1,0) = Palphabeta(0,1);
        
    //     // Compute the traction vectors
    //     rTraction[0] = rActualKinematic.a1[0]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
    //                  + rActualKinematic.a2[0]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
    //     rTraction[1] = rActualKinematic.a1[1]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
    //                  + rActualKinematic.a2[1]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
    //     rTraction[2] = rActualKinematic.a1[2]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
    //                  + rActualKinematic.a2[2]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
    // }

    // void CouplingNitsche6pCondition::CalculateFirstVariationStressCovariant(
    //     IndexType IntegrationPointIndex,
    //     Matrix& rFirstVariationStressCovariant,
    //     const KinematicVariables& rActualKinematic,
    //     ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
    //     const PatchType& rPatch)
    // {
    //      IndexType GeometryPart = (rPatch==PatchType::Master) ? 0 : 1;
    //     const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);
        
    //     const SizeType number_of_control_points = r_geometry.size();
    //     const SizeType mat_size = number_of_control_points * 3;
        
    //     const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

    //     //Compute the first variation of the Green-Lagrange strains
    //     Matrix dE_cartesian = ZeroMatrix(3, mat_size);
    //     Matrix T_patch = ZeroMatrix(3, 3);

    //     if (rPatch==PatchType::Master)
    //     {
    //         T_patch = m_T_vector_master[IntegrationPointIndex];
    //     }
    //     else
    //     {
    //         T_patch = m_T_vector_slave[IntegrationPointIndex];
    //     }

    //     for (IndexType r = 0; r < mat_size; r++)
    //     {
    //         // local node number kr and dof direction dirr
    //         IndexType kr = r / 3;
    //         IndexType dirr = r % 3;

    //         array_1d<double, 3> dE_curvilinear;
    //         // strain
    //         dE_curvilinear[0] = r_DN_De(kr, 0)*rActualKinematic.a1(dirr);
    //         dE_curvilinear[1] = r_DN_De(kr, 1)*rActualKinematic.a2(dirr);
    //         dE_curvilinear[2] = 0.5*(r_DN_De(kr, 0)*rActualKinematic.a2(dirr) + rActualKinematic.a1(dirr)*r_DN_De(kr, 1));

    //         dE_cartesian(0, r) = T_patch(0, 0)*dE_curvilinear[0] + T_patch(0, 1)*dE_curvilinear[1] + T_patch(0, 2)*dE_curvilinear[2];
    //         dE_cartesian(1, r) = T_patch(1, 0)*dE_curvilinear[0] + T_patch(1, 1)*dE_curvilinear[1] + T_patch(1, 2)*dE_curvilinear[2];
    //         dE_cartesian(2, r) = T_patch(2, 0)*dE_curvilinear[0] + T_patch(2, 1)*dE_curvilinear[1] + T_patch(2, 2)*dE_curvilinear[2];
    //     }

    //     //Compute the first variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
    //     Matrix first_variations_stress_cartesian = ZeroMatrix(3, mat_size);
    //     first_variations_stress_cartesian = prod(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix,dE_cartesian);

    //     //Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
    //     if (rPatch==PatchType::Master)
    //     {
    //         rFirstVariationStressCovariant = prod(m_T_hat_vector_master[IntegrationPointIndex], first_variations_stress_cartesian);
    //     }
    //     else
    //     {
    //         rFirstVariationStressCovariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], first_variations_stress_cartesian);
    //     }
    // }

    //     //Compute the first variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
    //     Matrix first_variations_stress_cartesian = ZeroMatrix(6, mat_size);
    //     first_variations_stress_cartesian = prod(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix,dE_cartesian);

    //     //Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
    //     if (rPatch==PatchType::Master)
    //     {
    //         rFirstVariationStressCovariant = prod(m_T_hat_vector_master[IntegrationPointIndex], first_variations_stress_cartesian);
    //     }
    //     else
    //     {
    //         rFirstVariationStressCovariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], first_variations_stress_cartesian);
    //     }
    // }
    
    void CouplingNitsche6pCondition::CalculateFirstVariationTraction(
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

    // void CouplingNitsche6pCondition::CalculateSecondVariationTractionProduct(
    //     IndexType IntegrationPointIndex,
    //     Matrix& rPi,
    //     const KinematicVariables& rActualKinematic,
    //     ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
    //     const PatchType& rPatch)
    // {
    //     array_1d<double, 2> n_contravariant_vector;

    //     if (rPatch==PatchType::Master)
    //     {
    //         rPi = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.ConstitutiveMatrix);
    //         rPi = prod(rPi, m_T_vector_master[IntegrationPointIndex]);

    //         n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
    //     }
    //     else
    //     {
    //         rPi = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.ConstitutiveMatrix);
    //         rPi = prod(rPi, m_T_vector_slave[IntegrationPointIndex]);

    //         n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
    //     }

    //     normal vector * covariant base vector
    //     Matrix n_a = ZeroMatrix(3, 3);

    //     for (IndexType r = 0; r < 3; r++)
    //     {
    //         n_a (r, 0) = rActualKinematic.a1[r] * n_contravariant_vector[0];
    //         n_a (r, 1) = rActualKinematic.a2[r] * n_contravariant_vector[1];
    //         n_a (r, 2) = rActualKinematic.a1[r] * n_contravariant_vector[1] + rActualKinematic.a2[r] * n_contravariant_vector[0];
    //     }

    //     rPi = prod(n_a, rPi);
    // }



    void CouplingNitsche6pCondition::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_control_points_master = r_geometry_master.size();
        const SizeType number_of_control_points_slave = r_geometry_slave.size();
        const SizeType mat_size = (number_of_control_points_master + number_of_control_points_slave) * 6;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points_master; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const array_1d<double, 3 >& rotation = r_geometry_master[i].FastGetSolutionStepValue(ROTATION, Step);
            IndexType index = i * 6;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
            rValues[index + 3] = rotation[0];
            rValues[index + 4] = rotation[1];
            rValues[index + 5] = rotation[2];
            
        }

        for (IndexType i = 0; i < number_of_control_points_slave; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const array_1d<double, 3 >& rotation = r_geometry_slave[i].FastGetSolutionStepValue(ROTATION, Step);
            IndexType index = 6 * (i + number_of_control_points_master) ;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
            rValues[index + 3] = rotation[0];
            rValues[index + 4] = rotation[1];
            rValues[index + 5] = rotation[2];
            
            
        }
    }

    void CouplingNitsche6pCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        if (rResult.size() != 6 * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(6 * (number_of_nodes_master + number_of_nodes_slave), false);

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry_master[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
            
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = 3 * (i + number_of_nodes_master);
            const auto& r_node = r_geometry_slave[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void CouplingNitsche6pCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * (number_of_nodes_master + number_of_nodes_slave));

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_X));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Y));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Z));
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_X));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Y));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos



