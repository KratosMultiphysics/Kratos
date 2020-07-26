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
#include "custom_conditions/nitsche_coupling_condition.h"

// Project includes

namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void NitscheCouplingCondition::Initialize()
    {
        KRATOS_TRY

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void NitscheCouplingCondition::InitializeMaterial()
    {
        KRATOS_TRY

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);
        const Properties& r_properties = GetProperties();
        const auto& r_N_master = r_geometry_master.ShapeFunctionsValues();
        const auto& r_N_slave = r_geometry_slave.ShapeFunctionsValues();

        const SizeType r_number_of_integration_points = r_geometry_master.IntegrationPointsNumber();

        //Constitutive Law initialisation
        if (mConstitutiveLawVector_master.size() != r_number_of_integration_points)
            mConstitutiveLawVector_master.resize(r_number_of_integration_points);
        if (mConstitutiveLawVector_slave.size() != r_number_of_integration_points)
            mConstitutiveLawVector_slave.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector_master.size(); ++point_number) {
            mConstitutiveLawVector_master[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector_master[point_number]->InitializeMaterial(r_properties, r_geometry_master, row(r_N_master, point_number));
        }
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector_slave.size(); ++point_number) {
            mConstitutiveLawVector_slave[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector_slave[point_number]->InitializeMaterial(r_properties, r_geometry_slave, row(r_N_slave, point_number));
        }

        KRATOS_CATCH("");
    }

    void NitscheCouplingCondition::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ConfigurationType& rConfiguration,
        const PatchType& rPatch
    )
    {
        IndexType GeometryPart;
        if (rPatch==PatchType::Master) GeometryPart = 0;
        if (rPatch==PatchType::Slave) GeometryPart = 1;

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
        if (rPatch==PatchType::Slave)
        {
            for (SizeType i=0;i<dimension*number_of_nodes;++i){
                current_displacement[i] = current_displacement_total[i+GetGeometry().GetGeometryPart(0).size()];
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
        rKinematicVariables.t = g1/norm_2(g1);
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
    void NitscheCouplingCondition::CalculateTransformation(
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

        // e * a_covariant
        Matrix G_hat = ZeroMatrix(2, 2);
        G_hat(0, 0) = inner_prod(e1, rKinematicVariables.a1);
        G_hat(0, 1) = inner_prod(e1, rKinematicVariables.a2);
        G_hat(1, 0) = inner_prod(e2, rKinematicVariables.a1);
        G_hat(1, 1) = inner_prod(e2, rKinematicVariables.a2);

        //Transformation matrix T_hat
        if (rT_hat.size1() != 3 && rT_hat.size2() != 3)
            rT_hat.resize(3, 3);
        noalias(rT_hat) = ZeroMatrix(3, 3);

        rT_hat(0, 0) = pow(G_hat(0, 0), 2);
        rT_hat(0, 1) = pow(G_hat(1, 0), 2);
        rT_hat(0, 2) = G_hat(0, 0) * G_hat(1, 0);

        rT_hat(1, 0) = pow(G_hat(0, 1), 2);
        rT_hat(1, 1) = pow(G_hat(1, 1), 2);
        rT_hat(1, 2) = G_hat(0, 1) * G_hat(1, 1);

        rT_hat(2, 0) = G_hat(0, 0) * G_hat(0, 1);
        rT_hat(2, 1) = G_hat(1, 0) * G_hat(1, 1);
        rT_hat(2, 2) = (G_hat(0, 0) * G_hat(1, 1) + G_hat(1, 0) * G_hat(0, 1));
    }

    void NitscheCouplingCondition::CalculateConstitutiveVariables(
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

            mConstitutiveLawVector_master[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        }
        if (rPatch==PatchType::Slave) 
        {
            array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector_slave[IntegrationPointIndex]);
            noalias(rThisConstitutiveVariablesMembrane.StrainVector) = prod(m_T_vector_slave[IntegrationPointIndex], strain_vector);

            // Constitive Matrices DMembrane
            rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
            rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
            rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

            mConstitutiveLawVector_slave[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        }

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
    }

    //Prestress Transformation Matrix
    void NitscheCouplingCondition::CalculateTransformationmatrixPrestress(
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

    void NitscheCouplingCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        const double penalty = GetProperties()[PENALTY_FACTOR];

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

        const IntegrationMethod integration_method_master = r_geometry_master.GetDefaultIntegrationMethod();
        const IntegrationMethod integration_method_slave = r_geometry_slave.GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_master = r_geometry_master.ShapeFunctionsLocalGradients(integration_method_master);
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_slave = r_geometry_slave.ShapeFunctionsLocalGradients(integration_method_slave);

        const SizeType r_number_of_integration_points = r_geometry_master.IntegrationPointsNumber();
        
        // Prepare memory
        if (m_A_ab_covariant_vector_master.size() != r_number_of_integration_points)
            m_A_ab_covariant_vector_master.resize(r_number_of_integration_points);
        if (m_A_ab_covariant_vector_slave.size() != r_number_of_integration_points)
            m_A_ab_covariant_vector_slave.resize(r_number_of_integration_points);
        if (m_dA_vector_master.size() != r_number_of_integration_points)
            m_dA_vector_master.resize(r_number_of_integration_points);
        if (m_dA_vector_slave.size() != r_number_of_integration_points)
            m_dA_vector_slave.resize(r_number_of_integration_points);
        if (m_T_vector_master.size() != r_number_of_integration_points)
            m_T_vector_master.resize(r_number_of_integration_points);
        if (m_T_vector_slave.size() != r_number_of_integration_points)
            m_T_vector_slave.resize(r_number_of_integration_points);
        if (m_T_hat_vector_master.size() != r_number_of_integration_points)
            m_T_hat_vector_master.resize(r_number_of_integration_points);
        if (m_T_hat_vector_slave.size() != r_number_of_integration_points)
            m_T_hat_vector_slave.resize(r_number_of_integration_points);
        if (m_reference_contravariant_base_master.size() != r_number_of_integration_points)
            m_reference_contravariant_base_master.resize(r_number_of_integration_points);
        if (m_reference_contravariant_base_slave.size() != r_number_of_integration_points)
            m_reference_contravariant_base_slave.resize(r_number_of_integration_points);
        if (m_n_contravariant_vector_master.size() != r_number_of_integration_points)
            m_n_contravariant_vector_master.resize(r_number_of_integration_points);
        if (m_n_contravariant_vector_slave.size() != r_number_of_integration_points)
            m_n_contravariant_vector_slave.resize(r_number_of_integration_points);

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
                r_geometry_master, GetProperties(), rCurrentProcessInfo);

            ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
                r_geometry_slave, GetProperties(), rCurrentProcessInfo);

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

            //Define Prestress
            array_1d<double, 3> prestress = GetProperties()[PRESTRESS]*GetProperties()[THICKNESS];

            PrestresstransVariables prestresstrans_variables_master(3);
            PrestresstransVariables prestresstrans_variables_slave(3);

            CalculateTransformationmatrixPrestress(
                kinematic_variables_master,
                prestresstrans_variables_master 
            );

            CalculateTransformationmatrixPrestress(
                kinematic_variables_slave,
                prestresstrans_variables_slave 
            );

            array_1d<double, 3> transformed_prestress_master = prod(prestresstrans_variables_master.Tpre, prestress);
            array_1d<double, 3> transformed_prestress_slave = prod(prestresstrans_variables_slave.Tpre, prestress);
            constitutive_variables_membrane_master.StressVector += transformed_prestress_master;
            constitutive_variables_membrane_slave.StressVector += transformed_prestress_slave;

            // calculate traction vectors
            array_1d<double, 3> traction_vector_master;
            array_1d<double, 3> traction_vector_slave;

            CalculateTraction(point_number, traction_vector_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateTraction(point_number, traction_vector_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            // calculate first variation of traction vectors
            Matrix d_traction_master = ZeroMatrix(3, 3*number_of_nodes_master);
            Matrix d_traction_slave = ZeroMatrix(3, 3*number_of_nodes_slave);

            CalculateDTraction(point_number, d_traction_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateDTraction(point_number, d_traction_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            Matrix d_traction = ZeroMatrix(3, mat_size);
            for (SizeType i=0;i<3 * number_of_nodes_master;++i){
                d_traction(0, i) = d_traction_master(0, i);
                d_traction(1, i) = d_traction_master(1, i);
                d_traction(2, i) = d_traction_master(2, i);
            }
            for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
                d_traction(0, i + 3 * number_of_nodes_master) = -d_traction_slave(0, i);
                d_traction(1, i + 3 * number_of_nodes_master) = -d_traction_slave(1, i);
                d_traction(2, i + 3 * number_of_nodes_master) = -d_traction_slave(2, i);
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

            //Get the displacement vectors of the previous iteration step
            Vector current_displacement_total = ZeroVector(mat_size);
            Vector current_displacement_master = ZeroVector(3 * number_of_nodes_master);
            Vector current_displacement_slave = ZeroVector(3 * number_of_nodes_slave);;

            GetValuesVector(current_displacement_total);

            for (SizeType i=0;i<3 * number_of_nodes_master;++i){
                current_displacement_master[i] = current_displacement_total[i];
            }
            for (SizeType i=0;i<3 * number_of_nodes_slave;++i){
                current_displacement_slave[i] = current_displacement_total[i + 3 * number_of_nodes_master];
            }

            array_1d<double, 3> disp_vector_master;
            array_1d<double, 3> disp_vector_slave;

            disp_vector_master = prod(r_N_master, current_displacement_master);
            disp_vector_slave = prod(r_N_slave, current_displacement_slave);

            //Compute the necessary products needed for the second variations of the traction vectors
            Matrix Pi_master = ZeroMatrix(3, 3);
            Matrix Pi_slave = ZeroMatrix(3, 3);

            CalculateDDTractionProduct(point_number, Pi_master, kinematic_variables_master, constitutive_variables_membrane_master, PatchType::Master);
            CalculateDDTractionProduct(point_number, Pi_slave, kinematic_variables_slave, constitutive_variables_membrane_slave, PatchType::Slave);

            array_1d<double, 3> dd_traction_product_vector_master;
            array_1d<double, 3> dd_traction_product_vector_slave;
            array_1d<double, 3> dd_traction_product_vector_master_slave;
            array_1d<double, 3> dd_traction_product_vector_slave_master;

            dd_traction_product_vector_master = prod(Pi_master, disp_vector_master);
            dd_traction_product_vector_slave = prod(Pi_slave, disp_vector_slave);
            dd_traction_product_vector_master_slave = prod(Pi_slave, disp_vector_master);
            dd_traction_product_vector_slave_master = prod(Pi_master, disp_vector_slave);

            // calculate second variation of traction vectors

            Matrix dd_traction_master = ZeroMatrix(3 * number_of_nodes_master, 3 * number_of_nodes_master);
            Matrix dd_traction_slave = ZeroMatrix(3 * number_of_nodes_slave, 3 * number_of_nodes_slave);

            CalculateDDTraction(point_number, dd_traction_master, kinematic_variables_master, d_traction_master, disp_vector_master, disp_vector_slave, dd_traction_product_vector_master, dd_traction_product_vector_master_slave, PatchType::Master);
            CalculateDDTraction(point_number, dd_traction_slave, kinematic_variables_slave, d_traction_slave, disp_vector_master, disp_vector_slave, dd_traction_product_vector_slave, dd_traction_product_vector_slave_master, PatchType::Slave);
            
            //Penalty part & RHS
            Matrix H = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                IndexType index = 3 * i;
                //if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    H(0, index) = N_master(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    H(1, index + 1) = N_master(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    H(2, index + 2) = N_master(point_number, i);
            }
            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                IndexType index = 3 * (i + number_of_nodes_master);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    H(0, index) = -N_slave(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    H(1, index + 1) = -N_slave(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    H(2, index + 2) = -N_slave(point_number, i);
            }

            //first variation part
            Matrix HH = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                IndexType index = 3 * i;
                //if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    HH(0, index) = N_master(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    HH(1, index + 1) = N_master(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    HH(2, index + 2) = N_master(point_number, i);
            }
            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                IndexType index = 3 * (i + number_of_nodes_master);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    HH(0, index) = N_slave(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    HH(1, index + 1) = N_slave(point_number, i);
                //if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    HH(2, index + 2) = N_slave(point_number, i);
            }

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = r_geometry_master.DeterminantOfJacobian(point_number);

            const double gammaTilde = 0.5;

            // Assembly
            if (CalculateStiffnessMatrixFlag) {

                noalias(rLeftHandSideMatrix) += (prod(trans(d_traction), HH) + prod(trans(HH), d_traction))
                    * integration_weight * determinat_jacobian * -gammaTilde;
                     
                for (IndexType i = 0; i < number_of_nodes_master; i++)
                {
                    for (IndexType j = 0; j < number_of_nodes_master; j++)
                    {
                        rLeftHandSideMatrix(i, j) += dd_traction_master(i, j);
                    }
                }

                for (IndexType i = 0; i < number_of_nodes_slave; i++)
                {
                    for (IndexType j = 0; j < number_of_nodes_slave; j++)
                    {
                        rLeftHandSideMatrix(i + number_of_nodes_slave, j + number_of_nodes_slave) += dd_traction_slave(i, j);
                    }
                }

                noalias(rLeftHandSideMatrix) += prod(trans(H), H)
                    * integration_weight * determinat_jacobian * penalty;
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
                    * integration_weight * determinat_jacobian * -gammaTilde;
                noalias(rRightHandSideVector) -= (prod(trans(d_traction), disp_vector_master) - prod(trans(d_traction), disp_vector_slave))
                    * integration_weight * determinat_jacobian * -gammaTilde;
                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
                    * integration_weight * determinat_jacobian * penalty;
            }
        }

        KRATOS_CATCH("")
    }

    void NitscheCouplingCondition::CalculateHessian(
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

    void NitscheCouplingCondition::CalculateTraction(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch)
    {
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;

        if (rPatch==PatchType::Master)
        {
            stress_vector_covariant = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
        }
        if (rPatch==PatchType::Master)
        {
            stress_vector_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
        }

        // Compute the stress components
        Matrix Palphabeta = ZeroMatrix(2, 2);
        Palphabeta(0,0) = stress_vector_covariant[0];
        Palphabeta(1,1) = stress_vector_covariant[1];
        Palphabeta(0,1) = stress_vector_covariant[2];
        Palphabeta(1,0) = Palphabeta(0,1);

        array_1d<double, 2> P_n;
        P_n = prod(Palphabeta, rActualKinematic.n_contravariant);
        
        // Compute the traction vectors
        rTraction[0] = rActualKinematic.a1[0]*P_n[0] + rActualKinematic.a2[0]*P_n[1];
        rTraction[1] = rActualKinematic.a1[1]*P_n[0] + rActualKinematic.a2[1]*P_n[1];
        rTraction[2] = rActualKinematic.a1[2]*P_n[0] + rActualKinematic.a2[2]*P_n[1];
    }

    void NitscheCouplingCondition::CalculateDDTractionProduct(
        IndexType IntegrationPointIndex,
        Matrix& rPi,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch)
    {
        Matrix a_u = ZeroMatrix(3, 3); //covariant times normal vector

        for (IndexType r = 0; r < 3; r++)
        {
            a_u (0, r) = rActualKinematic.a1[r] * rActualKinematic.n_contravariant[0];
            a_u (1, r) = rActualKinematic.a2[r] * rActualKinematic.n_contravariant[1];
            a_u (2, r) = rActualKinematic.a1[r] * rActualKinematic.n_contravariant[1] + rActualKinematic.a2[r] * rActualKinematic.n_contravariant[0];
        }

        if (rPatch==PatchType::Master)
        {
            rPi = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.ConstitutiveMatrix);
            rPi = prod(rPi, m_T_vector_master[IntegrationPointIndex]);
        }
        if (rPatch==PatchType::Master)
        {
            rPi = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.ConstitutiveMatrix);
            rPi = prod(rPi, m_T_vector_slave[IntegrationPointIndex]);
        }

        rPi = prod(a_u, rPi);
    }

    void NitscheCouplingCondition::CalculateDTraction(
        IndexType IntegrationPointIndex,
        Matrix& rDTraction,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch)
    {
        IndexType GeometryPart;
        if (rPatch==PatchType::Master) GeometryPart = 0;
        if (rPatch==PatchType::Slave) GeometryPart = 1;

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
        if (rPatch==PatchType::Slave)
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
        Matrix dN_cartesian = ZeroMatrix(3, mat_size);
        dN_cartesian = prod(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix,dE_cartesian);

        //Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
        Matrix dN_curvilinear = ZeroMatrix(3, mat_size);

        if (rPatch==PatchType::Master)
        {
            dN_curvilinear = prod(m_T_hat_vector_master[IntegrationPointIndex], dN_cartesian);
        }
        if (rPatch==PatchType::Slave)
        {
            dN_curvilinear = prod(m_T_hat_vector_slave[IntegrationPointIndex], dN_cartesian);
        }

        //Compute the first variation of the traction vectors
        Matrix n_a = ZeroMatrix(3, 3); //normal vector times covariant current

        n_a(0, 0) = rActualKinematic.n_contravariant[0]*rActualKinematic.a1[0];
        n_a(1, 0) = rActualKinematic.n_contravariant[0]*rActualKinematic.a1[1];
        n_a(2, 0) = rActualKinematic.n_contravariant[0]*rActualKinematic.a1[2];
        n_a(0, 1) = rActualKinematic.n_contravariant[1]*rActualKinematic.a2[0];
        n_a(1, 1) = rActualKinematic.n_contravariant[1]*rActualKinematic.a2[1];
        n_a(2, 1) = rActualKinematic.n_contravariant[1]*rActualKinematic.a2[2];
        n_a(0, 2) = rActualKinematic.n_contravariant[1]*rActualKinematic.a1[0] + rActualKinematic.n_contravariant[0]*rActualKinematic.a2[0];
        n_a(0, 2) = rActualKinematic.n_contravariant[1]*rActualKinematic.a1[1] + rActualKinematic.n_contravariant[0]*rActualKinematic.a2[1];
        n_a(0, 2) = rActualKinematic.n_contravariant[1]*rActualKinematic.a1[2] + rActualKinematic.n_contravariant[0]*rActualKinematic.a2[2];

        rDTraction = prod(n_a, dN_curvilinear); 

        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;

        if (rPatch==PatchType::Master)
        {
            stress_vector_covariant = prod(m_T_hat_vector_master[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
        }
        if (rPatch==PatchType::Slave)
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

        rDTraction += r_DN_Dxi*(rActualKinematic.n_contravariant[0]*stress_vector_covariant[0]+rActualKinematic.n_contravariant[1]*stress_vector_covariant[2])+
                      r_DN_Deta*(rActualKinematic.n_contravariant[1]*stress_vector_covariant[1]+rActualKinematic.n_contravariant[0]*stress_vector_covariant[2]);
    }

    void NitscheCouplingCondition::CalculateDDTraction(
        IndexType IntegrationPointIndex,
        Matrix& rDDTraction,
        const KinematicVariables& rActualKinematic,
        Matrix& rDTraction, 
        array_1d<double, 3>& rDispMaster,
        array_1d<double, 3>& rDispSlave,
        array_1d<double, 3>& rDDTractionProduct,
        array_1d<double, 3>& rDDTractionProductMasterSlave,
        const PatchType& rPatch)
    {
        IndexType GeometryPart;
        if (rPatch==PatchType::Master) GeometryPart = 0;
        if (rPatch==PatchType::Slave) GeometryPart = 1;

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

        Vector dNCovariant1 = ZeroVector(mat_size);
        Vector dNCovariant2 = ZeroVector(mat_size);
        Vector dNCovariant3 = ZeroVector(mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            dNCovariant1(r) = rDTraction(0, r);
            dNCovariant2(r) = rDTraction(1, r);
            dNCovariant3(r) = rDTraction(2, r);
        }

        Matrix dispdNCovariant1_master = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(1,:)
        Matrix dispdNCovariant2_master = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(2,:)
        Matrix dispdNCovariant3_master = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(3,:)

        Matrix dispdNCovariant1_slave = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(1,:)
        Matrix dispdNCovariant2_slave = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(2,:)
        Matrix dispdNCovariant3_slave = ZeroMatrix(3, mat_size); //dispVctI*dNCovariantI(3,:)

        for (IndexType r = 0; r < mat_size; r++)
        {
            dispdNCovariant1_master(0, r) = rDispMaster(0)*dNCovariant1(r);
            dispdNCovariant1_master(1, r) = rDispMaster(1)*dNCovariant1(r);
            dispdNCovariant1_master(2, r) = rDispMaster(2)*dNCovariant1(r);

            dispdNCovariant2_master(0, r) = rDispMaster(0)*dNCovariant2(r);
            dispdNCovariant2_master(1, r) = rDispMaster(1)*dNCovariant2(r);
            dispdNCovariant2_master(2, r) = rDispMaster(2)*dNCovariant2(r);

            dispdNCovariant3_master(0, r) = rDispMaster(0)*dNCovariant3(r);
            dispdNCovariant3_master(1, r) = rDispMaster(1)*dNCovariant3(r);
            dispdNCovariant3_master(2, r) = rDispMaster(2)*dNCovariant3(r);

            dispdNCovariant1_slave(0, r) = rDispSlave(0)*dNCovariant1(r);
            dispdNCovariant1_slave(1, r) = rDispSlave(1)*dNCovariant1(r);
            dispdNCovariant1_slave(2, r) = rDispSlave(2)*dNCovariant1(r);

            dispdNCovariant2_slave(0, r) = rDispSlave(0)*dNCovariant2(r);
            dispdNCovariant2_slave(1, r) = rDispSlave(1)*dNCovariant2(r);
            dispdNCovariant2_slave(2, r) = rDispSlave(2)*dNCovariant2(r);

            dispdNCovariant3_slave(0, r) = rDispSlave(0)*dNCovariant3(r);
            dispdNCovariant3_slave(1, r) = rDispSlave(1)*dNCovariant3(r);
            dispdNCovariant3_slave(2, r) = rDispSlave(2)*dNCovariant3(r);
        }

        if (rPatch==PatchType::Slave)
        {
            dispdNCovariant1_master *= -1.0;
            dispdNCovariant2_master *= -1.0;
            dispdNCovariant3_master *= -1.0;

            dispdNCovariant1_slave *= -1.0;
            dispdNCovariant2_slave *= -1.0;
            dispdNCovariant3_slave *= -1.0;
        }

        rDDTraction += prod(trans(r_DN_Dxi), r_DN_Dxi)*rDDTractionProduct(0);
        rDDTraction += prod(trans(r_DN_Deta), r_DN_Deta)*rDDTractionProduct(1);
        rDDTraction += 0.5*(prod(trans(r_DN_Dxi), r_DN_Deta) + prod(trans(r_DN_Deta), r_DN_Dxi))*rDDTractionProduct(2);

        rDDTraction += prod(trans(r_DN_Dxi), dispdNCovariant1_master)*rActualKinematic.n_contravariant(0);
        rDDTraction += prod(trans(r_DN_Deta), dispdNCovariant2_master)*rActualKinematic.n_contravariant(1);
        rDDTraction += 0.5*(prod(trans(r_DN_Dxi), dispdNCovariant3_master)*rActualKinematic.n_contravariant(0) + prod(trans(r_DN_Deta), dispdNCovariant3_master)*rActualKinematic.n_contravariant(1));
        rDDTraction += prod(trans(dispdNCovariant1_master), r_DN_Dxi)*rActualKinematic.n_contravariant(0);
        rDDTraction += prod(trans(dispdNCovariant2_master), r_DN_Deta)*rActualKinematic.n_contravariant(1);
        rDDTraction += 0.5*(prod(trans(dispdNCovariant3_master), r_DN_Dxi)*rActualKinematic.n_contravariant(0) + prod(trans(dispdNCovariant3_master), r_DN_Deta)*rActualKinematic.n_contravariant(1));
        
        rDDTraction -= prod(trans(r_DN_Dxi), r_DN_Dxi)*rDDTractionProductMasterSlave(0);
        rDDTraction -= prod(trans(r_DN_Deta), r_DN_Deta)*rDDTractionProductMasterSlave(1);
        rDDTraction -= 0.5*(prod(trans(r_DN_Dxi), r_DN_Deta) + prod(trans(r_DN_Deta), r_DN_Dxi))*rDDTractionProductMasterSlave(2);

        rDDTraction -= prod(trans(r_DN_Dxi), dispdNCovariant1_slave)*rActualKinematic.n_contravariant(0);
        rDDTraction -= prod(trans(r_DN_Deta), dispdNCovariant2_slave)*rActualKinematic.n_contravariant(1);
        rDDTraction -= 0.5*(prod(trans(r_DN_Dxi), dispdNCovariant3_slave)*rActualKinematic.n_contravariant(0) + prod(trans(r_DN_Deta), dispdNCovariant3_slave)*rActualKinematic.n_contravariant(1));
        rDDTraction -= prod(trans(dispdNCovariant1_slave), r_DN_Dxi)*rActualKinematic.n_contravariant(0);
        rDDTraction -= prod(trans(dispdNCovariant2_slave), r_DN_Deta)*rActualKinematic.n_contravariant(1);
        rDDTraction -= 0.5*(prod(trans(dispdNCovariant3_slave), r_DN_Dxi)*rActualKinematic.n_contravariant(0) + prod(trans(dispdNCovariant3_slave), r_DN_Deta)*rActualKinematic.n_contravariant(1));
    }

    void NitscheCouplingCondition::GetValuesVector(
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

    void NitscheCouplingCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
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

    void NitscheCouplingCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
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


