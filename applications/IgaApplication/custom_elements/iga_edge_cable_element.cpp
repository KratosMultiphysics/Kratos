//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Ricky Aristio
//                   Pia Halbich
//                   Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/iga_edge_cable_element.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void IgaEdgeCableElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (mReferenceArea.size() != r_number_of_integration_points)
            mReferenceArea.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
        {
            const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(point_number);

            array_1d<double, 3> tangents;
            GetGeometry().Calculate(LOCAL_TANGENT, tangents);
            const SizeType number_of_nodes = r_geometry.size();
            
            Matrix J = ZeroMatrix(3, 2);
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                for (unsigned int k = 0; k<3; k++)
                {
                    for (unsigned int m = 0; m<2; m++)
                    {
                        J(k, m) += (r_geometry[i]).Coordinates()[k] * r_DN_De(i, m);
                    }
                }
            }

            mReferenceArea[point_number] = norm_2(column(J, 0) * tangents[0] + column(J, 1) * tangents[1]);
        }



        KRATOS_CATCH("")
    }

    array_1d<double, 3> IgaEdgeCableElement::GetActualBaseVector(const Matrix& r_DN_De, const ConfigurationType& rConfiguration) 
    {
        const GeometryType& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        array_1d<double, 3> tangents;
        GetGeometry().Calculate(LOCAL_TANGENT, tangents);

        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        array_1d<double, 3> actual_base_vector = ZeroVector(dimension);

        Vector current_displacement = ZeroVector(dimension*number_of_nodes);
        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);

        //basis vectors g1 and g2
        Vector g1 = ZeroVector(dimension);
        Vector g2 = ZeroVector(dimension);

        for (SizeType i=0;i<number_of_nodes;++i){
            g1[0] += (GetGeometry().GetPoint( i ).X0()+current_displacement[i*dimension]) * r_DN_De(i, 0);
            g1[1] += (GetGeometry().GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * r_DN_De(i, 0);
            g1[2] += (GetGeometry().GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * r_DN_De(i, 0);

            g2[0] += (GetGeometry().GetPoint( i ).X0()+current_displacement[i*dimension]) * r_DN_De(i, 1);
            g2[1] += (GetGeometry().GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * r_DN_De(i, 1);
            g2[2] += (GetGeometry().GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * r_DN_De(i, 1);
        }

        actual_base_vector = g1 * tangents[0] + g2 * tangents[1];

        return actual_base_vector;
    }

    ///@}
    ///@name Assembly
    ///@{

    void IgaEdgeCableElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (mReferenceBaseVector.size() != r_number_of_integration_points)
            mReferenceBaseVector.resize(r_number_of_integration_points);
        
        //get properties
        array_1d<double, 3> tangents;
        GetGeometry().Calculate(LOCAL_TANGENT, tangents);
        const double E = GetProperties()[YOUNG_MODULUS];
        const double A = GetProperties()[CROSS_AREA];
        const double prestress = GetProperties()[PRESTRESS_CAUCHY];

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) 
        {   
            // get integration data
            const double& integration_weight = r_integration_points[point_number].Weight();
            const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(point_number);

            mReferenceBaseVector[point_number] = GetActualBaseVector(r_DN_De, ConfigurationType::Reference);
            const double reference_a = norm_2(mReferenceBaseVector[point_number]);    

            // compute base vectors
            const array_1d<double, 3> actual_base_vector = GetActualBaseVector(r_DN_De, ConfigurationType::Current);
    
            // green-lagrange strain
            const double e11_membrane = 0.5 * (inner_prod(actual_base_vector, actual_base_vector) - inner_prod(mReferenceBaseVector[point_number], mReferenceBaseVector[point_number]));

            // normal force reference_aa
            const double s11_membrane = prestress * A + e11_membrane * A * E / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

            for (IndexType r = 0; r < mat_size; r++)
            {
                // local node number kr and dof direction dirr
                IndexType kr = r / 3;
                IndexType dirr = r % 3;

                const double epsilon_var_r = actual_base_vector[dirr] *
                    (r_DN_De(kr, 0) * tangents[0] 
                    + r_DN_De(kr, 1) * tangents[1]) / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

                if (CalculateStiffnessMatrixFlag) {
                    for (IndexType s = 0; s < mat_size; s++)
                    {
                        // local node number ks and dof direction dirs
                        IndexType ks = s / 3;
                        IndexType dirs = s % 3;

                        const double epsilon_var_s =
                            actual_base_vector[dirs] *
                            (r_DN_De(ks, 0) * tangents[0]
                            + r_DN_De(ks, 1) * tangents[1])
                            / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

                        rLeftHandSideMatrix(r, s) = E * A * epsilon_var_r * epsilon_var_s * reference_a * integration_weight;

                        if (dirr == dirs) {
                            const double epsilon_var_rs =
                            (r_DN_De(kr, 0) * tangents[0] + r_DN_De(kr, 1) * tangents[1]) *
                            (r_DN_De(ks, 0) * tangents[0] + r_DN_De(ks, 1) * tangents[1]) /inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);
                     
                            rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs * reference_a * integration_weight; 
                        }
                    }
                }
                if (CalculateResidualVectorFlag) {
                    rRightHandSideVector[r] = -s11_membrane * epsilon_var_r * reference_a * integration_weight;
                }
            }
        }
        KRATOS_CATCH("");
    }

    void IgaEdgeCableElement::CalculateInitialStiffnessMatrix(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (mReferenceBaseVector.size() != r_number_of_integration_points)
            mReferenceBaseVector.resize(r_number_of_integration_points);
        
        //get properties
        array_1d<double, 3> tangents;
        GetGeometry().Calculate(LOCAL_TANGENT, tangents);
        const double E = GetProperties()[YOUNG_MODULUS];
        const double A = GetProperties()[CROSS_AREA];
        const double prestress = GetProperties()[PRESTRESS_CAUCHY];

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) 
        {   
            // get integration data
            const double& integration_weight = r_integration_points[point_number].Weight();
            const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(point_number);

            mReferenceBaseVector[point_number] = GetActualBaseVector(r_DN_De, ConfigurationType::Reference);
            const double reference_a = norm_2(mReferenceBaseVector[point_number]);    

            // normal forcereference_aa
            const double s11_membrane = prestress * A;

            for (IndexType r = 0; r < mat_size; r++)
            {
                // local node number kr and dof direction dirr
                IndexType kr = r / 3;
                IndexType dirr = r % 3;

                const double epsilon_var_r = mReferenceBaseVector[point_number][dirr] *
                    (r_DN_De(kr, 0) * tangents[0] 
                    + r_DN_De(kr, 1) * tangents[1]) / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

                for (IndexType s = 0; s < mat_size; s++)
                {
                    // local node number ks and dof direction dirs
                    IndexType ks = s / 3;
                    IndexType dirs = s % 3;

                    const double epsilon_var_s =
                        mReferenceBaseVector[point_number][dirs] *
                        (r_DN_De(ks, 0) * tangents[0]
                        + r_DN_De(ks, 1) * tangents[1])
                        / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

                    rLeftHandSideMatrix(r, s) = E * A * epsilon_var_r * epsilon_var_s * reference_a * integration_weight;

                    if (dirr == dirs) {
                        const double epsilon_var_rs =
                        (r_DN_De(kr, 0) * tangents[0] + r_DN_De(kr, 1) * tangents[1]) *
                        (r_DN_De(ks, 0) * tangents[0] + r_DN_De(ks, 1) * tangents[1]) /inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);
                    
                        rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs * reference_a * integration_weight; 
                    }
                }
            }
        }
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Implicit
    ///@{

    void IgaEdgeCableElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        if (rDampingMatrix.size1() != mat_size)
            rDampingMatrix.resize(mat_size, mat_size, false);

        noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);

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

        // Rayleigh Damping Matrix: alpha*M + beta*K

        // 2.-Calculate StiffnessMatrix:
        if (beta > 0.0)
        {
            //MatrixType StiffnessMatrix = Matrix();
            Element::MatrixType StiffnessMatrix;

            if (StiffnessMatrix.size1() != mat_size)
                StiffnessMatrix.resize(mat_size, mat_size);
            noalias(StiffnessMatrix) = ZeroMatrix(mat_size, mat_size);

            // // //VectorType ResidualVector = Vector();
            // Element::VectorType ResidualVector;

            // if (ResidualVector.size() != mat_size)
            //     ResidualVector.resize(mat_size);
            // noalias(ResidualVector) = ZeroVector(mat_size);

            //this->CalculateAll(StiffnessMatrix, ResidualVector, rCurrentProcessInfo, true, false);
            this->CalculateInitialStiffnessMatrix(StiffnessMatrix, rCurrentProcessInfo);


            noalias(rDampingMatrix) += beta * StiffnessMatrix;
        }

        // 3.-Calculate MassMatrix:
        if (alpha > 0.0)
        {
            MatrixType MassMatrix = Matrix();
            this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);
            noalias(rDampingMatrix) += alpha * MassMatrix;
        }

        KRATOS_CATCH("")
    }

    void IgaEdgeCableElement::CalculateMassMatrix(
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

            double area = this->GetProperties().GetValue(CROSS_AREA);
            double density = this->GetProperties().GetValue(DENSITY);

            double mass = area * density * mReferenceArea[point_number] * integration_weight;

            if (rMassMatrix.size1() != mat_size)
                rMassMatrix.resize(mat_size, mat_size, false);
                
            rMassMatrix = ZeroMatrix(mat_size, mat_size);

            for (unsigned int r = 0; r<number_of_nodes; r++)
            {
                for (unsigned int s = 0; s<number_of_nodes; s++)
                {
                    rMassMatrix(3 * s, 3 * r) = r_N(point_number, s)*r_N(point_number, r)*mass;
                    rMassMatrix(3 * s + 1, 3 * r + 1) = rMassMatrix(3 * s, 3 * r);
                    rMassMatrix(3 * s + 2, 3 * r + 2) = rMassMatrix(3 * s, 3 * r);
                }
            }
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Postprocessing
    ///@{

    void IgaEdgeCableElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rValues.size() != r_integration_points.size())
        {
            rValues.resize(r_integration_points.size());
        }

        //get properties
        const double E = GetProperties()[YOUNG_MODULUS];
        const double A = GetProperties()[CROSS_AREA];
        const double prestress = GetProperties()[PRESTRESS_CAUCHY];

        if (rVariable == CABLE_FORCE)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) 
            {
                // get integration data
                const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(point_number);

                const array_1d<double, 3> actual_base_vector = GetActualBaseVector(r_DN_De, ConfigurationType::Current);
                const double reference_a = norm_2(mReferenceBaseVector[point_number]);
                const double actual_a = norm_2(actual_base_vector);

                const double actual_aa = actual_a * actual_a;
                const double reference_aa = reference_a * reference_a;

                // green-lagrange strain
                const double e11_membrane = 0.5 * (actual_aa - reference_aa);

                // normal forcereference_aa
                double principal_stress = prestress * A + e11_membrane * A * E / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]); 

                rValues[point_number] = principal_stress;
            }   
        }
    }

    ///@}
    ///@name Dynamic Functions
    ///@{

    void IgaEdgeCableElement::GetValuesVector(
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

    void IgaEdgeCableElement::GetFirstDerivativesVector(
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

    void IgaEdgeCableElement::GetSecondDerivativesVector(
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

    void IgaEdgeCableElement::EquationIdVector(
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

    void IgaEdgeCableElement::GetDofList(
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

    int IgaEdgeCableElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        // Verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }
        else
        {
            // Verify that the constitutive law has the correct dimension
            // KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
            //     << "THICKNESS not provided for element " << this->Id() << std::endl;

            // Check strain size
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
                << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
                << this->Id() << std::endl;
        }

        return 0;
    }

    ///@}

} // Namespace Kratos