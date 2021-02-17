//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Michael Breitenberger
//                   Riccardo Rossi
//

// System includes

// External includes
#include "custom_conditions/lagrange_coupling_condition.h"

// Project includes

namespace Kratos
{

    void LagrangeCouplingCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        //const double penalty = GetProperties()[PENALTY_FACTOR];

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size = 6 * number_of_nodes_master + 3 * number_of_nodes_slave;

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
        Vector determinat_jacobian_vector_initial(integration_points.size());
        DeterminantOfJacobianInitial(r_geometry_master, determinat_jacobian_vector_initial);

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            Vector H;
            H.resize(number_of_nodes_master + number_of_nodes_slave);

            Vector HLambda;
            HLambda.resize(number_of_nodes_master);
            HLambda = ZeroVector(number_of_nodes_master);

            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                H[i] = N_master(point_number, i);
                HLambda[i] = N_master(point_number, i);
            }

            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                H[i + number_of_nodes_master] = -N_slave(point_number, i);
            }

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = determinat_jacobian_vector_initial[point_number];

            Matrix LHS = ZeroMatrix(mat_size, mat_size);

            for (IndexType i = 0; i < number_of_nodes_master; i++) //loop over Lagrange Multipliers
	        {
		        for (IndexType j = 0; j < number_of_nodes_master + number_of_nodes_slave; j++) // lopp over shape functions of displacements
		        {
                    double HH = H[j] * HLambda[i];
                        
                    const unsigned int ibase = i * 3 + 3 * (number_of_nodes_master + number_of_nodes_slave);
                    const unsigned int jbase = j * 3;

                    // Matrix in following shape:
			        // |0 H^T|
			        // |H 0  |

		            //lambda in X
		            LHS(ibase,     jbase)     = HH;
		            //lambda in Y
			        LHS(ibase + 1, jbase + 1) = HH;
			        //lambda in Z;
			        LHS(ibase + 2, jbase + 2) = HH;
		            //lambda in X
		            LHS(jbase,     ibase)     = HH;
		            //lambda in Y
			        LHS(jbase + 1, ibase + 1) = HH;
			        //lambda in Z;
		            LHS(jbase + 2, ibase + 2) = HH;

                    // // Depending on the solver if needed. 
		            // if (rLeftHandSideMatrix(3 * j, ibase) <= 0.000001)
			        // 	rLeftHandSideMatrix(i * 3, ibase) = 1e-6;
			        // if (rLeftHandSideMatrix(3 * j + 1, ibase + 1) <= 0.000001)
		        	// 	rLeftHandSideMatrix(i * 3 + 1, ibase + 1) = 1e-6;
		            // if (rLeftHandSideMatrix(3 * j + 2, ibase + 2) <= 0.000001)
		            // 	rLeftHandSideMatrix(i * 3 +2, ibase +2) = 1e-6;

                    // if (rLeftHandSideMatrix(ibase, 3 * j) <= 0.000001)
			        // 	rLeftHandSideMatrix(ibase, i * 3) = 1e-6;
			        // if (rLeftHandSideMatrix(ibase + 1, 3 * j + 1) <= 0.000001)
		            // 	rLeftHandSideMatrix(ibase + 1, i * 3 + 1) = 1e-6;
			        // if (rLeftHandSideMatrix(ibase + 2, 3 * j + 2) <= 0.000001)
			        // 	rLeftHandSideMatrix(ibase +2, i * 3 +2) = 1e-6;
                }
            }

            // Assembly (!: different than penalty)
            if (CalculateStiffnessMatrixFlag) {

                noalias(rLeftHandSideMatrix) += LHS
                    * integration_weight * determinat_jacobian;
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
                for (IndexType i = 0; i < number_of_nodes_master; i++)
                {
                    const array_1d<double, 3> disp = r_geometry_master[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                    IndexType index = 3 * (i + number_of_nodes_master + number_of_nodes_slave);
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }

                noalias(rRightHandSideVector) -= prod(LHS, u)
                     * integration_weight * determinat_jacobian; 
            }
        }

        KRATOS_CATCH("")
    }

    void LagrangeCouplingCondition::DeterminantOfJacobianInitial(
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

    void LagrangeCouplingCondition::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        // definition of problem size
        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size = 6 * number_of_nodes_master + 3 * number_of_nodes_slave;

        if (rDampingMatrix.size1() != mat_size)
            rDampingMatrix.resize(mat_size, mat_size, false);

        noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);

        // 1.-Get Damping Coeffitients (RAYLEIGH_BETA)

        double beta = 0.0;
        if (GetProperties().Has(RAYLEIGH_BETA))
            beta = GetProperties()[RAYLEIGH_BETA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
            beta = rCurrentProcessInfo[RAYLEIGH_BETA];

        //Rayleigh Damping Matrix: alpha*M + beta*K

        //2.-Calculate StiffnessMatrix:
        if (beta > 0.0)
        {
            //MatrixType StiffnessMatrix = Matrix();
            Condition::MatrixType StiffnessMatrix;

            if (StiffnessMatrix.size1() != mat_size)
                StiffnessMatrix.resize(mat_size, mat_size);
            noalias(StiffnessMatrix) = ZeroMatrix(mat_size, mat_size);

            //VectorType ResidualVector = Vector();
            Condition::VectorType ResidualVector;

            if (ResidualVector.size() != mat_size)
                ResidualVector.resize(mat_size);
            noalias(ResidualVector) = ZeroVector(mat_size);

            this->CalculateAll(StiffnessMatrix, ResidualVector, rCurrentProcessInfo, true, false);

            noalias(rDampingMatrix) += beta * StiffnessMatrix;

        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Dynamic Functions
    ///@{

    void LagrangeCouplingCondition::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_control_points_master = r_geometry_master.size();
        const SizeType number_of_control_points_slave = r_geometry_slave.size();
        const SizeType mat_size = 6* number_of_control_points_master + number_of_control_points_slave * 3;

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
            IndexType index = 3 * (i + number_of_control_points_master);

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }

        for (IndexType i = 0; i < number_of_control_points_master; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry_master[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, Step);
            IndexType index = 3 * (i + number_of_control_points_master + number_of_control_points_slave);

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void LagrangeCouplingCondition::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_control_points_master = r_geometry_master.size();
        const SizeType number_of_control_points_slave = r_geometry_slave.size();
        const SizeType mat_size = 6*number_of_control_points_master + number_of_control_points_slave*3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points_master; ++i)
        {
            const array_1d<double, 3 >& velocity = r_geometry_master[i].FastGetSolutionStepValue(VELOCITY, Step);
            IndexType index = i * 3;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }

        for (IndexType i = 0; i < number_of_control_points_slave; ++i)
        {
            const array_1d<double, 3 >& velocity = r_geometry_slave[i].FastGetSolutionStepValue(VELOCITY, Step);
            IndexType index = 3 * (i + number_of_control_points_master);

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }

        for (IndexType i = 0; i < number_of_control_points_master; ++i)
        {
            const array_1d<double, 3 >& velocity = r_geometry_master[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_VELOCITY, Step);
            IndexType index = 3 * (i + number_of_control_points_master + number_of_control_points_slave);

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void LagrangeCouplingCondition::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_control_points_master = r_geometry_master.size();
        const SizeType number_of_control_points_slave = r_geometry_slave.size();
        const SizeType mat_size = 6*number_of_control_points_master + number_of_control_points_slave*3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points_master; ++i)
        {
            const array_1d<double, 3 >& acceleration = r_geometry_master[i].FastGetSolutionStepValue(ACCELERATION, Step);
            IndexType index = i * 3;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }

        for (IndexType i = 0; i < number_of_control_points_slave; ++i)
        {
            const array_1d<double, 3 >& acceleration = r_geometry_slave[i].FastGetSolutionStepValue(ACCELERATION, Step);
            IndexType index = 3 * (i + number_of_control_points_master);

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }

        for (IndexType i = 0; i < number_of_control_points_master; ++i)
        {
            const array_1d<double, 3 >& acceleration = r_geometry_master[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_ACCELERATION, Step);
            IndexType index = 3 * (i + number_of_control_points_master + number_of_control_points_slave);

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void LagrangeCouplingCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        if (rResult.size() != 6 * number_of_nodes_master + 3* number_of_nodes_slave)
            rResult.resize(6 * number_of_nodes_master + 3* number_of_nodes_slave, false);

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

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = 3 * (i + number_of_nodes_master + number_of_nodes_slave);
            const auto& r_node = r_geometry_master[i];
            rResult[index]     = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
            rResult[index + 1] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void LagrangeCouplingCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * number_of_nodes_master + 3* number_of_nodes_slave);

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

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


