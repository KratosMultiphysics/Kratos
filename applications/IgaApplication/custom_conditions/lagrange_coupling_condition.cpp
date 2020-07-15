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
        ProcessInfo& rCurrentProcessInfo,
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
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            Vector N;
            N.resize(number_of_nodes_master + number_of_nodes_slave);

            Vector NLambda;
            NLambda.resize(number_of_nodes_master + number_of_nodes_slave);
            NLambda = ZeroVector(number_of_nodes_master + number_of_nodes_slave);

            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                N[i] = N_master(point_number, i);
                NLambda[i] = N_master(point_number, i);
            }

            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                N[i + number_of_nodes_master] = -N_slave(point_number, i);
            }

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = r_geometry_master.DeterminantOfJacobian(point_number);

            Matrix LHS = ZeroMatrix(mat_size, mat_size);

            for (IndexType i = 0; i < number_of_nodes_master + number_of_nodes_slave; i++) //loop over Lagrange Multipliers
	        {
		        for (IndexType j = 0; j < number_of_nodes_master + number_of_nodes_slave; j++) // lopp over shape functions of displacements
		        {
                    double NN = N[j] * NLambda[i];
                        
                    const unsigned int ibase = i * 3 + 3 * (number_of_nodes_master + number_of_nodes_slave);
                    const unsigned int jbase = j * 3;

                    // Matrix in following shape:
			        // |0 H^T|
			        // |H 0  |

		            //lambda in X
		            LHS(ibase,     jbase)     = NN;
		            //lambda in Y
			        LHS(ibase + 1, jbase + 1) = NN;
			        //lambda in Z;
			        LHS(ibase + 2, jbase + 2) = NN;
		            //lambda in X
		            LHS(jbase,     ibase)     = NN;
		            //lambda in Y
			        LHS(jbase + 1, ibase + 1) = NN;
			        //lambda in Z;
		            LHS(jbase + 2, ibase + 2) = NN;

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
                for (IndexType i = 0; i < number_of_nodes_slave; i++)
                {
                    const array_1d<double, 3> disp = r_geometry_slave[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                    IndexType index = 3 * (i + 2*number_of_nodes_master + number_of_nodes_slave);
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

    void LagrangeCouplingCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
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

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = 3 * (i + 2 * number_of_nodes_master + number_of_nodes_slave);
            const auto& r_node = r_geometry_slave[i];
            rResult[index]     = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
            rResult[index + 1] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void LagrangeCouplingCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
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

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


