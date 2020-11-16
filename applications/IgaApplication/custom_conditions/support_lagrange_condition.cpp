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
//                   Tobias Teschemacher
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_lagrange_condition.h"


namespace Kratos
{
    void SupportLagrangeCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = 2* r_geometry.WorkingSpaceDimension() * number_of_nodes;

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


        if (Has(DISPLACEMENT))
        {
            // Integration
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

            // initial determinant of jacobian 
            Vector determinat_jacobian_vector_initial(integration_points.size());
            DeterminantOfJacobianInitial(r_geometry, determinat_jacobian_vector_initial);

            for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
            {
                const Matrix& N = r_geometry.ShapeFunctionsValues();

                Vector H;
                H.resize(number_of_nodes);

                Vector HLambda;
                HLambda.resize(number_of_nodes);
                HLambda = ZeroVector(number_of_nodes);

                for (IndexType i = 0; i < number_of_nodes; i++)
                {
                    H[i] = N(point_number, i);
                    HLambda[i] = N(point_number, i);
                }

                // Differential area
                const double integration_weight = integration_points[point_number].Weight();
                const double determinat_jacobian = determinat_jacobian_vector_initial[point_number];

                Matrix LHS = ZeroMatrix(mat_size, mat_size);

                for (IndexType i = 0; i < number_of_nodes; i++) //loop over Lagrange Multipliers
                {
                    for (IndexType j = 0; j < number_of_nodes; j++) // lopp over shape functions of displacements
                    {
                        double HH = H[j] * HLambda[i];
                            
                        const unsigned int ibase = i * 3 + 3 * (number_of_nodes);
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

                if (CalculateStiffnessMatrixFlag) {

                    noalias(rLeftHandSideMatrix) += LHS* integration_weight * determinat_jacobian;
                }

                if (CalculateResidualVectorFlag) {

                    const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);

                    Vector u(mat_size);
                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                        IndexType index = 3 * i;
                        u[index]     = (disp[0] - displacement[0]);
                        u[index + 1] = (disp[1] - displacement[1]);
                        u[index + 2] = (disp[2] - displacement[2]);
                    }
                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                        IndexType index = 3 * (i + number_of_nodes);
                        u[index]     = (disp[0] - displacement[0]);
                        u[index + 1] = (disp[1] - displacement[1]);
                        u[index + 2] = (disp[2] - displacement[2]);
                    }

                    noalias(rRightHandSideVector) -= prod(LHS, u) * integration_weight * determinat_jacobian; 
                }
            }
        }
        KRATOS_CATCH("")
    }

    void SupportLagrangeCondition::DeterminantOfJacobianInitial(
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
            GetGeometry().Calculate(LOCAL_TANGENT, local_tangent);

            array_1d<double, 3> a_1 = column(J, 0);
            array_1d<double, 3> a_2 = column(J, 1);

            rDeterminantOfJacobian[pnt] = norm_2(a_1 * local_tangent[0] + a_2 * local_tangent[1]);
        }
    }


    int SupportLagrangeCondition::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        //KRATOS_WATCH(GetGeometry())
        return 0;
        KRATOS_CATCH("");
    }   

    void SupportLagrangeCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 6 * number_of_nodes)
            rResult.resize(6 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = 3 * (i + number_of_nodes);
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
            rResult[index + 1] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
        }
    }

    void SupportLagrangeCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
        }
    };

} // Namespace Kratos
