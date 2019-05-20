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

    void PenaltyCouplingCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        double augmention_factor = 0.0;
        if (Has(AUGMENTED_LAGRANGE_MULTIPLIER_FACTOR))
            augmention_factor = GetValue(AUGMENTED_LAGRANGE_MULTIPLIER_FACTOR);

        const auto& r_geometry_master = GetGeometry();
        const auto& r_geometry_slave = r_geometry_master.GetGeometryPart(1);

        // Size definitions
        const int number_of_nodes_master = r_geometry_master.size();
        const int number_of_nodes_slave = r_geometry_slave.size();

        const int mat_size = 3 * (2 * number_of_nodes_master + number_of_nodes_slave);

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
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints(integration_method);
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            Matrix N_master = r_geometry_master.ShapeFunctionsValues(integration_method);
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues(integration_method);

            //FOR DISPLACEMENTS
            Matrix K = ZeroMatrix(mat_size, mat_size);
            for (unsigned int i = 0; i < number_of_nodes_master; i++) //loop over Lagrange Multipliers
            {
                for (unsigned int j = 0; j < number_of_nodes_master; j++) // lopp over shape functions of displacements
                {
                    double NN = N_master[j] * N_master[i];

                    const unsigned int ibase = i * 3 + 3 * number_of_points;
                    const unsigned int jbase = j * 3;
                    const unsigned int ilambda = i * 3;

                    // Matrix in following shape:
                    // |0 H^T|
                    // |H 0  |
                    //lambda in X
                    rLeftHandSideMatrix(ibase, jbase)         = NN;// + Phi_r_Lambda(ilambda) *Phi_r(jbase);//Phi_r_Lambda((i - number_of_points)*3)*Phi_r(j*3); ShapeFunctionsN[i] * Phi_r(j * 3);// 
                    rLeftHandSideMatrix(ibase + 1, jbase + 1) = NN;// + Phi_r_Lambda(ilambda + 1) * Phi_r(jbase + 1);
                    rLeftHandSideMatrix(ibase + 2, jbase + 2) = NN;// + Phi_r_Lambda(ilambda + 2) * Phi_r(jbase + 2);

                    rLeftHandSideMatrix(jbase, ibase)         = NN;// + Phi_r_Lambda(ilambda) * Phi_r(jbase);
                    rLeftHandSideMatrix(jbase + 1, ibase + 1) = NN;// + Phi_r_Lambda(ilambda + 1) * Phi_r(jbase + 1);
                    rLeftHandSideMatrix(jbase + 2, ibase + 2) = NN;// + Phi_r_Lambda(ilambda + 2) * Phi_r(jbase + 2);

                    if (augmention_factor > 0.0)
                    {
                        //Depending on the solver if needed. 
                        if (rLeftHandSideMatrix(3 * j, ibase) <= 0.000001)
                            rLeftHandSideMatrix(i * 3, ibase) = augmention_factor;
                        if (rLeftHandSideMatrix(3 * j + 1, ibase + 1) <= 0.000001)
                            rLeftHandSideMatrix(i * 3 + 1, ibase + 1) = augmention_factor;
                        if (rLeftHandSideMatrix(3 * j + 2, ibase + 2) <= 0.000001)
                            rLeftHandSideMatrix(i * 3 +2, ibase +2) = augmention_factor;
                    }
                }
                for (unsigned int j = 0; j < number_of_nodes_slave; j++) // lopp over shape functions of displacements
                {
                    double NN = N_slave[j] * N_master[i];

                    const unsigned int ibase = i * 3 + 3 * number_of_points;
                    const unsigned int jbase = j * 3;
                    const unsigned int ilambda = i * 3;

                    // Matrix in following shape:
                    // |0 H^T|
                    // |H 0  |
                    //lambda in X
                    rLeftHandSideMatrix(ibase, jbase) = NN;// + Phi_r_Lambda(ilambda) *Phi_r(jbase);//Phi_r_Lambda((i - number_of_points)*3)*Phi_r(j*3); ShapeFunctionsN[i] * Phi_r(j * 3);// 
                    rLeftHandSideMatrix(ibase + 1, jbase + 1) = NN;// + Phi_r_Lambda(ilambda + 1) * Phi_r(jbase + 1);
                    rLeftHandSideMatrix(ibase + 2, jbase + 2) = NN;// + Phi_r_Lambda(ilambda + 2) * Phi_r(jbase + 2);

                    rLeftHandSideMatrix(jbase, ibase) = NN;// + Phi_r_Lambda(ilambda) * Phi_r(jbase);
                    rLeftHandSideMatrix(jbase + 1, ibase + 1) = NN;// + Phi_r_Lambda(ilambda + 1) * Phi_r(jbase + 1);
                    rLeftHandSideMatrix(jbase + 2, ibase + 2) = NN;// + Phi_r_Lambda(ilambda + 2) * Phi_r(jbase + 2);

                    if (augmention_factor > 0.0)
                    {
                        //Depending on the solver if needed. 
                        if (rLeftHandSideMatrix(3 * j, ibase) <= 0.000001)
                            rLeftHandSideMatrix(i * 3, ibase) = augmention_factor;
                        if (rLeftHandSideMatrix(3 * j + 1, ibase + 1) <= 0.000001)
                            rLeftHandSideMatrix(i * 3 + 1, ibase + 1) = augmention_factor;
                        if (rLeftHandSideMatrix(3 * j + 2, ibase + 2) <= 0.000001)
                            rLeftHandSideMatrix(i * 3 + 2, ibase + 2) = augmention_factor;
                    }
                }
            }

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = GetGeometry().DeterminatJacobian(point_number);

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += K * integration_weight * determinat_jacobian;
            }
            if (CalculateResidualVectorFlag) {

                Vector u(mat_size);
                for (unsigned int i = 0; i < number_of_nodes_master; i++)
                {
                    const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                    int index = 3 * i;
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }
                for (unsigned int i = 0; i < number_of_nodes_slave; i++)
                {
                    const array_1d<double, 3> disp = GetGeometry().GetSlaveGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                    int index = 3 * i + number_of_nodes_master;
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }
                for (unsigned int i = 0; i < number_of_nodes_master; i++)
                {
                    const array_1d<double, 3> lagrange_multiplier = GetGeometry()[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                    int index = 3 * i * 2 * number_of_nodes_master;
                    u[index]     = lagrange_multiplier[0];
                    u[index + 1] = lagrange_multiplier[1];
                    u[index + 2] = lagrange_multiplier[2];
                }

                noalias(rRightHandSideVector) -= prod(K, u)
                    * integration_weight * determinat_jacobian;
            }
        }

        KRATOS_CATCH("")
    }

    void PenaltyCouplingCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)

    {
        KRATOS_TRY;

        const int number_of_nodes_master = GetGeometry().size();
        const int number_of_nodes_slave = GetGeometry().GetSlaveGeometry().size();

        if (rResult.size() != 3 * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(3 * (number_of_nodes_master + number_of_nodes_slave), false);

        for (unsigned int i = 0; i < number_of_nodes_master; ++i) {
            const unsigned int index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }

        for (unsigned int i = 0; i < number_of_nodes_slave; ++i) {
            const unsigned int index = i * 3 + number_of_nodes_master;
            rResult[index]     = GetGeometry().GetSlaveGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry().GetSlaveGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry().GetSlaveGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }

        for (unsigned int i = 0; i < number_of_nodes_master; ++i) {
            const unsigned int index = i * 3 + 2 * number_of_nodes_master;
            rResult[index]     = GetGeometry().GetSlaveGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
            rResult[index + 1] = GetGeometry().GetSlaveGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
            rResult[index + 2] = GetGeometry().GetSlaveGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void PenaltyCouplingCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int number_of_nodes_master = GetGeometry().size();
        const int number_of_nodes_slave = GetGeometry().GetSlaveGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * (number_of_nodes_master + number_of_nodes_slave));

        for (unsigned int i = 0; i < number_of_nodes_master; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        for (unsigned int i = 0; i < number_of_nodes_slave; ++i) {
            rElementalDofList.push_back(GetGeometry().GetSlaveGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry().GetSlaveGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry().GetSlaveGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        for (unsigned int i = 0; i < number_of_nodes_master; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


