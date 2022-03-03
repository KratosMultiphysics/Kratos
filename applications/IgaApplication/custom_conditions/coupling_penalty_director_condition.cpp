//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "custom_conditions/coupling_penalty_director_condition.h"
#include "custom_elements/shell_5p_element.h"
#include "utilities/math_utils.h"

namespace Kratos
{

    BoundedVector<double, 3> EuclideanDistanceVector(const BoundedVector<double, 3>& a, const BoundedVector<double, 3>& b)
    {
        return a - b;
    }

    void CouplingPenaltyDirectorCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
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

        const SizeType mat_size = 2 * (number_of_nodes_master + number_of_nodes_slave);

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

        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        const bool integrate_conservative = GetProperties().Has(INTEGRATE_CONSERVATIVE)
            ? GetProperties()[INTEGRATE_CONSERVATIVE]
            : false;
        if (integrate_conservative) {
            DeterminantOfJacobianInitial(r_geometry_master, determinant_jacobian_vector);
        }
        else {
            r_geometry_master.DeterminantOfJacobian(determinant_jacobian_vector);
        }

        for (IndexType IntegrationPointIndex = 0; IntegrationPointIndex < integration_points.size(); IntegrationPointIndex++)
        {
            Vector N_master = row(r_geometry_master.ShapeFunctionsValues(),IntegrationPointIndex);
            Vector N_slave = row(r_geometry_slave.ShapeFunctionsValues(), IntegrationPointIndex);

            BoundedVector<double, 3> t_master=ZeroVector(3);
            for (IndexType i = 0; i < N_master.size(); ++i)
                t_master += N_master[i] * r_geometry_master[i].GetValue(DIRECTOR);

            BoundedVector<double, 3> t_slave = ZeroVector(3);
            for (IndexType i = 0; i < N_slave.size(); ++i)
                t_slave += N_slave[i] * r_geometry_slave[i].GetValue(DIRECTOR);

            const BoundedMatrix<double, 3, 3> Pd_master = (IdentityMatrix(3) - outer_prod(t_master, t_master)) / norm_2(t_master);
            const BoundedMatrix<double, 3, 3> Pd_slave = (IdentityMatrix(3) - outer_prod(t_slave, t_slave)) / norm_2(t_slave);
            const auto distanceVector = EuclideanDistanceVector(t_master, t_slave);
            const BoundedMatrix<double, 3, 3> S_master =  (inner_prod(distanceVector, t_master) * (3.0 * outer_prod(t_master, t_master) - IdentityMatrix(3)) - outer_prod(distanceVector,t_master) - trans(outer_prod(distanceVector, t_master)))/ norm_2_square(t_master);
            const BoundedMatrix<double, 3, 3> S_slave =  (inner_prod(distanceVector, t_slave) * (3.0 * outer_prod(t_slave, t_slave) - IdentityMatrix(3)) - outer_prod(distanceVector, t_slave) - trans(outer_prod(distanceVector, t_slave)))/ norm_2_square(t_slave);



            // Differential area
            const double penalty_integration = penalty * integration_points[IntegrationPointIndex].Weight() * determinant_jacobian_vector[IntegrationPointIndex];

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
                Matrix coupleStiffness = ZeroMatrix(mat_size, mat_size);
                for (IndexType i = 0; i < number_of_nodes_master; ++i)
                {
                    IndexType indexI = 2 * i;

                    for (IndexType j = 0; j < number_of_nodes_master; ++j)
                    {
                        IndexType indexJ = 2 * j;
                        const BoundedMatrix<double, 3, 2> BLA_masterI = r_geometry_master[i].GetValue(DIRECTORTANGENTSPACE);
                        const BoundedMatrix<double, 3, 2> BLA_masterJ = r_geometry_master[j].GetValue(DIRECTORTANGENTSPACE);
                        const BoundedMatrix<double, 3, 2> tmp3 = prod(Pd_master + S_master, BLA_masterJ); //useless temp due to nonworking ublas prod(prod())
                        const BoundedMatrix<double, 2, 2> tmp = N_master(j) * N_master(i) *  (prod(trans(BLA_masterI), tmp3) );

                        coupleStiffness(indexI,     indexJ   ) = tmp(0,0);
                        coupleStiffness(indexI + 1, indexJ   ) = tmp(1,0);
                        coupleStiffness(indexI,     indexJ +1) = tmp(0,1);
                        coupleStiffness(indexI + 1, indexJ +1) = tmp(1,1);
                    }
                }

                for (IndexType i = 0; i < number_of_nodes_slave; ++i)
                {
                    IndexType indexI = 2 * (i + number_of_nodes_master);

                    for (IndexType j = 0; j < number_of_nodes_slave; ++j)
                    {
                        IndexType indexJ = 2 * (j + number_of_nodes_master);
                        const BoundedMatrix<double, 3, 2> BLA_slaveI = r_geometry_slave[i].GetValue(DIRECTORTANGENTSPACE);
                        const BoundedMatrix<double, 3, 2> BLA_slaveJ = r_geometry_slave[j].GetValue(DIRECTORTANGENTSPACE);
                        const BoundedMatrix<double, 3, 2> tmp3 = prod(Pd_slave + S_slave, BLA_slaveJ); //useless temp due to nonworking ublas prod(prod())
                        const BoundedMatrix<double, 2, 2> tmp = N_slave(j) * N_slave(i) * (prod(trans(BLA_slaveI), tmp3));

                        coupleStiffness(indexI    , indexJ    ) = tmp(0, 0);
                        coupleStiffness(indexI + 1, indexJ    ) = tmp(1, 0);
                        coupleStiffness(indexI,     indexJ + 1) = tmp(0, 1);
                        coupleStiffness(indexI + 1, indexJ + 1) = tmp(1, 1);
                    }
                }


                for (IndexType i = 0; i < number_of_nodes_slave; ++i)
                {
                    IndexType indexI = 2 * (i + number_of_nodes_master);

                    for (IndexType j = 0; j < number_of_nodes_master; ++j)
                    {
                        IndexType indexJ = 2 * j;
                        const BoundedMatrix<double, 3, 2> BLA_slaveI = r_geometry_slave[i].GetValue(DIRECTORTANGENTSPACE);
                        const BoundedMatrix<double, 3, 2> BLA_masterJ = r_geometry_master[j].GetValue(DIRECTORTANGENTSPACE);
                        const BoundedMatrix<double, 3, 3> tmp2 = prod(Pd_slave, Pd_master); //useless temp due to nonworking ublas prod(prod())
                        const BoundedMatrix<double, 3, 2> tmp3 = prod(tmp2, BLA_masterJ); //useless temp due to nonworking ublas prod(prod())
                        const BoundedMatrix<double, 2, 2> tmp = N_master(j) * (-N_slave(i)) * prod(trans(BLA_slaveI), tmp3);

                        coupleStiffness(indexI,     indexJ    ) = tmp(0, 0);
                        coupleStiffness(indexI + 1, indexJ    ) = tmp(1, 0);
                        coupleStiffness(indexI,     indexJ + 1) = tmp(0, 1);
                        coupleStiffness(indexI + 1, indexJ + 1) = tmp(1, 1);
                    }
                }

                // this can be skipped later by transposing the thing before and inserting in the correct index for debuggin the follwing lines are currently used
                for (IndexType i = 0; i < number_of_nodes_master; ++i)
                {
                    IndexType indexI = 2 * i;

                    for (IndexType j = 0; j < number_of_nodes_slave; ++j)
                    {
                        IndexType indexJ = 2 * (j + number_of_nodes_master);
                        const BoundedMatrix<double, 3, 2> BLA_masterI = r_geometry_master[i].GetValue(DIRECTORTANGENTSPACE);
                        const BoundedMatrix<double, 3, 2> BLA_slaveJ = r_geometry_slave[j].GetValue(DIRECTORTANGENTSPACE);
                        const BoundedMatrix<double, 3, 3> tmp2 = prod(Pd_master, Pd_slave); //useless temp due to nonworking ublas prod(prod())
                        const BoundedMatrix<double, 3, 2> tmp3 = prod(tmp2, BLA_slaveJ); //useless temp due to nonworking ublas prod(prod())
                        const BoundedMatrix<double, 2, 2> tmp = -N_slave(j) * N_master(i) * (prod(trans(BLA_masterI), tmp3));

                        coupleStiffness(indexI,     indexJ    ) = tmp(0, 0);
                        coupleStiffness(indexI + 1, indexJ    ) = tmp(1, 0);
                        coupleStiffness(indexI,     indexJ + 1) = tmp(0, 1);
                        coupleStiffness(indexI + 1, indexJ + 1) = tmp(1, 1);
                    }
                }
                //std::cout << coupleStiffness << std::endl;
                noalias(rLeftHandSideMatrix) += coupleStiffness * penalty_integration;
            }
            if (CalculateResidualVectorFlag) {
                Vector coupleForce = ZeroVector(mat_size);
                for (IndexType i = 0; i < number_of_nodes_master; ++i)
                {
                    IndexType index = 2 * i;
                    const BoundedMatrix<double, 3, 2> BLA_master = r_geometry_master[i].GetValue(DIRECTORTANGENTSPACE);
                    const BoundedMatrix<double, 3, 2> tmp = prod(Pd_master, BLA_master);
                    const BoundedVector<double, 2> tmpVec = N_master(i) * prod(distanceVector, tmp);
                    coupleForce(index) = tmpVec(0);
                    coupleForce(index + 1) = tmpVec(1);
                }

                for (IndexType i = 0; i < number_of_nodes_slave; ++i)
                {
                    IndexType index = 2 * (i + number_of_nodes_master);
                    const BoundedMatrix<double, 3, 2> BLA_slave = r_geometry_slave[i].GetValue(DIRECTORTANGENTSPACE);
                    const BoundedMatrix<double, 3, 2> tmp = prod(Pd_slave, BLA_slave);
                    const BoundedVector<double, 2> tmpVec = -N_slave(i) * prod(distanceVector, tmp);
                    coupleForce(index) = tmpVec(0);
                    coupleForce(index + 1) = tmpVec(1);
                }

                noalias(rRightHandSideVector) -= coupleForce * penalty_integration;
            }
        }

        KRATOS_CATCH("")
    }

    void CouplingPenaltyDirectorCondition::DeterminantOfJacobianInitial(
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

    int CouplingPenaltyDirectorCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyCondition" << std::endl;
        return 0;
    }

    void CouplingPenaltyDirectorCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        if (rResult.size() != 2 * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(2 * (number_of_nodes_master + number_of_nodes_slave), false);

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry_master[i];
            rResult[index]     = r_node.GetDof(DIRECTORINC_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DIRECTORINC_Y).EquationId();
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = 2 * (i + number_of_nodes_master);
            const auto& r_node = r_geometry_slave[i];
            rResult[index]     = r_node.GetDof(DIRECTORINC_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DIRECTORINC_Y).EquationId();
        }

        KRATOS_CATCH("")
    }

    void CouplingPenaltyDirectorCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * (number_of_nodes_master + number_of_nodes_slave));

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DIRECTORINC_X));
            rElementalDofList.push_back(r_node.pGetDof(DIRECTORINC_Y));
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DIRECTORINC_X));
            rElementalDofList.push_back(r_node.pGetDof(DIRECTORINC_Y));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


