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
#include "custom_conditions/penalty_coupling_condition.h"

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
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            //FOR DISPLACEMENTS
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

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = r_geometry_master.DeterminantOfJacobian(point_number);

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
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

                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
                    * integration_weight * determinat_jacobian * penalty;
            }
        }

        KRATOS_CATCH("")
    }

    void PenaltyCouplingCondition::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        // definition of problem size
        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size = 3 * (number_of_nodes_master + number_of_nodes_slave);

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

    void PenaltyCouplingCondition::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

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
            IndexType index = 3 * (i + number_of_control_points_master);

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void PenaltyCouplingCondition::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_control_points_master = r_geometry_master.size();
        const SizeType number_of_control_points_slave = r_geometry_slave.size();
        const SizeType mat_size = (number_of_control_points_master + number_of_control_points_slave) * 3;

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
    }

    void PenaltyCouplingCondition::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_control_points_master = r_geometry_master.size();
        const SizeType number_of_control_points_slave = r_geometry_slave.size();
        const SizeType mat_size = (number_of_control_points_master + number_of_control_points_slave) * 3;

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
    }

    void PenaltyCouplingCondition::EquationIdVector(
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

    void PenaltyCouplingCondition::GetDofList(
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


