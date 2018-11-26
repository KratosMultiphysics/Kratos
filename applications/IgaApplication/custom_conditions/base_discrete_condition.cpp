//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// System includes

// External includes

// Project includes
#include "custom_conditions/base_discrete_condition.h"


namespace Kratos
{
    void BaseDiscreteCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(DofsPerNode() != 3) << "Calling BaseDiscreteCondition::EquationIdVector. Derived Condition does not has 3 DoFs." << std::endl;

        const int number_of_control_points = NumberOfNodes();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const unsigned int index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    }

    void BaseDiscreteCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(DofsPerNode() != 3) << "Calling BaseDiscreteCondition::GetDofList. Derived Condition does not has 3 DoFs." << std::endl;

        const int number_of_control_points = NumberOfNodes();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    };

    void BaseDiscreteCondition::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        MatrixType left_hand_side_matrix = Matrix(0, 0);

        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(number_of_dofs);
        }

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    void BaseDiscreteCondition::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        if (rLeftHandSideMatrix.size1() != number_of_dofs
            || rLeftHandSideMatrix.size2() != number_of_dofs) {
            rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
        }

        VectorType right_hand_side_vector = Vector(0);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    void BaseDiscreteCondition::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        if (rLeftHandSideMatrix.size1() != number_of_dofs
            || rLeftHandSideMatrix.size2() != number_of_dofs) {
            rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
        }

        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(number_of_dofs);
        }

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    void BaseDiscreteCondition::Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == VELOCITY) {
            const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);

            array_1d<double, 3> velocity = ZeroVector(3);
            for (SizeType i = 0; i < N.size(); i++)
            {
                const NodeType & iNode = GetGeometry()[i];
                const array_1d<double, 3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, 0);

                velocity[0] += N[i] * vel[0];
                velocity[1] += N[i] * vel[1];
                velocity[2] += N[i] * vel[2];
            }
            rOutput = velocity;
        }
        else if (rVariable == DISPLACEMENT) {
            Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

            array_1d<double, 3> displacements = ZeroVector(3);
            for (SizeType i = 0; i < N.size(); i++)
            {
                const NodeType& iNode = GetGeometry()[i];
                const array_1d<double, 3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT);

                displacements[0] += N[i] * disp[0];
                displacements[1] += N[i] * disp[1];
                displacements[2] += N[i] * disp[2];
            }
            rOutput = displacements;
        }
        else
        {
            Condition::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    void BaseDiscreteCondition::Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == COORDINATES) 
		{
            const int number_of_control_points = NumberOfNodes();
            const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
            if (rOutput.size() != 3)
                rOutput.resize(3);
            rOutput = ZeroVector(3);
            for (SizeType i = 0; i < N.size(); i++)
            {
                const NodeType & iNode = GetGeometry()[i];
                const array_1d<double, 3>& coords = iNode.Coordinates();

                rOutput[0] += N[i] * coords[0];
                rOutput[1] += N[i] * coords[1];
                rOutput[2] += N[i] * coords[2];
            }
        }
        else if (rVariable == EXTERNAL_FORCES_VECTOR) 
        {
            const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
            if (rOutput.size() != 3)
                rOutput.resize(3);
            rOutput = ZeroVector(3);
            for (SizeType i = 0; i < N.size(); i++)
            {
                const NodeType & iNode = GetGeometry()[i];
                const array_1d<double, 3>& forces = iNode.GetValue(EXTERNAL_FORCES_VECTOR);

                rOutput[0] += N[i] * forces[0];
                rOutput[1] += N[i] * forces[1];
                rOutput[2] += N[i] * forces[2];
            }
        }
        else if (rVariable == SURFACE_NORMAL) {
            const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            Matrix Jacobian(3,2);
            CalculateJacobian(DN_De, Jacobian, 3, 2);

            Vector g1 = ZeroVector(3);
            Vector g2 = ZeroVector(3);

            if (rOutput.size() != 3)
                rOutput.resize(3);
            rOutput = ZeroVector(3);

            g1[0] = Jacobian(0, 0);
            g2[0] = Jacobian(0, 1);
            g1[1] = Jacobian(1, 0);
            g2[1] = Jacobian(1, 1);
            g1[2] = Jacobian(2, 0);
            g2[2] = Jacobian(2, 1);

            MathUtils<double>::CrossProduct(rOutput, g1, g2);
        }
        else {
            Condition::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    void BaseDiscreteCondition::SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);

        if (rVariable == EXTERNAL_FORCES_VECTOR) {
            for (SizeType i = 0; i < N.size(); i++)
            {
                NodeType & iNode = GetGeometry()[i];

                Vector external_variable = N[i] * rValues[0] + iNode.GetValue(rVariable);
                iNode.SetValue(rVariable, external_variable);
            }
        }
    }

    void BaseDiscreteCondition::GetValuesVector(
        Vector& rValues,
        int Step)
    {
        KRATOS_ERROR_IF(DofsPerNode() != 3) << "Calling BaseDiscreteCondition::GetValuesVector. Derived Condition does not has 3 DoFs." << std::endl;

        const int& number_of_control_points = NumberOfNodes();
        const int& mat_size = NumberOfDofs();

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const int index = i * 3;

            rValues[index]     = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void BaseDiscreteCondition::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) 
    {
        KRATOS_ERROR_IF(DofsPerNode() != 3) << "Calling BaseDiscreteCondition::GetFirstDerivativesVector. Derived Condition does not has 3 DoFs." << std::endl;

        const int& number_of_control_points = NumberOfNodes();
        const int& mat_size = NumberOfDofs();

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * 3;

            rValues[index]     = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void BaseDiscreteCondition::GetSecondDerivativesVector(
        Vector& rValues,
        int Step)
    {
        KRATOS_ERROR_IF(DofsPerNode() != 3) << "Calling BaseDiscreteCondition::GetSecondDerivativesVector. Derived Condition does not has 3 DoFs." << std::endl;

        const int& number_of_control_points = NumberOfNodes();
        const int mat_size = NumberOfDofs();

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * 3;

            rValues[index]     = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void BaseDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_ERROR << "You have called to the CalculateAll() from the base class BaseDiscreteCondition" << std::endl;
    }

    void BaseDiscreteCondition::CalculateJacobian(const Matrix& DN_De,
        Matrix& Jacobian,
        const int rWorkingSpaceDimension,
        const int rLocalSpaceDimension)
    {
        Jacobian.resize(rWorkingSpaceDimension, rLocalSpaceDimension);

        Jacobian.clear();
        for (unsigned int i = 0; i < DN_De.size1(); i++)
        {
            for (unsigned int k = 0; k<rWorkingSpaceDimension; k++)
            {
                for (unsigned int m = 0; m<rLocalSpaceDimension; m++)
                {
                    Jacobian(k, m) += (GetGeometry()[i]).Coordinates()[k] * DN_De(i, m);
                }
            }
        }
    }
} // Namespace Kratos