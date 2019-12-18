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
#include "custom_conditions/iga_base_condition.h"


namespace Kratos
{
    void IgaBaseCondition::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        MatrixType left_hand_side_matrix = Matrix(0, 0);

        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(number_of_dofs);
        }
        rRightHandSideVector = ZeroVector(number_of_dofs);

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    void IgaBaseCondition::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        if (rLeftHandSideMatrix.size1() != number_of_dofs
            || rLeftHandSideMatrix.size2() != number_of_dofs) {
            rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
        }
        rLeftHandSideMatrix = ZeroMatrix(number_of_dofs, number_of_dofs);

        VectorType right_hand_side_vector = Vector(0);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    void IgaBaseCondition::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        if (rLeftHandSideMatrix.size1() != number_of_dofs
            || rLeftHandSideMatrix.size2() != number_of_dofs) {
            rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
        }
        rLeftHandSideMatrix = ZeroMatrix(number_of_dofs, number_of_dofs);

        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(number_of_dofs);
        }
        rRightHandSideVector = ZeroVector(number_of_dofs);

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    void IgaBaseCondition::Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == COORDINATES)
        {
            IgaGeometryUtilities::CalculateCoordinates(
                GetGeometry(),
                GetValue(SHAPE_FUNCTION_VALUES),
                3,
                rOutput
            );
        }
        else if (rVariable == VELOCITY || rVariable == DISPLACEMENT) {
            IgaGeometryUtilities::CalculateSolutionStepValue(
                rVariable,
                GetGeometry(),
                GetValue(SHAPE_FUNCTION_VALUES),
                3,
                rOutput,
                rCurrentProcessInfo.GetSolutionStepIndex()
            );
        }
        else
        {
            Condition::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    void IgaBaseCondition::Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == EXTERNAL_FORCES_VECTOR) 
        {
            IgaGeometryUtilities::CalculateValue(
                rVariable,
                GetGeometry(),
                GetValue(SHAPE_FUNCTION_VALUES),
                3,
                rOutput
            );
        }
        else {
            Condition::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    void IgaBaseCondition::SetValueOnIntegrationPoints(
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

    void IgaBaseCondition::GetValuesVector(
        Vector& rValues,
        int Step)
    {
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

    void IgaBaseCondition::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) 
    {
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

    void IgaBaseCondition::GetSecondDerivativesVector(
        Vector& rValues,
        int Step)
    {
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

    void IgaBaseCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_ERROR << "You have called to the CalculateAll() from the base class IgaBaseCondition" << std::endl;
    }
} // Namespace Kratos