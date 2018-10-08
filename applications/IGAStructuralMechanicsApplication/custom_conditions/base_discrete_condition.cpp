//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// System includes
#include "includes/define.h"
#include "utilities/math_utils.h"

// External includes

// Project includes
#include "custom_conditions/base_discrete_condition.h"

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "geometries/geometry.h"

namespace Kratos
{
    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const int number_of_control_points = GetGeometry().size();

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
    };

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const int number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    };

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo)
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = false;
        VectorType temp = Vector();

        CalculateAll(rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /***********************************************************************************/
    /// Calculate
    /***********************************************************************************/
    void BaseDiscreteCondition::Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rOutput.size() != 1)
            rOutput.resize(1);

        if (rVariable == VELOCITY) {
            const int& number_of_control_points = GetGeometry().size();
            const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);

            array_1d<double, 3> velocity = ZeroVector(3);
            for (SizeType i = 0; i < number_of_control_points; i++)
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
            const int& number_of_nodes = GetGeometry().size();
            Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

            array_1d<double, 3> displacements = ZeroVector(3);
            for (SizeType i = 0; i < number_of_nodes; i++)
            {
                const NodeType& iNode = GetGeometry()[i];
                const array_1d<double, 3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT);

                displacements[0] += N[i] * disp[0];
                displacements[1] += N[i] * disp[1];
                displacements[2] += N[i] * disp[2];
            }
            rOutput = displacements;
        }
        //else if (rVariable == REACTION) {
        //    const int& number_of_nodes = GetGeometry().size();
        //    Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

        //    array_1d<double, 3> reactions = ZeroVector(3);
        //    for (SizeType i = 0; i < number_of_nodes; i++)
        //    {
        //        const NodeType& iNode = GetGeometry()[i];
        //        const array_1d<double, 3>& reaction = iNode.FastGetSolutionStepValue(REACTION);

        //        reactions[0] += N[i] * reaction[0];
        //        reactions[1] += N[i] * reaction[1];
        //        reactions[2] += N[i] * reaction[2];
        //    }
        //    rOutput = reactions;
        //}
        else
        {
            Condition::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == COORDINATES) {
            const int& number_of_control_points = GetGeometry().size();
            const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
            Vector condition_coords = ZeroVector(3);
            for (SizeType i = 0; i < number_of_control_points; i++)
            {
                const NodeType & iNode = GetGeometry()[i];
                const array_1d<double, 3>& coords = iNode.Coordinates();

                condition_coords[0] += N[i] * coords[0];
                condition_coords[1] += N[i] * coords[1];
                condition_coords[2] += N[i] * coords[2];
            }
            rOutput = condition_coords;
        }
        else if (rVariable == EXTERNAL_FORCES_VECTOR) {
            const int& number_of_control_points = GetGeometry().size();
            const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
            Vector condition_coords = ZeroVector(3);
            for (SizeType i = 0; i < number_of_control_points; i++)
            {
                const NodeType & iNode = GetGeometry()[i];
                const array_1d<double, 3>& forces = iNode.GetValue(EXTERNAL_FORCES_VECTOR);

                KRATOS_WATCH(forces)

                condition_coords[0] += N[i] * forces[0];
                condition_coords[1] += N[i] * forces[1];
                condition_coords[2] += N[i] * forces[2];
            }
            rOutput = condition_coords;
        }
        else if (rVariable == SURFACE_NORMAL) {
            const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            Matrix Jacobian(3,2);
            CalculateJacobian(DN_De, Jacobian, 3, 2);

            Vector g1 = ZeroVector(3);
            Vector g2 = ZeroVector(3);
            Vector g3 = ZeroVector(3);

            g1[0] = Jacobian(0, 0);
            g2[0] = Jacobian(0, 1);
            g1[1] = Jacobian(1, 0);
            g2[1] = Jacobian(1, 1);
            g1[2] = Jacobian(2, 0);
            g2[2] = Jacobian(2, 1);

            MathUtils<double>::CrossProduct(g3, g1, g2);

            rOutput = g3;
        }
        else {
            Condition::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    /***********************************************************************************/
    /// SetValuesOnIntegrationPoints
    /***********************************************************************************/
    void BaseDiscreteCondition::SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const int number_of_control_points = GetGeometry().size();
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);

        if (rVariable == EXTERNAL_FORCES_VECTOR) {
            for (SizeType i = 0; i < number_of_control_points; i++)
            {
                NodeType & iNode = GetGeometry()[i];

                Vector external_variable = N[i] * rValues[0] + iNode.GetValue(rVariable);
                iNode.SetValue(rVariable, external_variable);
            }
        }
        else
        {
        }
    }


    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::GetValuesVector(
        Vector& rValues,
        int Step)
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

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

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::GetFirstDerivativesVector(
        Vector& rValues,
        int Step)
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

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

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteCondition::GetSecondDerivativesVector(
        Vector& rValues,
        int Step)
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

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

    //************************************************************************************/
    //************************************************************************************/
    void BaseDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_ERROR << "You have called to the CalculateAll() from the base class BaseDiscreteCondition" << std::endl;
    }

    //***********************************************************************************/
    //***********************************************************************************/
    void BaseDiscreteCondition::CalculateJacobian(const Matrix& DN_De,
        Matrix& Jacobian,
        const int rWorkingSpaceDimension,
        const int rLocalSpaceDimension)
    {
        const int number_of_control_points = GetGeometry().size();

        Jacobian.resize(rWorkingSpaceDimension, rLocalSpaceDimension);

        Jacobian.clear();
        for (unsigned int i = 0; i < number_of_control_points; i++)
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

/***********************************************************************************/
/***********************************************************************************/
} // Namespace Kratos