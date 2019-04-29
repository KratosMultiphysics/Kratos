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
#include "custom_elements/iga_base_element.h"


namespace Kratos
{
    void IgaBaseElement::CalculateRightHandSide(
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

    void IgaBaseElement::CalculateLeftHandSide(
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

    void IgaBaseElement::CalculateLocalSystem(
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

    void IgaBaseElement::Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == FORCE)
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
        else if (rVariable == STRESS_RESULTANT_FORCE)
        {

        }
        else
        {
            Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    void IgaBaseElement::Calculate(
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
            Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    void IgaBaseElement::SetValueOnIntegrationPoints(
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

    void IgaBaseElement::GetValuesVector(
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

    void IgaBaseElement::GetFirstDerivativesVector(
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

    void IgaBaseElement::GetSecondDerivativesVector(
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

    void IgaBaseElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;
        const unsigned int number_of_control_points = GetGeometry().size();

        // Resizing as needed the LHS
        const unsigned int mat_size = number_of_control_points * 3;

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
            MatrixType StiffnessMatrix = Matrix();
            VectorType ResidualVector = Vector();
            this->CalculateAll(StiffnessMatrix, ResidualVector, rCurrentProcessInfo, true, false);
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

    void IgaBaseElement::AddExplicitContribution(
        const VectorType &rRHSVector, const Variable<VectorType> &rRHSVariable,
        Variable<array_1d<double, 3>> &rDestinationVariable,
        const ProcessInfo &rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const int number_of_control_points = NumberOfNodes();
        const int local_size = LocalSize();
        const int mat_size = NumberOfDofs();

        Vector damping_residual_contribution = ZeroVector(mat_size);
        // Calculate damping contribution to residual -->
        if ((this->GetProperties().Has(RAYLEIGH_ALPHA) ||
            this->GetProperties().Has(RAYLEIGH_BETA)) &&
            (rDestinationVariable != NODAL_INERTIA)) {
            Vector current_nodal_velocities = ZeroVector(mat_size);
            this->GetFirstDerivativesVector(current_nodal_velocities);
            Matrix damping_matrix = ZeroMatrix(mat_size, mat_size);
            ProcessInfo temp_process_information; // cant pass const ProcessInfo
            this->CalculateDampingMatrix(damping_matrix, temp_process_information);
            // current residual contribution due to damping
            noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
        }

        if (rRHSVariable == RESIDUAL_VECTOR &&
            rDestinationVariable == FORCE_RESIDUAL) {

            for (IndexType i = 0; i < number_of_control_points; ++i) {
                const SizeType index = local_size * i;

                array_1d<double, 3> &r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

                for (IndexType j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                    r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
                }
            }
        }

        if (rRHSVariable == RESIDUAL_VECTOR &&
            rDestinationVariable == MOMENT_RESIDUAL) {

            for (IndexType i = 0; i < number_of_control_points; ++i) {
                const SizeType index = (local_size * i) + msDimension;

                array_1d<double, 3> &r_moment_residual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

                for (IndexType j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                    r_moment_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
                }
            }
        }

        if (rDestinationVariable == NODAL_INERTIA) {
            Matrix element_mass_matrix = ZeroMatrix(mat_size, mat_size);
            ProcessInfo temp_info; // Dummy
            this->CalculateMassMatrix(element_mass_matrix, temp_info);

            for (IndexType i = 0; i < number_of_control_points; ++i) {
                double aux_nodal_mass = 0.0;
                array_1d<double, 3> aux_nodal_inertia(3, 0.0);

                const SizeType index = i * local_size;

                for (IndexType j = 0; j < mat_size; ++j) {
                    aux_nodal_mass += element_mass_matrix(index, j);
                    for (IndexType k = 0; k < msDimension; ++k) {
                        aux_nodal_inertia[k] += element_mass_matrix(index + msDimension + k, j);
                    }
                }

                #pragma omp atomic
                GetGeometry()[i].GetValue(NODAL_MASS) += aux_nodal_mass;

                array_1d<double, 3>& r_nodal_inertia = GetGeometry()[i].GetValue(NODAL_INERTIA);
                for (IndexType k = 0; k < msDimension; ++k) {
                    #pragma omp atomic
                    r_nodal_inertia[k] += std::abs(aux_nodal_inertia[k]);
                }
            }
        }

        KRATOS_CATCH("")
    }

    void IgaBaseElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_ERROR << "You have called to the CalculateAll() from the base class IgaBaseElement" << std::endl;
    }
} // Namespace Kratos