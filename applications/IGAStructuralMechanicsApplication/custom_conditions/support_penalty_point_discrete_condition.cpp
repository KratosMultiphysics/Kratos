//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/support_penalty_point_discrete_condition.h"
#include "utilities/math_utils.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

#include "custom_utilities/iga_flags.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
    void SupportPenaltyPointDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        //resizing the system in case it does not have the right size
        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS

        //resizing as needed the RHS
        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size, false);
        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

        //Read in Data
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);

        double Penalty = GetProperties()[PENALTY_FACTOR];

        if (Has(DISPLACEMENT))
        {
            // DISPLACEMENTS
            Matrix Stiffness = ZeroMatrix(3, mat_size);
            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                if (Is(IGAFlags::FIX_DISPLACEMENT_X))
                {
                    Stiffness(0, 3 * i) = N[i];
                }
                if (Is(IGAFlags::FIX_DISPLACEMENT_Y))
                {
                    Stiffness(1, 3 * i + 1) = N[i];
                }
                if (Is(IGAFlags::FIX_DISPLACEMENT_Z))
                {
                    Stiffness(2, 3 * i + 2) = N[i];
                }
            }
            const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);

            Vector TDisplacements(mat_size);
            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                int index = 3 * i;
                TDisplacements[index] = (disp[0] - displacement[0]);
                TDisplacements[index + 1] = (disp[1] - displacement[1]);
                TDisplacements[index + 2] = (disp[2] - displacement[2]);
            }

            noalias(rLeftHandSideMatrix) += prod(trans(Stiffness), Stiffness);
            noalias(rRightHandSideVector) -= prod(prod(trans(Stiffness), Stiffness), TDisplacements);
        }

        //MAPPING
        rLeftHandSideMatrix  *= Penalty;
        rRightHandSideVector *= Penalty;

        KRATOS_CATCH("")
    }
} // Namespace Kratos


