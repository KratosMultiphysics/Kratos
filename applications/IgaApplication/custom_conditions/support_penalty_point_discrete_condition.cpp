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
#include "custom_conditions/support_penalty_point_discrete_condition.h"


namespace Kratos
{
    void SupportPenaltyPointDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const int number_of_control_points = NumberOfNodes();
        const int mat_size = NumberOfDofs();

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
                    Stiffness(0, 3 * i) = N[i];
                if (Is(IGAFlags::FIX_DISPLACEMENT_Y))
                    Stiffness(1, 3 * i + 1) = N[i];
                if (Is(IGAFlags::FIX_DISPLACEMENT_Z))
                    Stiffness(2, 3 * i + 2) = N[i];
            }
            const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);

            Vector TDisplacements(mat_size);
            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                int index = 3 * i;
                TDisplacements[index]     = (disp[0] - displacement[0]);
                TDisplacements[index + 1] = (disp[1] - displacement[1]);
                TDisplacements[index + 2] = (disp[2] - displacement[2]);
            }

            noalias(rLeftHandSideMatrix) += prod(trans(Stiffness), Stiffness);
            noalias(rRightHandSideVector) -= prod(prod(trans(Stiffness), Stiffness), TDisplacements);
        }

        rLeftHandSideMatrix  *= Penalty;
        rRightHandSideVector *= Penalty;

        KRATOS_CATCH("")
    }

    int SupportPenaltyPointDiscreteCondition::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        KRATOS_ERROR_IF(SHAPE_FUNCTION_VALUES.Key() == 0) << "SHAPE_FUNCTION_VALUES has Key zero! check if the application is correctly registered" << std::endl;
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR)) << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyPointDiscreteCondition" << std::endl;
        return 0;
        KRATOS_CATCH("");
    }
} // Namespace Kratos


