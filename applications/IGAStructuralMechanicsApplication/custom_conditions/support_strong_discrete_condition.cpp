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
#include "custom_conditions/support_strong_discrete_condition.h"
#include "utilities/math_utils.h"

#include "custom_utilities/iga_flags.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{
    void SupportStrongDiscreteCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        const int number_of_control_points = GetGeometry().size();

        if (Has(DISPLACEMENT))
        {
            const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);
            const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);

            // DISPLACEMENTS
            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                if (N[i] > 1e-7)
                {
                    if (Is(IGAFlags::FIX_DISPLACEMENT_X))
                    {
                        GetGeometry()[i].Free(DISPLACEMENT_X);
                        GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X) = displacement[0];
                        GetGeometry()[i].Fix(DISPLACEMENT_X);
                    }
                    if (Is(IGAFlags::FIX_DISPLACEMENT_Y))
                    {
                        GetGeometry()[i].Free(DISPLACEMENT_Y);
                        GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y) = displacement[1];
                        GetGeometry()[i].Fix(DISPLACEMENT_Y);
                    }
                    if (Is(IGAFlags::FIX_DISPLACEMENT_Z))
                    {
                        GetGeometry()[i].Free(DISPLACEMENT_Z);
                        GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z) = displacement[2];
                        GetGeometry()[i].Fix(DISPLACEMENT_Z);
                    }
                }
            }

            if (Has(ROTATION))
            {
                const array_1d<double, 3>& rotation = this->GetValue(ROTATION);
                const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

                for (unsigned int i = 0; i < number_of_control_points; i++)
                {
                    if (DN_De(i, 0) > 1e-7 || DN_De(i, 1) > 1e-7)
                    {
                        if (Is(IGAFlags::FIX_ROTATION_X))
                        {
                            GetGeometry()[i].SetValue(DISPLACEMENT_X, displacement[0]);
                            GetGeometry()[i].Fix(DISPLACEMENT_X);
                        }
                        if (Is(IGAFlags::FIX_ROTATION_Y))
                        {
                            GetGeometry()[i].SetValue(DISPLACEMENT_Y, displacement[1]);
                            GetGeometry()[i].Fix(DISPLACEMENT_Y);
                        }
                        if (Is(IGAFlags::FIX_ROTATION_Z))
                        {
                            GetGeometry()[i].SetValue(DISPLACEMENT_Z, displacement[2]);
                            GetGeometry()[i].Fix(DISPLACEMENT_Z);
                        }
                    }
                }
            }
        }
    }

//************************************************************************************
//************************************************************************************
    void SupportStrongDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
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
    }
} // Namespace Kratos


