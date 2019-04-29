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
#include "custom_conditions/support_penalty_curve_discrete_condition.h"


namespace Kratos
{
    void SupportPenaltyCurveDiscreteCondition::Initialize()
    {
        IgaSurfaceUtilities::CalculateBaseVectors(
            GetGeometry(), GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES),
            mG10, mG20, mG30);
    }

    void SupportPenaltyCurveDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        //Read in Data
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);

        double penalty = GetProperties()[PENALTY_FACTOR];

        //ROTATIONS
        if (Has(ROTATION))
        {
            const array_1d<double, 3>& rotation = this->GetValue(ROTATION);

            array_1d<double, 2> Phi = ZeroVector(2);
            Vector Phi_r = ZeroVector(mat_size);
            Matrix Phi_rs = ZeroMatrix(mat_size, mat_size);

            array_1d<double, 3> g1, g2, g3;
            IgaSurfaceUtilities::CalculateBaseVectors(
                GetGeometry(), DN_De,
                g1, g2, g3);

            IgaCurveOnSurfaceUtilities::CalculateVariationRotation(
                GetGeometry(), DN_De, GetValue(TANGENTS),
                mG10, mG20, mG30,
                g1, g2, g3,
                Phi, Phi_r, Phi_rs);

            for (unsigned int n = 0; n < mat_size; n++)
            {
                if (CalculateStiffnessMatrixFlag)
                {
                    for (unsigned int m = 0; m < mat_size; m++)
                    {
                        rLeftHandSideMatrix(n, m) += (Phi_r(n)*Phi_r(m) + Phi(0)*Phi_rs(n, m));
                    }
                }
                if (CalculateResidualVectorFlag)
                    rRightHandSideVector(n) -= Phi(0)*Phi_r(n);
            }
        }

        if (Has(DISPLACEMENT))
        {
            // DISPLACEMENTS
            Matrix Stiffness = ZeroMatrix(3, mat_size);
            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                {
                    Stiffness(0, 3 * i) = N[i];
                }
                if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                {
                    Stiffness(1, 3 * i + 1) = N[i];
                }
                if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
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
            if (CalculateStiffnessMatrixFlag)
                noalias(rLeftHandSideMatrix) += prod(trans(Stiffness), Stiffness);
            if (CalculateResidualVectorFlag)
                noalias(rRightHandSideVector) -= prod(prod(trans(Stiffness), Stiffness), TDisplacements);
        }

        array_1d<double, 3> t2 = ZeroVector(3);
        IgaCurveOnSurfaceUtilities::CalculateTangent(
            GetGeometry(),
            DN_De,
            GetValue(TANGENTS),
            t2);

        double weighting = integration_weight * norm_2(t2);

        //MAPPING
        if (CalculateStiffnessMatrixFlag)
            rLeftHandSideMatrix *= weighting * penalty;
        if (CalculateResidualVectorFlag)
            rRightHandSideVector *= weighting * penalty;

        KRATOS_CATCH("")
    }

    void SupportPenaltyCurveDiscreteCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int number_of_control_points = NumberOfNodes();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const unsigned int index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    }

    void SupportPenaltyCurveDiscreteCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int number_of_control_points = NumberOfNodes();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos