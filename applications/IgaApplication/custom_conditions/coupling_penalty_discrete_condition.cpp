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
#include "custom_conditions/coupling_penalty_discrete_condition.h"

// Project includes

namespace Kratos
{
    void CouplingPenaltyDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        const unsigned int number_of_control_points = NumberOfNodes();
        const unsigned int mat_size = NumberOfDofs();

        Vector N;
        GetShapeFunctions(N);

        const double Penalty = GetProperties()[PENALTY_FACTOR];
        const double Weighting = this->GetValue(INTEGRATION_WEIGHT);
        const Vector& localTrimTangents = this->GetValue(TANGENTS);
        const Matrix& ShapeFunctionDerivativesMaster = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        array_1d<double, 2> localTrimTangentsMaster;
        localTrimTangentsMaster[0] = localTrimTangents[0];
        localTrimTangentsMaster[1] = localTrimTangents[1];

        if (Is(IgaFlags::FIX_ROTATION_X))//(rot == 1)
        {
            Vector Phi_r = ZeroVector(mat_size);
            Vector Phi_r_Lambda = ZeroVector(mat_size);
            Matrix Phi_rs = ZeroMatrix(mat_size, mat_size);
            array_1d<double, 2> Diff_Phi;
            Diff_Phi.clear();

            CaculateRotationalShapeFunctions(Phi_r, Phi_r_Lambda, Phi_rs, Diff_Phi);

            for (unsigned int i = 0; i < mat_size; i++)
            {
                for (unsigned int j = 0; j < mat_size; j++)
                {
                    rLeftHandSideMatrix(i, j) = Phi_r(i)*Phi_r(j) + Diff_Phi(0)*Phi_rs(i, j);
                }
                rRightHandSideVector[i] = Diff_Phi(0)*Phi_r(i);
            }
        }



        //FOR DISPLACEMENTS
        Matrix Hcomplete = ZeroMatrix(3, mat_size);
        for (unsigned int i = 0; i < number_of_control_points; i++)
        {
            int index = 3 * i;
            if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                Hcomplete(0, index) = N[i];

            if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                Hcomplete(1, index + 1) = N[i];

            if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                Hcomplete(2, index + 2) = N[i];
        }

        Vector TDisplacements(mat_size);
        for (unsigned int i = 0; i < number_of_control_points; i++)
        {
            const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            int index = 3 * i;
            TDisplacements[index] = disp[0];
            TDisplacements[index + 1] = disp[1];
            TDisplacements[index + 2] = disp[2];
        }

        double JGeometrictoParameter;
        MappingGeometricToParameterMasterElement(ShapeFunctionDerivativesMaster, localTrimTangentsMaster, JGeometrictoParameter);

        noalias(rLeftHandSideMatrix) += prod(trans(Hcomplete), Hcomplete);
        noalias(rRightHandSideVector) += prod(prod(trans(Hcomplete), Hcomplete), TDisplacements);

        //Mapping:
        rLeftHandSideMatrix  *= (Weighting * JGeometrictoParameter * Penalty);
        rRightHandSideVector *= (Weighting * JGeometrictoParameter * Penalty);

        KRATOS_CATCH("")
    }

    void CouplingPenaltyDiscreteCondition::GetShapeFunctions(Vector& rShapeFunctions)
    {
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Vector& NSlave = this->GetValue(SHAPE_FUNCTION_VALUES_SLAVE);

        if (rShapeFunctions.size() != N.size() + NSlave.size())
            rShapeFunctions.resize(N.size() + NSlave.size());
        rShapeFunctions = ZeroVector(N.size() + NSlave.size());

        for (unsigned int i = 0; i < N.size(); i++)
        {
            rShapeFunctions[i] = N[i];
        }
        for (unsigned int i = 0; i < NSlave.size(); i++)
        {
            rShapeFunctions[i + N.size()] = -NSlave[i];
        }
    }

    void CouplingPenaltyDiscreteCondition::EquationIdVector(
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

    void CouplingPenaltyDiscreteCondition::GetDofList(
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


