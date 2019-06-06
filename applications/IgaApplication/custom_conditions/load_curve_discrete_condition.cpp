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
#include "custom_conditions/load_curve_discrete_condition.h"


namespace Kratos
{

    void LoadCurveDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        const int number_of_control_points = NumberOfNodes();
        const int mat_size = NumberOfDofs();

        Vector f_loads = ZeroVector(mat_size);

        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        double integration_weight = this->GetValue(INTEGRATION_WEIGHT);


        if (Has(TANGENTS))
        {
            array_1d<double, 3> t2 = ZeroVector(3);
            array_1d<double, 2> tangents = GetValue(TANGENTS);
            IgaCurveOnSurfaceUtilities::CalculateTangent(GetGeometry(), DN_De, tangents, t2);
            integration_weight = integration_weight * norm_2(t2);
        }

        // Line loads
        if (this->Has(LINE_LOAD))
        {
            const array_1d<double, 3> line_load = this->GetValue(LINE_LOAD);

            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                f_loads[index]     += line_load[0] * integration_weight * N[i];
                f_loads[index + 1] += line_load[1] * integration_weight * N[i];
                f_loads[index + 2] += line_load[2] * integration_weight * N[i];
            }
        }

        // Pressure loads
        if (this->Has(PRESSURE))
        {
            double pressure = this->GetValue(PRESSURE);

            array_1d<double, 3> t1 = ZeroVector(3);
            array_1d<double, 2> tangents = GetValue(TANGENTS);
            IgaCurveOnSurfaceUtilities::CalculateNormal(GetGeometry(), DN_De, tangents, t1);

            t1 = t1 / norm_2(t1);

            KRATOS_WATCH(t1)

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                rRightHandSideVector[index] = -t1[0] * pressure * integration_weight * N[i];
                rRightHandSideVector[index + 1] = -t1[1] * pressure * integration_weight * N[i];
                rRightHandSideVector[index + 2] = -t1[2] * pressure * integration_weight * N[i];
            }
        }

        if (this->Has(MOMENT))
        {
            const array_1d<double, 3> moment = GetValue(MOMENT);

            array_1d<double, 2> Phi = ZeroVector(2);
            Vector Phi_r = ZeroVector(mat_size);
            Matrix Phi_rs = ZeroMatrix(mat_size, mat_size);

            array_1d<double, 3> g1, g2, g3;
            IgaSurfaceUtilities::CalculateBaseVectors(
                GetGeometry(),
                DN_De,
                g1, g2, g3);

            array_1d<double, 3> mG10, mG20, mG30;
            IgaSurfaceUtilities::CalculateInitialBaseVectors(
                GetGeometry(),
                DN_De,
                mG10, mG20, mG30);

            IgaCurveOnSurfaceUtilities::CalculateVariationRotation(
                GetGeometry(), DN_De, GetValue(TANGENTS),
                mG10, mG20, mG30,
                g1, g2, g3,
                Phi, Phi_r, Phi_rs);

            //KRATOS_WATCH(Phi_r)

            for (unsigned int n = 0; n < mat_size; n++)
            {
                //for (unsigned int m = 0; m < mat_size; m++)
                //{
                    //rLeftHandSideMatrix(n, m) += (Phi_r(n)*Phi_r(m);// +Phi(0)*Phi_rs(n, m));
                //}
                rRightHandSideVector(n) += Phi_r(n) * moment[0];
                //rRightHandSideVector(n) -= Phi_r(n) * moment[1];
            }
        }

        noalias(rRightHandSideVector) += f_loads;
    }

    void LoadCurveDiscreteCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
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
    }

    void LoadCurveDiscreteCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        const int number_of_control_points = NumberOfNodes();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
    };

} // Namespace Kratos