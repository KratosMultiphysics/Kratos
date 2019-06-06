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
#include "custom_conditions/load_surface_discrete_condition.h"


namespace Kratos
{
    void LoadSurfaceDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        const int number_of_control_points = NumberOfNodes();

        const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        array_1d<double, 3> g3 = ZeroVector(3);
        IgaSurfaceUtilities::CalculateBaseVector(
            GetGeometry(),
            DN_De,
            g3);

        const double d_area = norm_2(g3);
        const double integration_weight_area = integration_weight * d_area;

        // Surface loads
        if (this->Has(SURFACE_LOAD))
        {
            Vector surface_load = this->GetValue(SURFACE_LOAD);

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                rRightHandSideVector[index]     = - surface_load[0] * integration_weight_area * N[i];
                rRightHandSideVector[index + 1] = - surface_load[1] * integration_weight_area * N[i];
                rRightHandSideVector[index + 2] = - surface_load[2] * integration_weight_area * N[i];
            }
        }

        // Pressure loads
        if (this->Has(PRESSURE))
        {
            double pressure = this->GetValue(PRESSURE);

            array_1d<double, 3> direction = g3 / d_area;

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                rRightHandSideVector[index]     = - direction[0] * pressure * integration_weight_area * N[i];
                rRightHandSideVector[index + 1] = - direction[1] * pressure * integration_weight_area * N[i];
                rRightHandSideVector[index + 2] = - direction[2] * pressure * integration_weight_area * N[i];
            }
        }
    }

    void LoadSurfaceDiscreteCondition::EquationIdVector(
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

    void LoadSurfaceDiscreteCondition::GetDofList(
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