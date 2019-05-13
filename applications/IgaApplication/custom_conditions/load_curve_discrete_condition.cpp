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
            // KRATOS_WATCH("LoadCurveDiscreteCondition: CalculateAll")
            const array_1d<double, 3> line_load = this->GetValue(LINE_LOAD);

            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                f_loads[index]     += line_load[0] * integration_weight * N[i];
                f_loads[index + 1] += line_load[1] * integration_weight * N[i];
                f_loads[index + 2] += line_load[2] * integration_weight * N[i];
            }
        }

        // Follower loads (pressure)
		if (Has(PRESSURE))
		{
            array_1d<double, 3> g3 = ZeroVector(3);
            IgaSurfaceUtilties::CalculateBaseVector(GetGeometry(), DN_De, g3);
            const double d_area = norm_2(g3);
            array_1d<double, 3> direction = g3 / d_area;
            
            const double integration_weight_area = integration_weight * d_area;
            
            double pressure = GetValue(PRESSURE);

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                f_loads[index]     =  direction[0] * pressure * integration_weight_area * N[i];
                f_loads[index + 1] =  direction[1] * pressure * integration_weight_area * N[i];
                f_loads[index + 2] =  direction[2] * pressure * integration_weight_area * N[i];
            }
        }		

        noalias(rRightHandSideVector) += f_loads;
    }

    void LoadCurveDiscreteCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        // KRATOS_WATCH("LoadCurveDiscreteCondition: EquationIdVector")
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