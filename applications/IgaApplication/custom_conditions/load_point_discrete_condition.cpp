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
#include "load_point_discrete_condition.h"



namespace Kratos
{
    void LoadPointDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        const int number_of_control_points = NumberOfNodes();
        const int mat_size = NumberOfDofs();
        Vector fLoads = ZeroVector(mat_size);

        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);

        // Point loads
        if (this->Has(POINT_LOAD))
        {
            const array_1d<double, 3 > PointLoad = this->GetValue(POINT_LOAD);

            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index]     = -PointLoad[0] * N[i];
                fLoads[index + 1] = -PointLoad[1] * N[i];
                fLoads[index + 2] = -PointLoad[2] * N[i];
            }
        }

        // Point loads
        if (this->Has(EXTERNAL_FORCES_VECTOR))
        {
            const array_1d<double, 3 > PointLoad = this->GetValue(EXTERNAL_FORCES_VECTOR);

            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index]     = -PointLoad[0] * N[i];
                fLoads[index + 1] = -PointLoad[1] * N[i];
                fLoads[index + 2] = -PointLoad[2] * N[i];
            }
        }

        // Pressure loads
        if (this->Has(PRESSURE))
        {
            double pressure = this->GetValue(PRESSURE);

            const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);

            Vector direction = ZeroVector(3);
            CalculateBaseVectorOnSurface(direction);

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index] = -direction[0] * pressure * integration_weight * N[i];
                fLoads[index + 1] = -direction[1] * pressure * integration_weight * N[i];
                fLoads[index + 2] = -direction[2] * pressure * integration_weight * N[i];
            }
        }

        noalias(rRightHandSideVector) -= fLoads;
    }
} // Namespace Kratos