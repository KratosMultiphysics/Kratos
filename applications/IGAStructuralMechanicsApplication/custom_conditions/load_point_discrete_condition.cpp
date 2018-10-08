//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// System includes


// External includes


// Project includes
#include "custom_conditions/load_point_discrete_condition.h"
#include "utilities/math_utils.h"
#include "includes/define.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"


namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    void LoadPointDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        Vector fLoads = ZeroVector(mat_size);

        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS

        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size, false);
        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

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
                fLoads[index] = PointLoad[0] * N[i];
                fLoads[index + 1] = PointLoad[1] * N[i];
                fLoads[index + 2] = PointLoad[2] * N[i];
            }

            //FinalizeSolutionStep(rCurrentProcessInfo);
        }

        // Pressure loads
        if (this->Has(PRESSURE))
        {
            const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
            const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

            Vector g3 = ZeroVector(3);
            CalculateBaseVector(g3, DN_De);

            const double dArea = norm_2(g3);
            const double integration_weight_area = integration_weight * dArea;

            double pressure = this->GetValue(PRESSURE);

            array_1d<double, 3> direction = g3 / dArea;

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index] = -direction[0] * pressure * integration_weight_area * N[i];
                fLoads[index + 1] = -direction[1] * pressure * integration_weight_area * N[i];
                fLoads[index + 2] = -direction[2] * pressure * integration_weight_area * N[i];
            }
        }

        noalias(rRightHandSideVector) -= fLoads;
    }

    void LoadPointDiscreteCondition::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        //std::cout << "arrived here" << std::endl;
        // Interface Condition
        if (this->Has(EXTERNAL_FORCES_VECTOR))
        {
            std::vector<Vector> PointLoads(1);
            PointLoads[0] = this->GetValue(EXTERNAL_FORCES_VECTOR);

            KRATOS_WATCH(PointLoads[0])

            SetValueOnIntegrationPoints(EXTERNAL_FORCES_VECTOR, PointLoads, rCurrentProcessInfo);
        }
    }
} // Namespace Kratos