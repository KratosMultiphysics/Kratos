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
#include "custom_conditions/load_curve_discrete_condition.h"
#include "utilities/math_utils.h"
#include "includes/define.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"


namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    void LoadCurveDiscreteCondition::CalculateAll(
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
        double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        Vector g3 = ZeroVector(3);
        if (Has(TANGENTS))
        {
            CalculateNormalVector(g3, DN_De);

            Vector t2 = ZeroVector(3);
            CalculateBaseVector(t2, DN_De);
            integration_weight = integration_weight * norm_2(t2);
            //KRATOS_WATCH(integration_weight)
        }
        // Edge loads
        if (this->Has(LINE_LOAD))
        {
            array_1d<double, 3> line_load = this->GetValue(LINE_LOAD);

            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index]     = - line_load[0] * integration_weight * N[i];
                fLoads[index + 1] = - line_load[1] * integration_weight * N[i];
                fLoads[index + 2] = - line_load[2] * integration_weight * N[i];
            }
        }

        // Pressure loads
        if (this->Has(PRESSURE))
        {
            double pressure = this->GetValue(PRESSURE);

            array_1d<double, 3> direction = g3 / norm_2(g3);

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index]     = - direction[0] * pressure * integration_weight * N[i];
                fLoads[index + 1] = - direction[1] * pressure * integration_weight * N[i];
                fLoads[index + 2] = - direction[2] * pressure * integration_weight * N[i];
            }
        }

        noalias(rRightHandSideVector) -= fLoads;
    }
} // Namespace Kratos