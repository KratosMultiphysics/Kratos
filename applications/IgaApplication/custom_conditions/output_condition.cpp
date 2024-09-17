//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes
#include "custom_conditions/output_condition.h"

// Project includes

namespace Kratos
{
    void OutputCondition::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType nb_nodes = r_geometry.size();

        // Integration Points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        // Shape function values
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        if (rOutput.size() != integration_points.size())
            rOutput.resize(integration_points.size());

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            rOutput[point_number] = 0.0;
            for (IndexType i = 0; i < nb_nodes; ++i)
            {
                double output_solution_step_value = r_geometry[i].FastGetSolutionStepValue(rVariable);
                rOutput[point_number] += r_N(point_number, i) * output_solution_step_value;
            }
        }
    }

    void OutputCondition::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const auto& r_geometry = GetGeometry();
        const SizeType nb_nodes = r_geometry.size();

        // Integration Points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        // Shape function values
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        if (rOutput.size() != integration_points.size())
            rOutput.resize(integration_points.size());

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            rOutput[point_number] = ZeroVector(3);
            for (IndexType i = 0; i < nb_nodes; ++i)
            {
                array_1d<double, 3> output_solution_step_value = r_geometry[i].FastGetSolutionStepValue(rVariable);
                rOutput[point_number] += r_N(point_number, i) * output_solution_step_value;
            }
        }
    }

    void OutputCondition::AddExplicitContribution(
        const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;

        KRATOS_CATCH( "" )
    }

} // Namespace Kratos
