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
        const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
        Vector determinants_of_jacobian;
        r_geometry.DeterminantOfJacobian(determinants_of_jacobian);
        // Shape function values
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        if (rOutput.size() != r_integration_points.size())
            rOutput.resize(r_integration_points.size());


        if (rVariable == DISPLACEMENT_X || rVariable == DISPLACEMENT_Y || rVariable == DISPLACEMENT_Z) {
            const auto& r_N = r_geometry.ShapeFunctionsValues();
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                rOutput[point_number] = 0.0;
                for (IndexType i = 0; i < r_geometry.size(); ++i) {
                    rOutput[point_number] += r_geometry[i].FastGetSolutionStepValue(rVariable) * r_N(point_number, i);
                }
            }
        }
        else if (rVariable == INTEGRATION_COORDINATES_X || rVariable == INTEGRATION_COORDINATES_Y) {
            const auto& r_N = r_geometry.ShapeFunctionsValues();
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                if (rVariable == INTEGRATION_COORDINATES_X) {
                    rOutput[point_number] = r_integration_points[point_number][0];
                }
                else {
                    rOutput[point_number] = r_integration_points[point_number][1];
                }
            }
        }
        else {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            {
                rOutput[point_number] = 0.0;
                for (IndexType i = 0; i < nb_nodes; ++i)
                {
                    double output_solution_step_value = r_geometry[i].FastGetSolutionStepValue(rVariable);
                    rOutput[point_number] += r_N(point_number, i) * output_solution_step_value;
                }
                rOutput[point_number] *= r_integration_points[point_number].Weight() * determinants_of_jacobian[point_number];
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
        const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
        Vector determinants_of_jacobian;
        r_geometry.DeterminantOfJacobian(determinants_of_jacobian);
        // Shape function values
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        if (rOutput.size() != r_integration_points.size())
            rOutput.resize(r_integration_points.size());

        if (rVariable == PENALTY_REACTION)
        {
            SupportPenaltyCondition::CalculatePenaltyReaction(rOutput,
                GetProperties()[PENALTY_FACTOR], GetProperties()[DISPLACEMENT],
                r_integration_points, r_N, r_geometry);
        }
        else if (rVariable == PENALTY_REACTION_FORCE)
        {
            SupportPenaltyCondition::CalculatePenaltyReaction(rOutput,
                GetProperties()[PENALTY_FACTOR], GetProperties()[DISPLACEMENT],
                r_integration_points, r_N, r_geometry);
            for (IndexType i = 0; i < r_integration_points.size(); ++i) {
                rOutput[i] *= determinants_of_jacobian[i];
            }
        }
        else if (rVariable == REACTION) {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            {
                rOutput[point_number] = ZeroVector(3);
                for (IndexType i = 0; i < nb_nodes; ++i)
                {
                    array_1d<double, 3> output_solution_step_value = r_geometry[i].FastGetSolutionStepValue(rVariable);
                    rOutput[point_number] += r_N(point_number, i) * output_solution_step_value;
                }
                rOutput[point_number] *= r_integration_points[point_number].Weight();
            }
        }
        else {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            {
                rOutput[point_number] = ZeroVector(3);
                for (IndexType i = 0; i < nb_nodes; ++i)
                {
                    array_1d<double, 3> output_solution_step_value = r_geometry[i].FastGetSolutionStepValue(rVariable);
                    rOutput[point_number] += r_N(point_number, i) * output_solution_step_value;
                }
                rOutput[point_number] *= r_integration_points[point_number].Weight() * determinants_of_jacobian[point_number];
            }
        }
    }

} // Namespace Kratos
