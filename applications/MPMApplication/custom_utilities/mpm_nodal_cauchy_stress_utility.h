//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Laura Moreno Martínez

#pragma once

// System includes
#include <algorithm>
#include <limits>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/parallel_utilities.h"
#include "mpm_application_variables.h"

namespace Kratos
{

class MPMNodalCauchyStressUtility
{
public:
    static void CalculateNodalCauchyStress(
        ModelPart& rMaterialPointModelPart,
        ModelPart& rGridModelPart)
    {
        const ProcessInfo& r_current_process_info = rMaterialPointModelPart.GetProcessInfo();
        std::size_t max_node_id = 0;
        for (const auto& r_node : rGridModelPart.Nodes()) {
            max_node_id = std::max(max_node_id, r_node.Id());
        }
        Vector nodal_stress_weights(max_node_id + 1, 0.0);

        ResetNodalCauchyStress(rGridModelPart, r_current_process_info);
        AddMaterialPointStressContributions(rMaterialPointModelPart, r_current_process_info, nodal_stress_weights);
        NormalizeNodalCauchyStress(rGridModelPart, nodal_stress_weights);
    }

private:
    static std::size_t GetStressSize(const ProcessInfo& rCurrentProcessInfo)
    {
        const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC))
            ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
            : false;

        return (rCurrentProcessInfo[DOMAIN_SIZE] == 3) ? 6 : (is_axisymmetric ? 4 : 3);
    }

    static void ResetNodalCauchyStress(
        ModelPart& rGridModelPart,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t stress_size = GetStressSize(rCurrentProcessInfo);

        block_for_each(rGridModelPart.Nodes(), [&](Node& rNode) {
            auto& r_nodal_cauchy_stress_vector = rNode.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR);
            if (r_nodal_cauchy_stress_vector.size() != stress_size) {
                r_nodal_cauchy_stress_vector.resize(stress_size, false);
            }
            noalias(r_nodal_cauchy_stress_vector) = ZeroVector(stress_size);
        });
    }

    static void AddMaterialPointStressContributions(
        ModelPart& rMaterialPointModelPart,
        const ProcessInfo& rCurrentProcessInfo,
        Vector& rNodalStressWeights)
    {
        auto& r_elements = rMaterialPointModelPart.Elements();

        IndexPartition<std::size_t>(r_elements.size()).for_each([&](std::size_t i_elem) {
            auto it_elem = r_elements.begin() + i_elem;
            auto& r_geometry = it_elem->GetGeometry();
            const std::size_t number_of_nodes = r_geometry.PointsNumber();
            const std::size_t number_of_integration_points = r_geometry.IntegrationPointsNumber();

            std::vector<Vector> cauchy_stress_values;
            it_elem->CalculateOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, cauchy_stress_values, rCurrentProcessInfo);
            if (cauchy_stress_values.empty() || cauchy_stress_values[0].size() == 0) {
                return;
            }

            std::vector<double> mass_values;
            it_elem->CalculateOnIntegrationPoints(MP_MASS, mass_values, rCurrentProcessInfo);
            if (mass_values.empty()) {
                return;
            }

            const Vector& r_cauchy_stress_vector = cauchy_stress_values[0];
            const double material_point_mass = mass_values[0];
            const std::size_t stress_size = r_cauchy_stress_vector.size();

            for (std::size_t int_p = 0; int_p < number_of_integration_points; ++int_p) {
                const double weight = (number_of_integration_points > 1)
                    ? r_geometry.IntegrationPoints()[int_p].Weight()
                    : 1.0;

                for (std::size_t i = 0; i < number_of_nodes; ++i) {
                    const double shape_function_value = r_geometry.ShapeFunctionValue(int_p, i);
                    if (shape_function_value < 0.0) {
                        continue;
                    }

                    Vector nodal_cauchy_stress_vector = r_cauchy_stress_vector;
                    const double stress_weight = shape_function_value * material_point_mass * weight;
                    nodal_cauchy_stress_vector *= stress_weight;

                    r_geometry[i].SetLock();
                    auto& r_nodal_stress = r_geometry[i].FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR);
                    if (r_nodal_stress.size() != stress_size) {
                        r_nodal_stress.resize(stress_size, false);
                        noalias(r_nodal_stress) = ZeroVector(stress_size);
                    }
                    noalias(r_nodal_stress) += nodal_cauchy_stress_vector;
                    rNodalStressWeights[r_geometry[i].Id() - 1] += stress_weight;
                    r_geometry[i].UnSetLock();
                }
            }
        });
    }

    static void NormalizeNodalCauchyStress(
        ModelPart& rGridModelPart,
        const Vector& rNodalStressWeights)
    {
        block_for_each(rGridModelPart.Nodes(), [&](Node& rNode) {
            const double nodal_stress_weight = rNodalStressWeights[rNode.Id() - 1];
            if (nodal_stress_weight > std::numeric_limits<double>::epsilon()) {
                rNode.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR) /= nodal_stress_weight;
            }
        });
    }
};

} // namespace Kratos
