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
#include <limits>
#include <unordered_map>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/atomic_utilities.h"
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

        std::unordered_map<ModelPart::IndexType, std::size_t> nodal_stress_weight_indices;
        nodal_stress_weight_indices.reserve(rGridModelPart.NumberOfNodes());
        std::size_t node_counter = 0;
        for (const auto& r_node : rGridModelPart.Nodes()) {
            nodal_stress_weight_indices[r_node.Id()] = node_counter++;
        }

        Vector nodal_stress_weights(rGridModelPart.NumberOfNodes(), 0.0);

        ResetNodalCauchyStress(rGridModelPart, r_current_process_info);
        AddMaterialPointStressContributions(
            rMaterialPointModelPart,
            r_current_process_info,
            nodal_stress_weight_indices,
            nodal_stress_weights);
        NormalizeNodalCauchyStress(
            rGridModelPart,
            nodal_stress_weight_indices,
            nodal_stress_weights);
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
        const std::unordered_map<ModelPart::IndexType, std::size_t>& rNodalStressWeightIndices,
        Vector& rNodalStressWeights)
    {
        auto& r_elements = rMaterialPointModelPart.Elements();

        IndexPartition<std::size_t>(r_elements.size()).for_each([&](std::size_t i_elem) {
            auto it_elem = r_elements.begin() + i_elem;
            auto& r_geometry = it_elem->GetGeometry();
            const std::size_t number_of_nodes = r_geometry.PointsNumber();
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();

            std::vector<Vector> cauchy_stress_values;
            it_elem->CalculateOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, cauchy_stress_values, rCurrentProcessInfo);

            std::vector<double> mass_values;
            it_elem->CalculateOnIntegrationPoints(MP_MASS, mass_values, rCurrentProcessInfo);

            for (std::size_t i = 0; i < number_of_nodes; ++i) {

                Vector nodal_cauchy_stress_vector = cauchy_stress_values[0];
                const double stress_weight = r_N(0,i) * mass_values[0];
                nodal_cauchy_stress_vector *= stress_weight;

                AtomicAddVector(r_geometry[i].FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR), nodal_cauchy_stress_vector);
                const auto nodal_stress_weight_index = rNodalStressWeightIndices.at(r_geometry[i].Id());
                AtomicAdd(rNodalStressWeights[nodal_stress_weight_index], stress_weight);
            }
        });
    }

    static void NormalizeNodalCauchyStress(
        ModelPart& rGridModelPart,
        const std::unordered_map<ModelPart::IndexType, std::size_t>& rNodalStressWeightIndices,
        const Vector& rNodalStressWeights)
    {
        block_for_each(rGridModelPart.Nodes(), [&](Node& rNode) {
            const auto nodal_stress_weight_index = rNodalStressWeightIndices.at(rNode.Id());
            const double nodal_stress_weight = rNodalStressWeights[nodal_stress_weight_index];
            if (nodal_stress_weight > std::numeric_limits<double>::epsilon()) {
                rNode.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR) /= nodal_stress_weight;
            }
        });
    }
};

} // namespace Kratos
