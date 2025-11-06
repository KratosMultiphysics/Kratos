// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Gennady Markelov
//

#pragma once

#include "custom_utilities/linear_nodal_extrapolator.h"
#include "includes/kratos_export_api.h"

#include <algorithm>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ExtrapolationUtilities
{
public:
    [[nodiscard]] static Matrix CalculateExtrapolationMatrix(const Geometry<Node>& rGeometry,
                                                             const GeometryData::IntegrationMethod& rIntegrationMethod)
    {
        KRATOS_TRY

        const auto extrapolator = LinearNodalExtrapolator{};
        return extrapolator.CalculateElementExtrapolationMatrix(rGeometry, rIntegrationMethod);

        KRATOS_CATCH("")
    }

    [[nodiscard]] static std::vector<std::optional<Vector>> CalculateNodalStresses(
        const std::vector<std::size_t>&        node_ids,
        const Geometry<Node>&                  rGeometry,
        const GeometryData::IntegrationMethod& rIntegrationMethod,
        const std::vector<Vector>&             rIntegrationPointStresses)
    {
        const auto               number_of_nodes = rGeometry.size();
        std::vector<std::size_t> element_node_ids(number_of_nodes);
        std::transform(rGeometry.begin(), rGeometry.end(), element_node_ids.begin(),
                       [](const auto& node) { return node.Id(); });

        const auto extrapolation_matrix = CalculateExtrapolationMatrix(rGeometry, rIntegrationMethod);

        std::vector<Vector> nodal_stresses(number_of_nodes, ZeroVector(rIntegrationPointStresses[0].size()));
        for (unsigned int node_index = 0; node_index < number_of_nodes; ++node_index) {
            for (unsigned int gp_index = 0; gp_index < rIntegrationPointStresses.size(); ++gp_index) {
                nodal_stresses[node_index] +=
                    extrapolation_matrix(node_index, gp_index) * rIntegrationPointStresses[gp_index];
            }
        }

        std::vector<std::optional<Vector>> result;
        for (unsigned int i = 0; i < node_ids.size(); ++i) {
            if (std::find(element_node_ids.begin(), element_node_ids.end(), node_ids[i]) !=
                element_node_ids.end()) {
                result.emplace_back(nodal_stresses[i]);
            } else {
                result.emplace_back(std::nullopt);
            }
        }

        return result;
    }
};

} // namespace Kratos