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

#include "custom_utilities/extrapolation_utilities.h"
#include "custom_utilities/linear_nodal_extrapolator.h"
#include "includes/node.h"

namespace Kratos
{

Matrix ExtrapolationUtilities::CalculateExtrapolationMatrix(const Geometry<Node>& rGeometry,
                                                            const GeometryData::IntegrationMethod& rIntegrationMethod,
                                                            size_t ElementId)
{
    KRATOS_TRY

    const auto extrapolator = LinearNodalExtrapolator{};

    const auto result = extrapolator.CalculateElementExtrapolationMatrix(rGeometry, rIntegrationMethod);

    KRATOS_ERROR_IF_NOT(result.size1() == rGeometry.size())
        << "An extrapolation matrix size1 " << result.size1() << " is not equal to a number of nodes "
        << rGeometry.size() << " for element id " << ElementId << std::endl;

    KRATOS_ERROR_IF_NOT(result.size2() == rGeometry.IntegrationPoints(rIntegrationMethod).size())
        << "An extrapolation matrix size2 " << result.size2() << " is not equal to a number of integration points "
        << rGeometry.IntegrationPoints(rIntegrationMethod).size() << " for element id " << ElementId
        << std::endl;

    return result;

    KRATOS_CATCH("")
}

std::vector<std::optional<Vector>> ExtrapolationUtilities::CalculateNodalStresses(
    const std::vector<std::size_t>&        node_ids,
    const Geometry<Node>&                  rGeometry,
    const GeometryData::IntegrationMethod& rIntegrationMethod,
    const std::vector<Vector>&             rIntegrationPointStresses,
    size_t                                 ElementId)
{
    const auto               number_of_nodes = rGeometry.size();
    std::vector<std::size_t> element_node_ids(number_of_nodes);
    std::transform(rGeometry.begin(), rGeometry.end(), element_node_ids.begin(),
                   [](const auto& node) { return node.Id(); });

    const auto extrapolation_matrix = CalculateExtrapolationMatrix(rGeometry, rIntegrationMethod, ElementId);

    KRATOS_ERROR_IF_NOT(extrapolation_matrix.size2() == rIntegrationPointStresses.size())
        << "An extrapolation matrix size " << extrapolation_matrix.size2()
        << " is not equal to given stress vectors size " << rIntegrationPointStresses.size()
        << " for element Id " << ElementId << std::endl;

    std::vector<Vector> nodal_stresses(number_of_nodes, ZeroVector(rIntegrationPointStresses[0].size()));
    for (unsigned int node_index = 0; node_index < number_of_nodes; ++node_index) {
        for (unsigned int integration_point = 0;
             integration_point < rIntegrationPointStresses.size(); ++integration_point) {
            nodal_stresses[node_index] += extrapolation_matrix(node_index, integration_point) *
                                          rIntegrationPointStresses[integration_point];
        }
    }

    std::vector<std::optional<Vector>> result;
    for (const auto& node_id : node_ids) {
        auto it = std::find(element_node_ids.begin(), element_node_ids.end(), node_id);
        if (it != element_node_ids.end()) {
            result.emplace_back(nodal_stresses[std::distance(element_node_ids.begin(), it)]);
        } else {
            result.emplace_back(std::nullopt);
        }
    }

    return result;
}

} // namespace Kratos