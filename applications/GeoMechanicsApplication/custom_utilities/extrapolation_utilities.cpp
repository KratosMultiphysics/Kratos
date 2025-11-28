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
#include "geometry_utilities.h"
#include "includes/element.h"
#include "includes/node.h"

namespace Kratos
{

Matrix ExtrapolationUtilities::CalculateExtrapolationMatrix(const Element& rElement)
{
    const auto& r_geometry         = rElement.GetGeometry();
    const auto  integration_method = rElement.GetIntegrationMethod();
    const auto  element_id         = rElement.Id();

    const auto extrapolator = LinearNodalExtrapolator{};
    const auto result       = extrapolator.CalculateElementExtrapolationMatrix(rElement);

    KRATOS_ERROR_IF_NOT(result.size1() == r_geometry.PointsNumber())
        << "A number of extrapolation matrix rows " << result.size1() << " is not equal to a number of nodes "
        << r_geometry.PointsNumber() << " for element id " << element_id << std::endl;

    KRATOS_ERROR_IF_NOT(result.size2() == r_geometry.IntegrationPointsNumber(integration_method))
        << "A number of extrapolation matrix columns " << result.size2()
        << " is not equal to a number of integration points "
        << r_geometry.IntegrationPointsNumber(integration_method) << " for element id "
        << element_id << std::endl;

    return result;
}

std::vector<std::optional<Vector>> ExtrapolationUtilities::CalculateNodalVectors(
    const std::vector<std::size_t>& rNodeIds, const Element& rElement, const std::vector<Vector>& rVectorsAtIntegrationPoints)
{
    const auto& r_geometry = rElement.GetGeometry();
    const auto  element_id = rElement.Id();

    const auto element_node_ids     = GeometryUtilities::GetNodeIdsFromGeometry(r_geometry);
    const auto extrapolation_matrix = CalculateExtrapolationMatrix(rElement);

    KRATOS_ERROR_IF_NOT(extrapolation_matrix.size2() == rVectorsAtIntegrationPoints.size())
        << "An extrapolation matrix size " << extrapolation_matrix.size2()
        << " is not equal to given stress vectors size " << rVectorsAtIntegrationPoints.size()
        << " for element Id " << element_id << std::endl;

    const auto          number_of_nodes = r_geometry.PointsNumber();
    std::vector<Vector> extrapolated_vectors_at_nodes(
        number_of_nodes, ZeroVector(rVectorsAtIntegrationPoints[0].size()));
    for (unsigned int node_index = 0; node_index < number_of_nodes; ++node_index) {
        for (unsigned int integration_point = 0;
             integration_point < rVectorsAtIntegrationPoints.size(); ++integration_point) {
            extrapolated_vectors_at_nodes[node_index] += extrapolation_matrix(node_index, integration_point) *
                                                         rVectorsAtIntegrationPoints[integration_point];
        }
    }

    std::vector<std::optional<Vector>> result;
    for (const auto& node_id : rNodeIds) {
        auto it = std::find(element_node_ids.begin(), element_node_ids.end(), node_id);
        if (it != element_node_ids.end()) {
            result.emplace_back(extrapolated_vectors_at_nodes[std::distance(element_node_ids.begin(), it)]);
        } else {
            result.emplace_back(std::nullopt);
        }
    }

    return result;
}

} // namespace Kratos