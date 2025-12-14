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

#include "custom_elements/interface_element.h"
#include "custom_utilities/generic_utilities.h"
#include "custom_utilities/linear_nodal_extrapolator.h"
#include "element_utilities.hpp"
#include "geometries/geometry.h"
#include "geometry_utilities.h"
#include "includes/element.h"

namespace Kratos
{

Matrix ExtrapolationUtilities::CalculateExtrapolationMatrix(const Element& rElement)
{
    auto        p_interface_element = dynamic_cast<const InterfaceElement*>(&rElement);
    const auto& r_geometry_for_extrapolation =
        p_interface_element ? p_interface_element->GetMidGeometry() : rElement.GetGeometry();
    const auto integration_points = GeoElementUtilities::GetIntegrationPointsOf(rElement);

    const auto extrapolator = LinearNodalExtrapolator{};
    auto result = extrapolator.CalculateElementExtrapolationMatrix(r_geometry_for_extrapolation, integration_points);

    if (p_interface_element) {
        // The extrapolation matrix has been calculated for the mid-geometry. Since both sides of
        // the interface element need to use that matrix, we have to expand the resultant matrix.
        const auto number_of_rows_per_side = result.size1();
        result.resize(2 * number_of_rows_per_side, result.size2());
        subrange(result, number_of_rows_per_side, result.size1(), 0, result.size2()) =
            subrange(result, 0, number_of_rows_per_side, 0, result.size2());
    }

    KRATOS_ERROR_IF_NOT(result.size1() == rElement.GetGeometry().PointsNumber())
        << "A number of extrapolation matrix rows " << result.size1() << " is not equal to a number of nodes "
        << rElement.GetGeometry().PointsNumber() << " for element id " << rElement.Id() << std::endl;

    KRATOS_ERROR_IF_NOT(result.size2() == integration_points.size())
        << "A number of extrapolation matrix columns " << result.size2()
        << " is not equal to a number of integration points " << integration_points.size()
        << " for element id " << rElement.Id() << std::endl;

    return result;
}

std::vector<std::optional<Vector>> ExtrapolationUtilities::CalculateNodalVectors(
    const std::vector<std::size_t>& rNodeIds, const Element& rElement, const std::vector<Vector>& rVectorsAtIntegrationPoints)
{
    const auto& r_geometry           = rElement.GetGeometry();
    const auto  element_node_ids     = GenericUtilities::GetIdsFromEntityContents(r_geometry);
    const auto  extrapolation_matrix = CalculateExtrapolationMatrix(rElement);

    KRATOS_ERROR_IF_NOT(extrapolation_matrix.size2() == rVectorsAtIntegrationPoints.size())
        << "An extrapolation matrix size " << extrapolation_matrix.size2()
        << " is not equal to given stress vectors size " << rVectorsAtIntegrationPoints.size()
        << " for element Id " << rElement.Id() << std::endl;

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
