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
    template <unsigned int TNumNodes>
    static void CalculateExtrapolationMatrix(const Geometry<Node>& rGeometry,
                                             const GeometryData::IntegrationMethod& rIntegrationMethod,
                                             BoundedMatrix<double, TNumNodes, TNumNodes>& rExtrapolationMatrix)
    {
        KRATOS_TRY

        const auto extrapolator = LinearNodalExtrapolator{};
        const auto result = extrapolator.CalculateElementExtrapolationMatrix(rGeometry, rIntegrationMethod);
        KRATOS_ERROR_IF_NOT(result.size1() == TNumNodes)
            << "Extrapolation matrix has unexpected number of rows: " << result.size1()
            << " (expected " << TNumNodes << ")" << std::endl;
        KRATOS_ERROR_IF_NOT(result.size2() == TNumNodes)
            << "Extrapolation matrix has unexpected number of columns: " << result.size2()
            << " (expected " << TNumNodes << ")" << std::endl;
        rExtrapolationMatrix = result;

        KRATOS_CATCH("")
    }

    template <unsigned int TNumNodes>
    static std::vector<Vector> CalculateNodalStresses(const std::vector<std::size_t>& node_ids,
                                                      const Geometry<Node>&           rGeometry,
                                                      const GeometryData::IntegrationMethod& rIntegrationMethod,
                                                      std::vector<Vector>& rIntegrationPointStresses)
    {
        std::vector<std::size_t> element_node_ids(TNumNodes);
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            element_node_ids[i] = rGeometry[i].Id();
        }

        BoundedMatrix<double, TNumNodes, TNumNodes> extrapolation_matrix;
        CalculateExtrapolationMatrix(rGeometry, rIntegrationMethod, extrapolation_matrix);

        std::vector<Vector> nodal_stresses(TNumNodes, ZeroVector(rIntegrationPointStresses[0].size()));
        for (unsigned int node_index = 0; node_index < TNumNodes; ++node_index) {
            for (unsigned int gp_index = 0; gp_index < rIntegrationPointStresses.size(); ++gp_index) {
                nodal_stresses[node_index] +=
                    extrapolation_matrix(node_index, gp_index) * rIntegrationPointStresses[gp_index];
            }
        }

        std::vector<Vector> result;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            std::size_t node_id = node_ids[i];
            if (std::find(node_ids.begin(), node_ids.end(), node_id) != node_ids.end()) {
                result.emplace_back(nodal_stresses[i]);
            }
        }

        return result;
    }
};

} // namespace Kratos