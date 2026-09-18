// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse,
//                   Wijtze Pieter Kikstra

#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "includes/define.h"
#include "includes/dof.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos::Geo
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SeepageBoundaryUtilities
{
public:
    // Nodal water flow, keyed by node id.
    using NodalFlowMap = std::unordered_map<std::size_t, double>;

    /**
     * @brief Returns the nodal water flow for every unique node of the given elements, assembled
     * from the right-hand side of every element. For a Pw element that right-hand side is exactly
     * the sum of the permeability flow, the compressibility flow and the fluid body flow.
     */
    static NodalFlowMap CalculateNodalWaterFlows(ModelPart::ElementsContainerType& rElements,
                                                 const ProcessInfo&                rProcessInfo);

    /**
     * @brief Writes the nodal water flows onto the NODAL_WATER_FLOW solution-step variable of the
     * model part. Every node is set to zero first, so nodes absent from rNodalFlows (e.g. nodes
     * without a WATER_PRESSURE degree of freedom) hold a defined value rather than stale data. This
     * makes the assembled flow visualisable through the normal nodal output path.
     */
    static void AssignNodalWaterFlows(ModelPart& rModelPart, const NodalFlowMap& rNodalFlows);

    /**
     * @brief Returns the distinct nodes of every GeoSeepageCondition in the model part, sorted
     * ascending by node id. Nodes shared by adjacent conditions appear exactly once.
     *
     * The model part is non-const because the returned nodes must be mutable: the strategy fixes
     * and frees their WATER_PRESSURE degree of freedom.
     */
    static std::vector<Node*> CollectSeepageNodes(ModelPart& rModelPart);

    /**
     * @brief Switches at most one seepage node between a Dirichlet and a zero-flux Neumann
     * boundary, and returns whether it switched anything.
     *
     * A free node with negative WATER_PRESSURE is selected for a zero-pressure Dirichlet boundary;
     * the node with most-negative water pressure is chosen. Otherwise, the fixed node with the
     * largest inflow is released. Fixing takes precedence over releasing, and when equally suitable
     * candidates are found, the one with the lowest node ID will be selected, just so the result is
     * reproducible.
     */
    static bool SwitchOneSeepageNodeIfNeeded(const std::vector<Node*>& rSeepageNodes,
                                             const NodalFlowMap&       rNodalFlows,
                                             int                       EchoLevel = 0);
};

} // namespace Kratos::Geo
