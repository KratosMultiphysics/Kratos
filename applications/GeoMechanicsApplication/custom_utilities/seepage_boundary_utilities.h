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

    // Adds the entries of an element's right-hand side that belong to WATER_PRESSURE degrees of
    // freedom onto their nodes.
    //
    // rDofs must be the element's own degrees of freedom, in the same order as rRightHandSide.
    static void AccumulateWaterPressureEntries(const std::vector<Dof<double>*>& rDofs,
                                               const Vector&                    rRightHandSide,
                                               NodalFlowMap&                    rNodalFlows);

    // Returns the nodal water flow for every unique node of the given elements, assembled from the
    // right-hand side of every element. For a Pw element that right-hand side is exactly the sum of
    // the permeability flow, the compressibility flow and the fluid body flow.
    static NodalFlowMap CalculateNodalWaterFlows(ModelPart::ElementsContainerType& rElements,
                                                 const ProcessInfo&                rProcessInfo);

    // Writes the nodal water flows onto the NODAL_WATER_FLOW solution-step variable of the model
    // part. Every node is set to zero first, so nodes absent from rNodalFlows (e.g. nodes without a
    // WATER_PRESSURE degree of freedom) hold a defined value rather than stale data. This makes the
    // assembled flow visualisable through the normal nodal output path.
    static void AssignNodalWaterFlows(ModelPart& rModelPart, const NodalFlowMap& rNodalFlows);

    // Returns the distinct nodes of every GeoSeepageCondition in the model part, sorted ascending
    // by node id. Nodes shared by adjacent conditions appear exactly once.
    //
    // The model part is non-const because the returned nodes must be mutable: the strategy fixes
    // and frees their WATER_PRESSURE degree of freedom.
    static std::vector<Node*> CollectSeepageNodes(ModelPart& rModelPart);

    // Switches at most one seepage node between a Dirichlet and a zero-flux Neumann boundary, and
    // returns whether it switched anything.
    //
    // A free node whose WATER_PRESSURE is positive is unsaturated, so it should not be draining:
    // the highest-pressure such node is fixed at zero pressure. Otherwise the fixed node with the
    // largest outflow is released. Fixing takes precedence over releasing, and ties are broken by
    // the lowest node id so the result is reproducible.
    static bool SwitchOneSeepageNodeIfNeeded(const std::vector<Node*>& rSeepageNodes,
                                             const NodalFlowMap&       rNodalFlows,
                                             int                       EchoLevel = 0);
};

} // namespace Kratos::Geo
