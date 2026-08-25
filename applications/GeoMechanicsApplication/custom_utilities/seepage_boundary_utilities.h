// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "includes/define.h"
#include "includes/dof.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos::Geo::SeepageBoundaryUtilities
{

// Nodal water flow, keyed by node id.
using NodalFlowMap = std::unordered_map<std::size_t, double>;

// Adds the entries of an element's right-hand side that belong to WATER_PRESSURE degrees of freedom
// onto their nodes.
//
// rDofs must be the element's own degrees of freedom, in the same order as rRightHandSide. For a
// U-Pw element the displacement and pressure entries are interleaved, so the entries cannot be
// matched to nodes positionally; the degrees of freedom are what makes the mapping correct.
void KRATOS_API(GEO_MECHANICS_APPLICATION) AccumulateWaterPressureEntries(
    const std::vector<Dof<double>*>& rDofs, const Vector& rRightHandSide, NodalFlowMap& rNodalFlows);

// Returns the nodal water flow for every node in the model part, assembled from the right-hand side
// of every element. For a Pw element that right-hand side is exactly the permeability flow plus the
// compressibility flow plus the fluid body flow.
NodalFlowMap KRATOS_API(GEO_MECHANICS_APPLICATION)
    CalculateNodalWaterFlows(ModelPart& rModelPart, const ProcessInfo& rProcessInfo);

// Returns the distinct nodes of every GeoSeepageCondition in the model part, sorted ascending by
// node id. Nodes shared by adjacent conditions appear exactly once.
//
// The model part is non-const because the returned nodes must be mutable: the strategy fixes and
// frees their WATER_PRESSURE degree of freedom.
std::vector<Node*> KRATOS_API(GEO_MECHANICS_APPLICATION) CollectSeepageNodes(ModelPart& rModelPart);

// Returns true if a Dirichlet seepage node carrying this nodal flow should be released to a
// zero-flux Neumann boundary.
//
// This sign is the key modelling assumption of the prototype: a positive nodal flow means outflow,
// and the node carrying the largest outflow is released. Confirming or inverting it is what the
// Muskat validation case is for, so it deliberately lives in one place.
[[nodiscard]] bool KRATOS_API(GEO_MECHANICS_APPLICATION) ShouldReleaseToNeumann(double NodalFlow);

// Switches at most one seepage node between a Dirichlet and a zero-flux Neumann boundary, and
// returns whether it switched anything.
//
// A free node whose WATER_PRESSURE is positive is unsaturated, so it should not be draining: the
// highest-pressure such node is fixed at zero pressure. Otherwise the fixed node with the largest
// outflow is released. Fixing takes precedence over releasing, and ties are broken by the lowest
// node id so the result is reproducible.
bool KRATOS_API(GEO_MECHANICS_APPLICATION)
    SwitchOneSeepageNode(const std::vector<Node*>& rSeepageNodes, const NodalFlowMap& rNodalFlows);

} // namespace Kratos::Geo::SeepageBoundaryUtilities



