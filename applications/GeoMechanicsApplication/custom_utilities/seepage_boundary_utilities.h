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

} // namespace Kratos::Geo::SeepageBoundaryUtilities


