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

// Returns the distinct nodes of every GeoSeepageCondition in the model part, sorted ascending by
// node id. Nodes shared by adjacent conditions appear exactly once.
//
// The model part is non-const because the returned nodes must be mutable: the strategy fixes and
// frees their WATER_PRESSURE degree of freedom.
std::vector<Node*> KRATOS_API(GEO_MECHANICS_APPLICATION) CollectSeepageNodes(ModelPart& rModelPart);

} // namespace Kratos::Geo::SeepageBoundaryUtilities

