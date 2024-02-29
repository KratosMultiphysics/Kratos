// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <functional>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "geometry_utilities.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "utilities/builtin_timer.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "includes/global_pointer_variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/atomic_utilities.h"
#include "geometries/geometry_data.h"
#include "processes/calculate_nodal_area_process.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/heat_method_utilities.h"




// ==============================================================================

namespace Kratos
{

void HeatMethodUtilities::ComputeGeodesicDistance()
{

    KRATOS_TRY;

    // Non parallalized loop
    // for (auto& r_node_i : mrModelPart.Nodes()) {

    //     array_3d& heat_diffusion = r_node_i.FastGetSolutionStepValue(HEAT_DIFFUSION);
    //     heat_diffusion[0] = r_node_i.Id();
    //     heat_diffusion[1] = r_node_i.Id() + 1;
    //     heat_diffusion[2] = r_node_i.Id() + 2;

    //     double& geodesic_distance = r_node_i.FastGetSolutionStepValue(GEODESIC_DISTANCE);
    //     geodesic_distance = r_node_i.Id();
    // }

    // Parallalized loop
    block_for_each(mrModelPart.Nodes(), [&](Node& rNode)
    {
        array_3d& heat_diffusion = rNode.FastGetSolutionStepValue(HEAT_DIFFUSION);
        heat_diffusion[0] = rNode.Id();
        heat_diffusion[1] = rNode.Id() + 1;
        heat_diffusion[2] = rNode.Id() + 2;

        double& geodesic_distance = rNode.FastGetSolutionStepValue(GEODESIC_DISTANCE);
        geodesic_distance = rNode.Id();
    });


    // mrSolver.Solve(K, U, U0);
    KRATOS_CATCH("");
}

}  // namespace Kratos.