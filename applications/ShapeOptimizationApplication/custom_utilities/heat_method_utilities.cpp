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

void HeatMethodUtilities::ComputeLaplacian() {
    KRATOS_TRY

    // TODO
    KRATOS_CATCH("");
}

}  // namespace Kratos.
