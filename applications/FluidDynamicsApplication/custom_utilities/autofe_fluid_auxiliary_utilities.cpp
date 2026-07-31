//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolas Sibuet
//

// System includes


// External includes


// Project includes
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/variable_utils.h"

// Application includes
#include "autofe_fluid_auxiliary_utilities.h"

namespace Kratos
{

void AutoFEFluidAuxiliaryUtilities::PostprocessP2P1ContinuousScalarVariable(
    ModelPart& rModelPart,
    const Variable<double>& rVariable
)
{
    // Reset VISITED flag to indicate the already postprocessed nodes
    VariableUtils().SetFlag(VISITED, false, rModelPart.Nodes());

    // Loop the P2P1 elements to postprocess the pressure in edge midpoint nodes
    // Note that there is no need to do MPI synchronization as we are updating the ghost nodes too
    if (rModelPart.GetCommunicator().LocalMesh().NumberOfElements() != 0) {
        // Check DOMAIN_SIZE
        KRATOS_ERROR_IF_NOT(rModelPart.GetProcessInfo().Has(DOMAIN_SIZE))
            << "DOMAIN_SIZE cannot be found in '" << rModelPart.FullName() << "' ProcessInfo container."<<  std::endl;
        // Loop the elements to assign the edge midpoint nodes PRESSURE
        if (rModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2) {
            block_for_each(rModelPart.Elements(), [&rVariable](auto& rElement){
                auto& r_geometry = rElement.GetGeometry();
                KRATOS_ERROR_IF_NOT(r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 3, 0, 1);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 4, 1, 2);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 5, 0, 2);
            });
        } else {
            block_for_each(rModelPart.Elements(), [&rVariable](auto& rElement){
                auto& r_geometry = rElement.GetGeometry();
                KRATOS_ERROR_IF_NOT(r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 4, 0, 1);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 5, 1, 2);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 6, 0, 2);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 7, 0, 3);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 8, 1, 3);
                PostprocessP2P1NodeScalarVariable(r_geometry, rVariable, 9, 2, 3);
            });
        }
    }
}

} // namespace Kratos
