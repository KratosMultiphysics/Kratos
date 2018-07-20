// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "custom_processes/solid_shell_thickness_compute_process.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

void SolidShellThickComputeProcess::Execute()
{
    KRATOS_TRY

    // We set to zero the thickness
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(THICKNESS, 0.0, nodes_array);

    // Connectivity map
    std::unordered_map<IndexType, IndexType> connectivity_map;

    // Now we iterate over the elements to create a connectivity map
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

//     #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i){
        const auto it_elem = elements_array.begin() + i;

        // We get the condition geometry
        GeometryType& r_this_geometry = it_elem->GetGeometry();

        if (r_this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) {
            connectivity_map.insert({r_this_geometry[0].Id(), r_this_geometry[3].Id()});
            connectivity_map.insert({r_this_geometry[1].Id(), r_this_geometry[4].Id()});
            connectivity_map.insert({r_this_geometry[2].Id(), r_this_geometry[5].Id()});
        } else if (r_this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) {
            // NOTE: It will be assumed that it keeps the default connectivity of GiD, please check it or it will return whatever
//             KRATOS_WARNING_ONCE("SolidShellThickComputeProcess") << "It will be assumed that it keeps the default connectivity of GiD, please check it or it will return whatever" << std::endl;
            connectivity_map.insert({r_this_geometry[0].Id(), r_this_geometry[4].Id()});
            connectivity_map.insert({r_this_geometry[1].Id(), r_this_geometry[5].Id()});
            connectivity_map.insert({r_this_geometry[2].Id(), r_this_geometry[6].Id()});
            connectivity_map.insert({r_this_geometry[3].Id(), r_this_geometry[7].Id()});
        } else {
            KRATOS_ERROR << "You are not using a geometry that is supposed to be a solid-shell. Check that you assigned a proper model part (just containing the solid-shells elements)" << std::endl;
        }
    }

    // Now that the connectivity has been constructed
    double distance;
    for (auto& pair : connectivity_map) {
        auto pnode1 = mrThisModelPart.pGetNode(pair.first);
        auto pnode2 = mrThisModelPart.pGetNode(pair.second);
        distance = norm_2(pnode1->Coordinates() - pnode2->Coordinates());

        const double previous_thickness_1 = pnode1->GetValue(THICKNESS);
        const double previous_thickness_2 = pnode2->GetValue(THICKNESS);

        // We set the thickness on the firt node
        if (previous_thickness_1 > 0.0) {
            pnode1->SetValue(THICKNESS, (previous_thickness_1 + distance));
        } else {
            pnode1->SetValue(THICKNESS, distance);
        }
        // We set the thickness on the second node
        if (previous_thickness_2 > 0.0) {
            pnode2->SetValue(THICKNESS, (previous_thickness_2 + distance));
        } else {
            pnode2->SetValue(THICKNESS, distance);
        }
    }

    KRATOS_CATCH("")
} // class SolidShellThickComputeProcess
} // namespace Kratos
