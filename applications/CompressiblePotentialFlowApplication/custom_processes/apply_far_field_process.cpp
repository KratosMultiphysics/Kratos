//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez
//

// Project includes
#include "apply_far_field_process.h"

namespace Kratos {

// Constructor for ApplyFarFieldProcess Process
ApplyFarFieldProcess::ApplyFarFieldProcess(ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart)
{}

void ApplyFarFieldProcess::Execute()
{
    double min_projection = std::numeric_limits<double>::epsilon();
    const auto free_stream_velocity = mrModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY];
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        ModelPart::NodeIterator it_node = mrModelPart.NodesBegin() + i;
        const auto& r_coordinates = it_node->Coordinates();
        double distance_projection = inner_prod(r_coordinates,free_stream_velocity);
        if (distance_projection<min_projection){
            min_projection=distance_projection;
        }

    }
    std::cout << min_projection << std::endl;
}
} // namespace Kratos.
