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
        auto it_node = mrModelPart.NodesBegin() + i;
        const auto& r_coordinates = it_node->Coordinates();
        double distance_projection = inner_prod(r_coordinates,free_stream_velocity);
        if (distance_projection<min_projection){
            min_projection=distance_projection;
            mrReferenceNode = *it_node;
        }

    }

    for (int i = 0; i < static_cast<int>(mrModelPart.Conditions().size()); i++) {
        auto it_cond = mrModelPart.ConditionsBegin() + i;
        auto& r_geometry = it_cond->GetGeometry();
        array_1d<double,3> aux_coordinates;
        r_geometry.PointLocalCoordinates(aux_coordinates, r_geometry.Center());
        array_1d<double,3> normal = it_cond->GetGeometry().Normal(aux_coordinates);

        double velocity_projection = inner_prod(normal, free_stream_velocity);
        if (velocity_projection < 0.0) {
            for (std::size_t i_node = 0; i_node < r_geometry.size(); i_node++){
                array_1d<double,3> relative_coordinates = r_geometry[i_node].Coordinates() - mrReferenceNode.Coordinates();
                double inlet_potential = inner_prod(relative_coordinates, free_stream_velocity);
                r_geometry[i_node].Fix(VELOCITY_POTENTIAL);
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = inlet_potential + mReferencePotential;
                if (r_geometry[i_node].SolutionStepsDataHas(ADJOINT_VELOCITY_POTENTIAL)){
                    r_geometry[i_node].Fix(ADJOINT_VELOCITY_POTENTIAL);
                    r_geometry[i_node].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL) = 0.0;
                }
            }
        }
        else {
            it_cond->SetValue(FREE_STREAM_VELOCITY, free_stream_velocity);
        }

    }
}


} // namespace Kratos.
