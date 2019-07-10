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
    mFreeStreamVelocity = mrModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY];

    std::size_t num_threads = omp_get_max_threads();
    std::vector<double> min_projections(num_threads, std::numeric_limits<double>::epsilon());
    std::vector<std::size_t> nodes_id_list(num_threads, 0);

    double distance_projection=0.0;

    #pragma omp parallel for firstprivate (distance_projection)
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = mrModelPart.NodesBegin() + i;
        std::size_t thread_id = omp_get_thread_num();
        const auto& r_coordinates = it_node->Coordinates();
        distance_projection = inner_prod(r_coordinates,mFreeStreamVelocity);
        if (distance_projection < min_projections[thread_id]){
            min_projections[thread_id] = distance_projection;
            nodes_id_list[thread_id] = it_node->Id();
        }
    }

    std::size_t minimum_node_thread_id = 0;
    for (std::size_t i_thread = 0; i_thread<num_threads; i_thread++){
        if (min_projections[i_thread] < min_projections[minimum_node_thread_id]){
            minimum_node_thread_id = i_thread;
        }
    }

    mpReferenceNode = mrModelPart.pGetNode(nodes_id_list[minimum_node_thread_id]);

    double velocity_projection = 0.0;
    double inlet_potential = 0.0;
    array_1d<double,3> normal;
    array_1d<double,3> aux_coordinates;
    array_1d<double,3> relative_coordinates;

    #pragma omp parallel for firstprivate(velocity_projection, normal, aux_coordinates, inlet_potential, relative_coordinates)
    for (int i = 0; i < static_cast<int>(mrModelPart.Conditions().size()); i++) {
        auto it_cond = mrModelPart.ConditionsBegin() + i;
        auto& r_geometry = it_cond->GetGeometry();
        r_geometry.PointLocalCoordinates(aux_coordinates, r_geometry.Center());
        normal = it_cond->GetGeometry().Normal(aux_coordinates);

        velocity_projection = inner_prod(normal, mFreeStreamVelocity);
        if (velocity_projection < 0.0) {
            for (std::size_t i_node = 0; i_node < r_geometry.size(); i_node++){
                relative_coordinates = r_geometry[i_node].Coordinates() - mpReferenceNode->Coordinates();
                inlet_potential = inner_prod(relative_coordinates, mFreeStreamVelocity);
                r_geometry[i_node].Fix(VELOCITY_POTENTIAL);
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = inlet_potential + mReferencePotential;
                if (r_geometry[i_node].SolutionStepsDataHas(ADJOINT_VELOCITY_POTENTIAL)){
                    r_geometry[i_node].Fix(ADJOINT_VELOCITY_POTENTIAL);
                    r_geometry[i_node].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL) = 0.0;
                }
            }
        }
        else {
            it_cond->SetValue(FREE_STREAM_VELOCITY, mFreeStreamVelocity);
        }

    }
}

void ApplyFarFieldProcess::InitializeFlowField()
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();
    double inlet_potential = 0.0;
    array_1d<double,3> relative_coordinates;

    #pragma omp parallel for firstprivate(inlet_potential,relative_coordinates)
    for (int i = 0; i < static_cast<int>(root_model_part.Nodes().size()); i++) {
        auto it_node = root_model_part.NodesBegin() + i;

        relative_coordinates = it_node->Coordinates() - mpReferenceNode->Coordinates();
        inlet_potential = inner_prod(relative_coordinates, mFreeStreamVelocity);

        it_node->FastGetSolutionStepValue(VELOCITY_POTENTIAL) = inlet_potential + mReferencePotential;
        it_node->FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = inlet_potential + mReferencePotential;
    }
}

} // namespace Kratos.
