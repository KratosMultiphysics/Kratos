//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "mesh_velocity_calculation.h"

namespace Kratos {
namespace MeshVelocityCalculation {

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF1& rBDF)
{
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

    const auto coeffs = rBDF.ComputeBDFCoefficients(delta_time);

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        const auto it_node  = nodes_begin + i;
        auto& r_mesh_v0       = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
        noalias(r_mesh_v0)  = coeffs[0] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        noalias(r_mesh_v0) += coeffs[1] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
    }

    rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
}

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF2& rBDF)
{
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    const auto& r_current_process_info = rModelPart.GetProcessInfo();

    const double delta_time = r_current_process_info[DELTA_TIME];
    const double previous_delta_time = r_current_process_info.GetPreviousTimeStepInfo(1)[DELTA_TIME];

    const auto coeffs = rBDF.ComputeBDFCoefficients(delta_time, previous_delta_time);

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        const auto it_node  = nodes_begin + i;
        auto& r_mesh_v0 = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
        noalias(r_mesh_v0)  = coeffs[0] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        noalias(r_mesh_v0) += coeffs[1] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
        noalias(r_mesh_v0) += coeffs[2] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 2);
    }

    rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
}

} // namespace MeshVelocityCalculation.
} // namespace Kratos.
