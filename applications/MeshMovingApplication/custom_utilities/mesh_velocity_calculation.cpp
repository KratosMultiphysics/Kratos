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
#include "includes/model_part.h"
#include "includes/mesh_moving_variables.h"
#include "mesh_velocity_calculation.h"

namespace Kratos {
namespace MeshVelocityCalculation {

namespace { // helpers-namespace
void CalculateMeshVelocitiesGeneralizedAlpha(ModelPart& rModelPart,
                                             const double Beta,
                                             const double Gamma)
{
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();
    const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
    const double const_u = Gamma / (delta_time * Beta);
    const double const_v = 1.0 - Gamma / Beta;
    const double const_a = delta_time * (1.0 - Gamma / (2.0 * Beta));

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        const auto it_node  = nodes_begin + i;
        const auto& r_mesh_u0 = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        auto&       r_mesh_v0 = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
        auto&       r_mesh_a0 = it_node->FastGetSolutionStepValue(MESH_ACCELERATION);

        const auto& r_mesh_u1 = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
        const auto& r_mesh_v1 = it_node->FastGetSolutionStepValue(MESH_VELOCITY,     1);
        const auto& r_mesh_a1 = it_node->FastGetSolutionStepValue(MESH_ACCELERATION, 1);

        r_mesh_v0 = const_u * (r_mesh_u0 - r_mesh_u1) + const_v * r_mesh_v1 + const_a * r_mesh_a1;
        r_mesh_a0 = (1.0 / (delta_time * Gamma)) * (r_mesh_v0 - r_mesh_v1) - ((1 - Gamma) / Gamma) * r_mesh_a1;

    }

    rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
    rModelPart.GetCommunicator().SynchronizeVariable(MESH_ACCELERATION);
}
}// helpers-namespace

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF1& rBDF)
{
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    const auto coeffs = rBDF.ComputeBDFCoefficients(rModelPart.GetProcessInfo());

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        const auto it_node  = nodes_begin + i;
        auto& r_mesh_v0       = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
        noalias(r_mesh_v0) = coeffs[0] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        noalias(r_mesh_v0) += coeffs[1] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
    }

    rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
}

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF2& rBDF)
{
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    const auto coeffs = rBDF.ComputeBDFCoefficients(rModelPart.GetProcessInfo());

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        const auto it_node  = nodes_begin + i;
        auto& r_mesh_v0 = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
        noalias(r_mesh_v0) = coeffs[0] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        noalias(r_mesh_v0) += coeffs[1] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
        noalias(r_mesh_v0) += coeffs[2] * it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 2);
    }

    rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
}

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::Newmark& rGenAlpha)
{
    CalculateMeshVelocitiesGeneralizedAlpha(rModelPart,
                                            rGenAlpha.GetBeta(),
                                            rGenAlpha.GetGamma());
}

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::Bossak& rGenAlpha)
{
    CalculateMeshVelocitiesGeneralizedAlpha(rModelPart,
                                            rGenAlpha.GetBeta(),
                                            rGenAlpha.GetGamma());
}

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::GeneralizedAlpha& rGenAlpha)
{
    CalculateMeshVelocitiesGeneralizedAlpha(rModelPart,
                                            rGenAlpha.GetBeta(),
                                            rGenAlpha.GetGamma());
}

} // namespace MeshVelocityCalculation.
} // namespace Kratos.
