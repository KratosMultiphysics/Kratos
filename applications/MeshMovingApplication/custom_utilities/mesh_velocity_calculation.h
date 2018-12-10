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

#if !defined(KRATOS_MESH_VELOCITY_COMPUTATION_H_INCLUDED)
#define KRATOS_MESH_VELOCITY_COMPUTATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/mesh_moving_variables.h"
#include "utilities/time_discretization.h"

namespace Kratos {
namespace MeshVelocityCalculation {

template<class TGenAlpha>
void CalculateMeshVelocitiesGeneralizedAlpha(ModelPart& rModelPart,
                                             const TGenAlpha& rGenAlphaHelper)
{
    const double beta = rGenAlphaHelper.GetBeta();
    const double gamma = rGenAlphaHelper.GetGamma();
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();
    const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
    const double const_u = gamma / (delta_time * beta);
    const double const_v = 1.0 - gamma / beta;
    const double const_a = delta_time * (1.0 - gamma / (2.0 * beta));

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
        r_mesh_a0 = (1.0 / (delta_time * gamma)) * (r_mesh_v0 - r_mesh_v1) - ((1 - gamma) / gamma) * r_mesh_a1;
    }

    rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
    rModelPart.GetCommunicator().SynchronizeVariable(MESH_ACCELERATION);
}

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF1& rBDF);

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF2& rBDF);

inline void CalculateMeshVelocities(ModelPart& rModelPart,
                                    const TimeDiscretization::Newmark& rGenAlpha)
{
    CalculateMeshVelocitiesGeneralizedAlpha(rModelPart, rGenAlpha);
}

inline void CalculateMeshVelocities(ModelPart& rModelPart,
                                    const TimeDiscretization::Bossak& rGenAlpha)
{
    CalculateMeshVelocitiesGeneralizedAlpha(rModelPart, rGenAlpha);
}

inline void CalculateMeshVelocities(ModelPart& rModelPart,
                                    const TimeDiscretization::GeneralizedAlpha& rGenAlpha)
{
    CalculateMeshVelocitiesGeneralizedAlpha(rModelPart, rGenAlpha);
}


} // namespace MeshVelocityCalculation
} // namespace Kratos.

#endif // KRATOS_MESH_VELOCITY_COMPUTATION_H_INCLUDED  defined
