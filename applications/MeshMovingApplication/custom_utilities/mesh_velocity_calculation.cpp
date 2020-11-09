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
#include "utilities/parallel_utilities.h"

namespace Kratos {
namespace MeshVelocityCalculation {

namespace { // helpers-namespace
void CalculateMeshVelocitiesGeneralizedAlpha(ModelPart& rModelPart,
                                             const double Beta,
                                             const double Gamma)
{
    const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
    const double const_u = Gamma / (delta_time * Beta);
    const double const_v = 1.0 - Gamma / Beta;
    const double const_a = delta_time * (1.0 - Gamma / (2.0 * Beta));

    block_for_each(rModelPart.GetCommunicator().LocalMesh().Nodes(),
        [&](Node<3>& rNode){
            const auto& r_mesh_u0 = rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT);
            auto&       r_mesh_v0 = rNode.FastGetSolutionStepValue(MESH_VELOCITY);
            auto&       r_mesh_a0 = rNode.FastGetSolutionStepValue(MESH_ACCELERATION);

            const auto& r_mesh_u1 = rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
            const auto& r_mesh_v1 = rNode.FastGetSolutionStepValue(MESH_VELOCITY,     1);
            const auto& r_mesh_a1 = rNode.FastGetSolutionStepValue(MESH_ACCELERATION, 1);

            r_mesh_v0 = const_u * (r_mesh_u0 - r_mesh_u1) + const_v * r_mesh_v1 + const_a * r_mesh_a1;
            r_mesh_a0 = (1.0 / (delta_time * Gamma)) * (r_mesh_v0 - r_mesh_v1) - ((1 - Gamma) / Gamma) * r_mesh_a1;
    });

    rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
    rModelPart.GetCommunicator().SynchronizeVariable(MESH_ACCELERATION);
}
}// helpers-namespace

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF1& rBDF)
{
    const auto coeffs = rBDF.ComputeBDFCoefficients(rModelPart.GetProcessInfo());

    block_for_each(rModelPart.GetCommunicator().LocalMesh().Nodes(),
        [&]( Node<3>& rNode ){
            auto& r_mesh_v0     = rNode.FastGetSolutionStepValue(MESH_VELOCITY);
            noalias(r_mesh_v0)  = coeffs[0] * rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT);
            noalias(r_mesh_v0) += coeffs[1] * rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
    });

    rModelPart.GetCommunicator().SynchronizeVariable(MESH_VELOCITY);
}

void CalculateMeshVelocities(ModelPart& rModelPart,
                             const TimeDiscretization::BDF2& rBDF)
{
    const auto coeffs = rBDF.ComputeBDFCoefficients(rModelPart.GetProcessInfo());

    block_for_each(rModelPart.GetCommunicator().LocalMesh().Nodes(),
        [&](Node<3>& rNode){
            auto& r_mesh_v0     = rNode.FastGetSolutionStepValue(MESH_VELOCITY);
            noalias(r_mesh_v0)  = coeffs[0] * rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT);
            noalias(r_mesh_v0) += coeffs[1] * rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
            noalias(r_mesh_v0) += coeffs[2] * rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT, 2);
    });

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
