//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//	                Kratos default license: kratos/license.txt
//
//  Main Authors:   Máté Kelemen
//

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/force_and_torque_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{


std::array<array_1d<double,3>,2> ForceAndTorqueUtils::SumForceAndTorque(const ModelPart& rModelPart,
                                                                        const array_1d<double,3>& rReferencePoint)
{
    KRATOS_TRY

    using CompoundArray = array_1d<double,6>;

    CompoundArray force_and_moment = block_for_each<SumReduction<CompoundArray>>(
        rModelPart.GetCommunicator().LocalMesh().Nodes(),
        [&rReferencePoint](const Node<3>& rNode) -> CompoundArray
        {
            array_1d<double,3> force, moment;

            force = rNode.GetSolutionStepValue(REACTION);

            if (rNode.SolutionStepsDataHas(MOMENT)) {
                moment = rNode.GetSolutionStepValue(MOMENT);
            }
            else {
                moment[0] = 0.0;
                moment[1] = 0.0;
                moment[2] = 0.0;
            }

            moment += MathUtils<double>::CrossProduct<array_1d<double,3>>(
                rNode-rReferencePoint,
                force
            );

            return CompoundArray {force[0], force[1], force[2],
                                  moment[0], moment[1], moment[2]};
        }
    );

    std::array<array_1d<double,3>,2> output;
    output[0][0] = force_and_moment[0];
    output[0][1] = force_and_moment[1];
    output[0][2] = force_and_moment[2];
    output[1][0] = force_and_moment[3];
    output[1][1] = force_and_moment[4];
    output[1][2] = force_and_moment[5];

    output[0] = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(output[0]);
    output[1] = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(output[1]);

    return output;

    KRATOS_CATCH("")
}


} /* namespace Kratos */