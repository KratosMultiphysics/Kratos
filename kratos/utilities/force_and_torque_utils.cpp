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
#include "utilities/variable_utils.h"

namespace Kratos
{


array_1d<double,3> ForceAndTorqueUtils::SumForce(
    const ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rForceVariable)
{
    KRATOS_TRY

    return VariableUtils().SumHistoricalVariable<array_1d<double,3>>(
        rForceVariable,
        rModelPart);

    KRATOS_CATCH("")
}


std::array<array_1d<double,3>,2> ForceAndTorqueUtils::SumForceAndTorque(
    const ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rForceVariable,
    const Variable<array_1d<double,3>>& rTorqueVariable)
{
    KRATOS_TRY

    using CompoundArray = std::array<array_1d<double,3>,2>;
    using CompoundTuple = std::tuple<array_1d<double,3>,array_1d<double,3>>;
    using CompoundReduction = CombinedReduction<
        SumReduction<array_1d<double,3>>,
        SumReduction<array_1d<double,3>>
    >;

    auto force_and_moment = block_for_each<CompoundReduction>(
        rModelPart.GetCommunicator().LocalMesh().Nodes(),
        [&](const Node<3>& rNode) -> CompoundTuple {
            // {{fx, fy, fz},{mx, my, mz}}
            CompoundTuple force_and_moment_nodal;
            array_1d<double,3>& r_force = std::get<0>(force_and_moment_nodal);
            array_1d<double,3>& r_moment = std::get<1>(force_and_moment_nodal);
            
            r_force = rNode.FastGetSolutionStepValue(rForceVariable);
            r_moment = rNode.FastGetSolutionStepValue(rTorqueVariable);

            return force_and_moment_nodal;
        }
    );

    array_1d<double,3>& r_force = std::get<0>(force_and_moment);
    array_1d<double,3>& r_moment = std::get<1>(force_and_moment);

    r_force = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(r_force);
    r_moment = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(r_moment);

    return CompoundArray {r_force, r_moment};

    KRATOS_CATCH("")
}


std::array<array_1d<double,3>,2> ForceAndTorqueUtils::ComputeEquivalentForceAndTorque(
    const ModelPart& rModelPart,
    const array_1d<double,3>& rReferencePoint,
    const Variable<array_1d<double,3>>& rForceVariable,
    const Variable<array_1d<double,3>>& rTorqueVariable)
{
    KRATOS_TRY

    using CompoundArray = std::array<array_1d<double,3>,2>;
    using CompoundTuple = std::tuple<array_1d<double,3>,array_1d<double,3>>;
    using CompoundReduction = CombinedReduction<
        SumReduction<array_1d<double,3>>,
        SumReduction<array_1d<double,3>>
    >;

    std::function<CompoundTuple(const Node<3>&)> compute_force_and_moment;

    // Check whether nodes have the torque variable and set the lambda accordingly.
    // Note: multiple node types within the model part are not supported.
    //       In that case, this check would have to performed for each node.
    if (rModelPart.HasNodalSolutionStepVariable(rTorqueVariable)) {
        compute_force_and_moment = [&](const Node<3>& rNode) -> CompoundTuple
        {
            // {{fx, fy, fz},{mx, my, mz}}
            CompoundTuple force_and_moment_nodal;
            array_1d<double,3>& r_force = std::get<0>(force_and_moment_nodal);
            array_1d<double,3>& r_moment = std::get<1>(force_and_moment_nodal);
            
            r_force = rNode.FastGetSolutionStepValue(rForceVariable);
            r_moment = rNode.FastGetSolutionStepValue(rTorqueVariable);

            r_moment += MathUtils<double>::CrossProduct<array_1d<double,3>>(
                rNode-rReferencePoint,
                r_force
            );

            return force_and_moment_nodal;
        };
    }
    else {
        compute_force_and_moment = [&](const Node<3>& rNode) -> CompoundTuple
        {
            // {{fx, fy, fz},{mx, my, mz}}
            CompoundTuple force_and_moment_nodal;
            array_1d<double,3>& r_force = std::get<0>(force_and_moment_nodal);
            array_1d<double,3>& r_moment = std::get<1>(force_and_moment_nodal);

            r_force = rNode.FastGetSolutionStepValue(rForceVariable);
            r_moment = MathUtils<double>::CrossProduct<array_1d<double,3>>(
                rNode-rReferencePoint,
                r_force
            );

            return force_and_moment_nodal;
        };
    }

    auto force_and_moment = block_for_each<CompoundReduction>(
        rModelPart.GetCommunicator().LocalMesh().Nodes(),
        compute_force_and_moment
    );

    array_1d<double,3>& r_force = std::get<0>(force_and_moment);
    array_1d<double,3>& r_moment = std::get<1>(force_and_moment);

    r_force = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(r_force);
    r_moment = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(r_moment);

    return CompoundArray {r_force, r_moment};

    KRATOS_CATCH("")
}


} /* namespace Kratos */
