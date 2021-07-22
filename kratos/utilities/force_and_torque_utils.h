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

#ifndef KRATOS_FORCE_AND_TORQUE_UTILS
#define KRATOS_FORCE_AND_TORQUE_UTILS

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * 
 */
class KRATOS_API(KRATOS_CORE) ForceAndTorqueUtils
{
public:
    ///@name Operations
    ///@{

    /** Sum forces on all nodes
     *  @param rModelPart model part containing all nodes to perform the sum on
     *  @param rForceVariable nodal force variable to be summed up (example: REACTION)
     *  @note this function is identical to @ref{VariableUtils::SumHistoricalVariable} and calls it internally
     */
    static array_1d<double,3> SumForce(
        const ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rForceVariable);

    /** Sum forces and torques on all nodes
     *  @param rModelPart model part containing all nodes to perform the sum on
     *  @param rForceVariable nodal force variable to be summed up (example: REACTION)
     *  @param rTorqueVariable nodal torque variable to be summed up (example: MOMENT)
     */
    static std::array<array_1d<double,3>,2> SumForceAndTorque(
        const ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rForceVariable,
        const Variable<array_1d<double,3>>& rTorqueVariable);

    /** Reduce total force and moment exerted by the model part on a reference point
     *  @param rModelPart model part containing all nodes to perform the sum on
     *  @param rReferencePoint reference point for computing moments
     *  @param rForceVariable nodal force variable to be summed up (example: REACTION)
     *  @param rTorqueVariable nodal torque variable to be summed up (example: MOMENT)
     */
    static std::array<array_1d<double,3>,2> ComputeEquivalentForceAndTorque(
        const ModelPart& rModelPart,
        const array_1d<double,3>& rReferencePoint,
        const Variable<array_1d<double,3>>& rForceVariable,
        const Variable<array_1d<double,3>>& rTorqueVariable);

    ///@}
}; /* class ForceAndTorqueUtils */

///@}
} /* namespace Kratos */

#endif /* KRATOS_FORCE_AND_TORQUE_UTILS  defined */