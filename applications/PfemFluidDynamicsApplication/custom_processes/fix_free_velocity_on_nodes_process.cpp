//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics PfemFluidDynamics Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Carlos Eulogio Flores
//

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "custom_processes/fix_free_velocity_on_nodes_process.h"
#include "utilities/variable_utils.h"

namespace Kratos {

PFEMFixFreeVelocityOnNodesProcess::PFEMFixFreeVelocityOnNodesProcess(
    ModelPart& rModelPart,
    const bool rFreeOrFix)
    : mrModelPart(rModelPart),
      mFreeOrFix(rFreeOrFix)
{
}

/***********************************************************************************/
/***********************************************************************************/

void PFEMFixFreeVelocityOnNodesProcess::Execute()
{
    VariableUtils().ApplyFixity(VELOCITY_X, mFreeOrFix, mrModelPart.Nodes());
    VariableUtils().ApplyFixity(VELOCITY_Y, mFreeOrFix, mrModelPart.Nodes());
    VariableUtils().ApplyFixity(VELOCITY_Z, mFreeOrFix, mrModelPart.Nodes());
}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos
