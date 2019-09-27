//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/fix_free_velocity_on_nodes_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

FixFreeVelocityOnNodesProcess::FixFreeVelocityOnNodesProcess(
    ModelPart& rModelPart,
    const std::string& rFreeOrFix)
    : mrModelPart(rModelPart),
      mFreeOrFix(rFreeOrFix)
{
}

/***********************************************************************************/
/***********************************************************************************/

void FixFreeVelocityOnNodesProcess::Execute() 
{

}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos