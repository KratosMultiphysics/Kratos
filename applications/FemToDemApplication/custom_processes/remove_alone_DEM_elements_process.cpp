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

#include "custom_processes/remove_alone_DEM_elements_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

RemoveAloneDEMElementsProcess::RemoveAloneDEMElementsProcess(
    ModelPart& rModelPart,
    const std::string& rFreeOrFix)
    : mrModelPart(rModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

void RemoveAloneDEMElementsProcess::Execute() 
{

}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos