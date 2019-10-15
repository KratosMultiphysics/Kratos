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

#include "custom_processes/compute_initial_volume_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

ComputeInitialVolumeProcess::ComputeInitialVolumeProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart) 
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    mDimension = r_process_info[DOMAIN_SIZE];
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeInitialVolumeProcess::Execute() 
{

}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos