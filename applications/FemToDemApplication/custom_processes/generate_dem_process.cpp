//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/generate_dem_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {


GenerateDemProcess::GenerateDemProcess(
    ModelPart& r_model_part,
    ModelPart& r_dem_model_part)
    : mrModelPart(r_model_part),
      mrDEMModelPart(r_dem_model_part)
{
}

/***********************************************************************************/
/***********************************************************************************/

void GenerateDemProcess::Execute() 
{

}


}  // namespace Kratos