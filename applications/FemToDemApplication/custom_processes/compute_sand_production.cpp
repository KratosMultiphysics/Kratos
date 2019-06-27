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

#include "custom_processes/compute_sand_production.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {


ComputeSandProduction::ComputeSandProduction(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart) 
{
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeSandProduction::Execute() 
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    double erased_volume = r_process_info[ERASED_VOLUME];

    const auto it_element_begin = mrModelPart.ElementsBegin();
    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {

        auto it_elem = it_element_begin + i;
        bool is_active = true;
        if (it_elem->IsDefined(ACTIVE))
            is_active = it_elem->Is(ACTIVE);

        if (!is_active) { // it is going to be removed inside GenerateDEM
            erased_volume += it_elem->GetGeometry().Volume();
        }
    }
    r_process_info[ERASED_VOLUME] = erased_volume;
}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos