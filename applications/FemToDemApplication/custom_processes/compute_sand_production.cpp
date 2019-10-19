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
        std::vector<double> damage;  
        it_elem->GetValueOnIntegrationPoints(DAMAGE_ELEMENT, damage, r_process_info);
        if (damage[0] >= 0.7 && !it_elem->GetValue(VOLUME_COUNTED)) {
            erased_volume += it_elem->GetGeometry().Volume();
            it_elem->SetValue(VOLUME_COUNTED, true);
        }
    }
    r_process_info[ERASED_VOLUME] = erased_volume;
}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos