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

#include "custom_processes/update_pressure_volume_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/global_variables.h"

namespace Kratos {

UpdatePressureVolumeProcess::UpdatePressureVolumeProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart) 
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    std::size_t dimension = r_process_info[DOMAIN_SIZE];
    mDimension = dimension;    
    mPressureName = (dimension == 2) ? "Normal_Load" : "Pressure_Load";
}

/***********************************************************************************/
/***********************************************************************************/

void UpdatePressureVolumeProcess::Execute() 
{
    const auto it_elem_begin = mrModelPart.ElementsBegin();
    //#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = it_elem_begin + i;

        bool condition_is_active = true;
        if (it_elem->IsDefined(ACTIVE)) {
            condition_is_active = it_elem->Is(ACTIVE);
        }
        if (!condition_is_active) { // if the elem is removed, we add its vol to the nodes of the submodel
            double elem_vol;
            if (mDimension == 2)
                elem_vol = it_elem->GetGeometry().Area();
            else
                elem_vol = it_elem->GetGeometry().Volume();

            const int pressure_id = this->GetPressureId(it_elem);
            if (pressure_id != 0) {
                std::string sub_model_name;
                sub_model_name = mPressureName + "-auto-" + std::to_string(pressure_id);
                auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);
                const auto it_node_begin = r_sub_model_part.NodesBegin();
                #pragma omp parallel for
                for (int j = 0; j < static_cast<int>(r_sub_model_part.Nodes().size()); j++) {
                    auto it_node = it_node_begin + j;
                    const double previous_volume = it_node->GetValue(PRESSURE_VOLUME);
                    const double new_volume = elem_vol + previous_volume;
                    it_node->SetValue(PRESSURE_VOLUME, new_volume);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int UpdatePressureVolumeProcess::GetPressureId(
    ModelPart::ElementIterator itElem
    )
{
    auto& r_geometry = (itElem)->GetGeometry();
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        const int pressure_id = r_geometry[i].GetValue(PRESSURE_ID);
        if (pressure_id != 0) {
            return pressure_id;
        }
    }
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos