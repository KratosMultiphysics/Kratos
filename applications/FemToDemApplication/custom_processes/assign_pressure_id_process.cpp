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

#include "custom_processes/assign_pressure_id_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {


AssignPressureIdProcess::AssignPressureIdProcess(
    ModelPart& r_model_part)
    : mrModelPart(r_model_part) 
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    const std::size_t dimension = r_process_info[DOMAIN_SIZE];
    mPressureName = (dimension == 2) ? "Normal_Load" : "Pressure_Load";
}

/***********************************************************************************/
/***********************************************************************************/

void AssignPressureIdProcess::Execute() 
{
    std::vector<std::string> submodel_parts_names = mrModelPart.GetSubModelPartNames();
    std::vector<std::string> pressure_sub_models;

    auto& r_process_info = mrModelPart.GetProcessInfo();
    const std::size_t dimension = r_process_info[DOMAIN_SIZE];

    for (IndexType i = 0; i < submodel_parts_names.size(); ++i) {
        const IndexType string_size = submodel_parts_names[i].size();
        KRATOS_WATCH(submodel_parts_names[i].substr(0, 11))
        KRATOS_WATCH(mPressureName.substr(0, 11))
        if (submodel_parts_names[i].substr(0, 11) == mPressureName.substr(0, 11)) {

            int pressure_id;

            if (dimension == 2) { // 2D
                if (submodel_parts_names[i].size() == 18) { // from 1 to 9
                    pressure_id = std::stoi(submodel_parts_names[i].substr(string_size - 1, string_size));
                } else { // from 10 to 99
                    pressure_id = std::stoi(submodel_parts_names[i].substr(string_size - 2, string_size));
                }            
            } else { // 3D

            }

            this->AssignPressureIdToNodes(submodel_parts_names[i], pressure_id);
        }
    }
}

void AssignPressureIdProcess::AssignPressureIdToNodes(
    std::string rSubModelPartName, 
    const int PressureId
    )
{
    // Loop over nodes of that submodel to assign prressure id
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.GetSubModelPart(rSubModelPartName).Nodes().size()); i++) {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->SetValue(PRESSURE_ID, PressureId);
    }
}

}  // namespace Kratos