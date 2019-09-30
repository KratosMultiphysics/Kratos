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

#include "custom_processes/assign_pressure_id_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {


AssignPressureIdProcess::AssignPressureIdProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart) 
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
    auto& r_process_info = mrModelPart.GetProcessInfo();
    const std::size_t dimension = r_process_info[DOMAIN_SIZE];
    const IndexType ref_string_size = (dimension == 2) ? 18 : 20;

    for (IndexType i = 0; i < submodel_parts_names.size(); ++i) {
        const IndexType string_size = submodel_parts_names[i].size();
        if (submodel_parts_names[i].substr(0, 11) == mPressureName.substr(0, 11)) {
            int pressure_id;
            if (submodel_parts_names[i].size() == ref_string_size) { // from 1 to 9
                pressure_id = std::stoi(submodel_parts_names[i].substr(string_size - 1, string_size));
            } else { // from 10 to 99
                pressure_id = std::stoi(submodel_parts_names[i].substr(string_size - 2, string_size));
            }
            this->AssignPressureIdToNodes(submodel_parts_names[i], pressure_id);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AssignPressureIdProcess::AssignPressureIdToNodes(
    std::string rSubModelPartName, 
    const int PressureId
    )
{
    auto& r_submodel = mrModelPart.GetSubModelPart(rSubModelPartName);
    const auto it_node_begin = r_submodel.NodesBegin();
    // Loop over nodes of that submodel to assign prressure id
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_submodel.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        it_node->SetValue(PRESSURE_ID, PressureId);
    }
}

}  // namespace Kratos