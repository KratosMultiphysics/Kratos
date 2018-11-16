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
    ModelPart &r_model_part)
    : mr_model_part(r_model_part) 
{
}

/***********************************************************************************/
/***********************************************************************************/

void AssignPressureIdProcess::Execute() 
{
    std::vector<std::string> submodel_parts_names = mr_model_part.GetSubModelPartNames();
    std::vector<std::string> pressure_sub_models;

    for (IndexType i = 0; i < submodel_parts_names.size(); ++i) {
        const IndexType string_size = submodel_parts_names[i].size();
        if (submodel_parts_names[i].substr(0, 11) == "Normal_Load") {

            int pressure_id;
            if (submodel_parts_names[i].size() == 18) { // from 1 to 9
                pressure_id = std::stoi(submodel_parts_names[i].substr(string_size - 1, string_size));
            } else { // from 10 to 99
                pressure_id = std::stoi(submodel_parts_names[i].substr(string_size - 2, string_size));
            }
            this->AssignPressureIdToNodes(submodel_parts_names[i], pressure_id);
        }
    }
}

void AssignPressureIdProcess::AssignPressureIdToNodes(std::string rSubModelPartName, const int PressureId)
{
    // Loop over nodes of that submodel to assign prressure id
    for (ModelPart::NodeIterator it = mr_model_part.GetSubModelPart(rSubModelPartName).NodesBegin();
         it != mr_model_part.GetSubModelPart(rSubModelPartName).NodesEnd();
         ++it) {

		(*it).SetValue(PRESSURE_ID, PressureId);
    }
}

}  // namespace Kratos