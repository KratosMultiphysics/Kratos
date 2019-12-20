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
#include "includes/global_variables.h"

namespace Kratos {

ComputeInitialVolumeProcess::ComputeInitialVolumeProcess(
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

void ComputeInitialVolumeProcess::Execute() 
{
    const std::vector<std::string>& submodel_parts_names = mrModelPart.GetSubModelPartNames();
    for (IndexType i = 0; i < submodel_parts_names.size(); ++i) {
        if (submodel_parts_names[i].substr(0, 8) == mPressureName.substr(0, 8)) { // It is a normal pressure
            auto& r_sub_model_part = mrModelPart.GetSubModelPart(submodel_parts_names[i]);
            const double initial_volume = this->ComputeInitialVolumeSubModel(r_sub_model_part); // Computes the volume inside the SubModel
            this->AssignInitialVolumeToNodes(r_sub_model_part, initial_volume);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

double ComputeInitialVolumeProcess::ComputeInitialVolumeSubModel(
    const ModelPart& rSubModel
)
{
    double max_distance = 0.0, distance;
    const auto& it_node_begin = rSubModel.NodesBegin();
    if (mDimension == 2) {
        // #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rSubModel.Nodes().size()); i++) {
            auto it_node = it_node_begin + i;
            const auto dist_increment = it_node->GetInitialPosition() - 
                it_node_begin->GetInitialPosition();
            distance = inner_prod(dist_increment, dist_increment);
            max_distance = (distance > max_distance) ? distance : max_distance;
        }
        max_distance = std::sqrt(max_distance);
        return Globals::Pi * std::pow(max_distance, 2) * 0.25;
    } else { // 3D
        // TODO
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeInitialVolumeProcess::AssignInitialVolumeToNodes(
    const ModelPart& rSubModel,
    const double InitialVolume
)
{
    const auto& it_node_begin = rSubModel.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rSubModel.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        it_node->SetValue(PRESSURE_INITIAL_VOLUME, InitialVolume);
        it_node->SetValue(PRESSURE_VOLUME, InitialVolume);
    }
}

/***********************************************************************************/
/***********************************************************************************/

int ComputeInitialVolumeProcess::GetPressureIdSubModel(
    const std::string& rSubModelName 
)
{
    const IndexType ref_string_size = (mDimension == 2) ? 18 : 20;
    const IndexType string_size = rSubModelName.size();
    
    if (rSubModelName.size() == ref_string_size) { // from 1 to 9
        return std::stoi(rSubModelName.substr(string_size - 1, string_size));
    } else if (rSubModelName.size() == ref_string_size + 1) { // from 10 to 99
        return std::stoi(rSubModelName.substr(string_size - 2, string_size));
    } else { // from 100 to 999
        return std::stoi(rSubModelName.substr(string_size - 3, string_size));
    }
}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos