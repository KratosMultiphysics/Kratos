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
    ModelPart& rDemModelPart)
    : mrModelPart(rModelPart),
      mrDEMModelPart(rDemModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

void RemoveAloneDEMElementsProcess::Execute() 
{
    const auto it_node_begin = mrModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        it_node->SetValue(NUMBER_OF_ACTIVE_ELEMENTS, 0);
    }

    const auto it_elem_begin = mrModelPart.ElementsBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = it_elem_begin + i;
        auto& r_geometry = it_elem->GetGeometry();
        for (SizeType i_node = 0; i_node < r_geometry.PointsNumber(); i_node++) {
            auto& r_node = r_geometry[i];
            int& r_number_of_active_elements = r_node.GetValue(NUMBER_OF_ACTIVE_ELEMENTS);
            #pragma omp atomic
            r_number_of_active_elements++;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        const int& r_number_of_active_elements = it_node->GetValue(NUMBER_OF_ACTIVE_ELEMENTS);
        if (r_number_of_active_elements == 0) {
            auto& p_associated_DEM = it_node->GetValue(DEM_PARTICLE_POINTER);
            const int DEM_id = p_associated_DEM->Id();
            mrDEMModelPart.GetNode(DEM_id).Set(TO_ERASE, true);
        }
    }
    mrDEMModelPart.RemoveElementsFromAllLevels(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos