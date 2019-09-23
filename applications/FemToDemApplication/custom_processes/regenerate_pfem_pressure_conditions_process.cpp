//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:             BSD License
//                      Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/regenerate_pfem_pressure_conditions_process.h"
// #include "processes/find_elements_neighbours_process.h"

namespace Kratos {

template <SizeType TDim>
RegeneratePfemPressureConditionsProcess<TDim>::RegeneratePfemPressureConditionsProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void RegeneratePfemPressureConditionsProcess<TDim>::Execute()
{
    // We search the neighbours for the generation of line loads
    auto find_neigh = FindElementalNeighboursProcess(mrModelPart, TDim, 5);
    find_neigh.Execute();
    auto& r_process_info = mrModelPart.GetProcessInfo();

    // Remove previous line loads-> Only the 1st iteration
    if (r_process_info[PFEM_PRESSURE_ITERATION] == 1) {
        // We fill the properties vectors to be reassigned afterwards
        this->SavePreviousProperties();
        this->RemovePreviousLineLoads();
        this->ResetFlagOnElements();
    }
    // Generate the new ones
    this->CreateNewConditions();
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void RegeneratePfemPressureConditionsProcess<TDim>::GetMaximumConditionIdOnSubmodelPart(
    int& rMaximumConditionId
)
{
    rMaximumConditionId = 0;
    for (auto it_cond = mrModelPart.ConditionsBegin(); it_cond != mrModelPart.ConditionsEnd(); it_cond++) {
        rMaximumConditionId = (((*it_cond)).Id() > rMaximumConditionId) ? ((*it_cond)).Id() : rMaximumConditionId;
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<3>::CreatePressureLoads(
    const int Id1,
    const int Id2,
    const int Id3,
	ModelPart::ElementsContainerType::ptr_iterator itElem,
	ModelPart& rSubModelPart,
    ModelPart::PropertiesType::Pointer pProperties,
    int& rMaximumConditionId
    )
{
    auto& r_geom = (*itElem)->GetGeometry();
    std::vector<IndexType> condition_nodes_id(3);
    condition_nodes_id[0] = r_geom[Id1].Id();
    condition_nodes_id[1] = r_geom[Id2].Id();
    condition_nodes_id[2] = r_geom[Id3].Id();
    rMaximumConditionId++;

    // Adding the nodes to the SubModelPart
	rSubModelPart.AddNodes(condition_nodes_id);

    // We create the Line Load Condition
    const auto p_pressure_condition = rSubModelPart.CreateNewCondition(
                                        "SurfaceLoadCondition3D3N",
                                        rMaximumConditionId,
                                        condition_nodes_id,
                                        pProperties, 0);

    // Adding the conditions to the computing model part
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(p_pressure_condition); 
}

/***********************************************************************************/
/***********************************************************************************/

template class RegeneratePfemPressureConditionsProcess<3>;

}  // namespace Kratos