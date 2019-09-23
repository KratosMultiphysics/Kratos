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
    int max_id;
    this->GetMaximumConditionIdOnSubmodelPart(max_id);
    const auto it_cond_begin = mrModelPart.ConditionsBegin();
    ModelPart::PropertiesType::Pointer pProperties = it_cond_begin->pGetProperties();

    std::vector<IndexType> condition_nodes_id(3);
    auto &elem = mrModelPart.GetElement(194);
    condition_nodes_id[0] = elem.GetGeometry()[0].Id();
    condition_nodes_id[1] = elem.GetGeometry()[1].Id();
    condition_nodes_id[2] = elem.GetGeometry()[2].Id();
    max_id++;
    max_id = 1000000;
    mrModelPart.CreateSubModelPart("dummy_pressure");
    auto& r_sub_model = mrModelPart.GetSubModelPart("dummy_pressure");

    r_sub_model.AddNodes(condition_nodes_id);
    // We create the Line Load Condition
    const auto p_pressure_condition = r_sub_model.CreateNewCondition(
                                        "SurfaceLoadCondition3D3N",
                                        max_id,
                                        condition_nodes_id,
                                        pProperties, 0);
    p_pressure_condition->SetValue(POSITIVE_FACE_PRESSURE, -1000);
    // Adding the conditions to the computing model part
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(p_pressure_condition);
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

template class RegeneratePfemPressureConditionsProcess<3>;

}  // namespace Kratos