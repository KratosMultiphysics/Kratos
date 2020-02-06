//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes

// Project includes
#include "includes/define.h"

// Application incldues
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_find_condition_parent_process.h"

namespace Kratos
{
RansFindConditionParentProcess::RansFindConditionParentProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"                : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"                     : 0
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mModelPartName = mrParameters["model_part_name"].GetString();
    mEchoLevel = mrParameters["echo_level"].GetInt();

    KRATOS_CATCH("");
}

int RansFindConditionParentProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    return 0;

    KRATOS_CATCH("");
}

void RansFindConditionParentProcess::ExecuteInitialize()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
        SetConditionParent(*(r_model_part.ConditionsBegin() + i_condition));

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Found parents for conditions in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansFindConditionParentProcess::Info() const
{
    return std::string("RansFindConditionParentProcess");
}

void RansFindConditionParentProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansFindConditionParentProcess::PrintData(std::ostream& rOStream) const
{
}

void RansFindConditionParentProcess::SetConditionParent(Condition& rCondition)
{
    KRATOS_TRY

    GlobalPointersVector<Element> element_candidates;
    const Condition::GeometryType& r_condition_geometry = rCondition.GetGeometry();
    const int number_of_condition_nodes = r_condition_geometry.PointsNumber();

    std::vector<int> node_ids(number_of_condition_nodes), element_node_ids;

    for (int i_node = 0; i_node < number_of_condition_nodes; ++i_node)
    {
        const NodeType& r_node = r_condition_geometry[i_node];
        const GlobalPointersVector<Element>& r_node_element_candidates =
            r_node.GetValue(NEIGHBOUR_ELEMENTS);

        KRATOS_ERROR_IF(r_node_element_candidates.size() == 0)
            << "No neighbour elements found for node with id=" << r_node.Id()
            << " in condition with condition id=" << rCondition.Id()
            << " belongs to " << mModelPartName << ".\n";

        for (int j_element = 0;
             j_element < static_cast<int>(r_node_element_candidates.size()); ++j_element)
        {
            element_candidates.push_back(r_node_element_candidates(j_element));
        }
        node_ids[i_node] = r_node.Id();
    }

    std::sort(node_ids.begin(), node_ids.end());

    for (int i_element = 0; i_element < static_cast<int>(element_candidates.size()); ++i_element)
    {
        const Element::GeometryType& r_geometry =
            element_candidates[i_element].GetGeometry();
        const int number_of_element_candidate_nodes = r_geometry.PointsNumber();
        if (static_cast<int>(element_node_ids.size()) != number_of_element_candidate_nodes)
            element_node_ids.resize(number_of_element_candidate_nodes);

        for (int i_node = 0; i_node < number_of_element_candidate_nodes; ++i_node)
        {
            element_node_ids[i_node] = r_geometry[i_node].Id();
        }

        std::sort(element_node_ids.begin(), element_node_ids.end());
        if (std::includes(element_node_ids.begin(), element_node_ids.end(),
                          node_ids.begin(), node_ids.end()))
        {
            GlobalPointersVector<Element> parent_element;
            parent_element.clear();
            parent_element.push_back(element_candidates(i_element));
            rCondition.SetValue(NEIGHBOUR_ELEMENTS, parent_element);
            return;
        }
    }

    KRATOS_ERROR << "Parent element for condition id=" << rCondition.Id()
                 << " not found in " << mModelPartName << ".\n";
    KRATOS_CATCH("");
}

} // namespace Kratos.
