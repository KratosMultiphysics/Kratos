//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <functional>
#include <sstream>

// External includes

// Project includes
#include "containers/model.h"
#include "containers/global_pointers_vector.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"

// Include base h
#include "test_utilities.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{

ModelPart& CreateTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFuncion,
    const std::function<void(ModelPart::NodeType&)>& rAddDofsFunction,
    const std::function<void(Properties&)>& rSetProperties,
    const int BufferSize)
{
    auto& r_model_part = rModel.CreateModelPart("test", BufferSize);
    rAddNodalSolutionStepVariablesFuncion(r_model_part);

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        rAddDofsFunction(r_node);
    }

    Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("RansNewtonian2DLaw").Clone());
    rSetProperties(*p_elem_prop);

    using nid_list = std::vector<ModelPart::IndexType>;

    r_model_part.CreateNewElement(rElementName, 1, nid_list{3, 2, 1}, p_elem_prop);
    auto& r_element = r_model_part.Elements().front();

    r_model_part.CreateNewCondition(rConditionName, 1, nid_list{1, 2}, p_elem_prop)->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>{&r_element});
    r_model_part.CreateNewCondition(rConditionName, 2, nid_list{2, 3}, p_elem_prop)->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>{&r_element});
    r_model_part.CreateNewCondition(rConditionName, 3, nid_list{3, 1}, p_elem_prop)->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>{&r_element});

    RansVariableUtilities::SetElementConstitutiveLaws(r_model_part.Elements());

    return r_model_part;
}

ModelPart& CreateScalarVariableTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFuncion,\
    const std::function<void(Properties&)>& rSetProperties,
    const Variable<double>& rDofVariable,
    const int BufferSize,
    const bool DoInitializeElements,
    const bool DoInitializeConditions)
{
    auto& r_model_part = CreateTestModelPart(
        rModel, rElementName, rConditionName, rAddNodalSolutionStepVariablesFuncion,
        [&rDofVariable](ModelPart::NodeType& rNode) {
            rNode.AddDof(rDofVariable).SetEquationId(rNode.Id());
        },
        rSetProperties,
        BufferSize);

    if (DoInitializeElements) {
        r_model_part.Elements().front().Initialize(r_model_part.GetProcessInfo());
    }

    if (DoInitializeConditions) {
        r_model_part.Conditions().front().Initialize(r_model_part.GetProcessInfo());
    }

    return r_model_part;
}

template <class TContainerType>
void TestEquationIdVector(
    ModelPart& rModelPart)
{
    auto eqn_ids = std::vector<IndexType>{};
    for (const auto& r_item : RansCalculationUtilities::GetContainer<TContainerType>(rModelPart)) {
        r_item.EquationIdVector(eqn_ids, rModelPart.GetProcessInfo());
        KRATOS_CHECK_EQUAL(eqn_ids.size(), r_item.GetGeometry().PointsNumber());
        const auto& r_geometry = r_item.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            KRATOS_CHECK_EQUAL(eqn_ids[i_node], r_geometry[i_node].Id());
        }
    }
}

template <class TContainerType>
void TestGetDofList(
    ModelPart& rModelPart,
    const Variable<double>& rVariable)
{
    auto dofs = Element::DofsVectorType{};
    for (const auto& r_item : RansCalculationUtilities::GetContainer<TContainerType>(rModelPart)) {
        r_item.GetDofList(dofs, rModelPart.GetProcessInfo());
        KRATOS_CHECK_EQUAL(dofs.size(), r_item.GetGeometry().PointsNumber());
        const auto& r_geometry = r_item.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            KRATOS_CHECK_EQUAL(dofs[i_node]->GetVariable(), rVariable);
            KRATOS_CHECK_EQUAL(dofs[i_node]->EquationId(), r_geometry[i_node].Id());
        }
    }
}

void CheckElementsAndConditions(
    const ModelPart& rModelPart)
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    for (const auto& r_element : rModelPart.Elements()) {
        r_element.Check(r_process_info);
    }

    for (const auto& r_condition : rModelPart.Conditions()) {
        r_condition.Check(r_process_info);
    }
}

// template instantiations
template void TestEquationIdVector<ModelPart::ElementsContainerType>(
    ModelPart&);
template void TestEquationIdVector<ModelPart::ConditionsContainerType>(
    ModelPart&);

template void TestGetDofList<ModelPart::ElementsContainerType>(
    ModelPart&,
    const Variable<double>&);

template void TestGetDofList<ModelPart::ConditionsContainerType>(
    ModelPart&,
    const Variable<double>&);

} // namespace RansApplicationTestUtilities
} // namespace Kratos