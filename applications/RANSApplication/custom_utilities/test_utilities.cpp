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
#include <random>
#include <sstream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "test_utilities.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{
template <>
void AssignRandomValues(
    double& rValue,
    const std::string& rSeed,
    const double MinValue,
    const double MaxValue)
{
    std::seed_seq seed(rSeed.begin(), rSeed.end());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(MinValue, MaxValue);

    rValue = distribution(generator);
}

template <>
void AssignRandomValues(
    array_1d<double, 3>& rValue,
    const std::string& rSeed,
    const double MinValue,
    const double MaxValue)
{
    std::seed_seq seed(rSeed.begin(), rSeed.end());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(MinValue, MaxValue);

    rValue[0] = distribution(generator);
    rValue[1] = distribution(generator);
    rValue[2] = distribution(generator);
}

ModelPart& CreateTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFuncion,
    const std::function<void(ModelPart::NodeType&)>& rAddDofsFunction,
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

    using nid_list = std::vector<ModelPart::IndexType>;

    r_model_part.CreateNewElement(rElementName, 1, nid_list{3, 2, 1}, p_elem_prop);

    r_model_part.CreateNewCondition(rConditionName, 1, nid_list{1, 2}, p_elem_prop);
    r_model_part.CreateNewCondition(rConditionName, 2, nid_list{2, 3}, p_elem_prop);
    r_model_part.CreateNewCondition(rConditionName, 3, nid_list{3, 1}, p_elem_prop);

    r_model_part.Elements().front().Check(r_model_part.GetProcessInfo());
    r_model_part.Conditions().front().Check(r_model_part.GetProcessInfo());

    return r_model_part;
}

ModelPart& CreateScalarVariableTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFuncion,
    const Variable<double>& rDofVariable,
    const int BufferSize,
    const bool DoInitializeElements,
    const bool DoInitializeConditions)
{
    auto& r_model_part = CreateTestModelPart(
        rModel, rElementName, rConditionName, rAddNodalSolutionStepVariablesFuncion,
        [rDofVariable](ModelPart::NodeType& rNode) {
            rNode.AddDof(rDofVariable).SetEquationId(rNode.Id());
        },
        BufferSize);

    if (DoInitializeElements) {
        r_model_part.Elements().front().Initialize();
    }

    if (DoInitializeConditions) {
        r_model_part.Conditions().front().Initialize();
    }

    return r_model_part;
}

template <class TDataType>
void RandomFillNodalHistoricalVariable(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const double MinValue,
    const double MaxValue,
    const int Step)
{
    for (auto& node : rModelPart.Nodes()) {
        std::stringstream seed;
        seed << node.Id() << "_HistoricalV_" << rVariable.Name();
        TDataType& value = node.FastGetSolutionStepValue(rVariable, Step);
        AssignRandomValues<TDataType>(value, seed.str(), MinValue, MaxValue);
    }
}

template <class TContainerType, class TDataType>
void RandomFillContainerVariable(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const double MinValue,
    const double MaxValue)
{
    auto& container = RansCalculationUtilities::GetContainer<TContainerType>(rModelPart);
    for (auto& item : container) {
        std::stringstream seed;
        seed << item.Id() << "_NonHistoricalV_" << rVariable.Name();
        TDataType value = rVariable.Zero();
        AssignRandomValues<TDataType>(value, seed.str(), MinValue, MaxValue);
        item.SetValue(rVariable, value);
    }
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

template void RandomFillNodalHistoricalVariable<double>(
    ModelPart&,
    const Variable<double>&,
    const double,
    const double,
    const int);

template void RandomFillNodalHistoricalVariable<array_1d<double, 3>>(
    ModelPart&,
    const Variable<array_1d<double, 3>>&,
    const double,
    const double,
    const int);

template void RandomFillContainerVariable<ModelPart::NodesContainerType, array_1d<double, 3>>(
    ModelPart&,
    const Variable<array_1d<double, 3>>&,
    const double,
    const double);

template void RandomFillContainerVariable<ModelPart::ConditionsContainerType, array_1d<double, 3>>(
    ModelPart&,
    const Variable<array_1d<double, 3>>&,
    const double,
    const double);

template void RandomFillContainerVariable<ModelPart::ElementsContainerType, array_1d<double, 3>>(
    ModelPart&,
    const Variable<array_1d<double, 3>>&,
    const double,
    const double);

template void RandomFillContainerVariable<ModelPart::NodesContainerType, double>(
    ModelPart&,
    const Variable<double>&,
    const double,
    const double);

template void RandomFillContainerVariable<ModelPart::ConditionsContainerType, double>(
    ModelPart&,
    const Variable<double>&,
    const double,
    const double);

template void RandomFillContainerVariable<ModelPart::ElementsContainerType, double>(
    ModelPart&,
    const Variable<double>&,
    const double,
    const double);

} // namespace RansApplicationTestUtilities
} // namespace Kratos