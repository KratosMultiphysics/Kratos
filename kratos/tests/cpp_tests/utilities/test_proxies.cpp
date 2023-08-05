//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Máté Kelemen
//

// Project includes
#include "includes/global_variables.h"
#include "includes/variables.h"
#include "testing/testing.h"
#include "utilities/proxies.h"
#include "includes/model_part_io.h"
#include "containers/model.h"


namespace Kratos::Testing {


auto MakeProxyTestModel()
{
    auto p_model = std::make_unique<Model>();
    auto& r_model_part = p_model->CreateModelPart("root");
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    auto p_properties = Properties::Pointer(new Properties());

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.GetSolutionStepValue(PRESSURE) = r_node.Id();
        r_node.GetValue(PRESSURE) = r_node.Id();
    }

    r_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);
    r_model_part.CreateNewElement("Element2D3N", 2, {2, 3, 1}, p_properties);
    for (auto& r_element : r_model_part.Elements()) {
        r_element.GetValue(PRESSURE) = r_element.Id();
    }

    r_model_part.CreateNewCondition("LineCondition2D2N", 1, std::vector<ModelPart::IndexType> {1, 2}, p_properties);
    r_model_part.CreateNewCondition("LineCondition2D2N", 2, std::vector<ModelPart::IndexType> {3, 4}, p_properties);
    for (auto& r_condition : r_model_part.Conditions()) {
        r_condition.GetValue(PRESSURE) = r_condition.Id();
    }

    return p_model;
}


template <Globals::DataLocation TLocation, class TContainer>
void TestEntityProxy(TContainer& rEntities)
{
    KRATOS_CHECK_IS_FALSE(rEntities.empty());
    const TContainer& r_immutable_entities = rEntities;

    // Check immutable EntityProxy
    for (const auto& r_entity : r_immutable_entities) {
        const auto id = r_entity.Id();
        auto proxy = MakeProxy<TLocation>(r_entity);

        // Check EntityProxy::GetValue
        KRATOS_CHECK_EQUAL(proxy.GetValue(PRESSURE), double(id));
    }

    // Check mutable EntityProxy
    for (auto& r_entity : rEntities) {
        const auto id = r_entity.Id();
        auto proxy = MakeProxy<TLocation>(r_entity);

        // Check EntityProxy::GetValue
        KRATOS_CHECK_EQUAL(proxy.GetValue(PRESSURE), double(id));

        // Check mutable EntityProxy::GetValue
        proxy.GetValue(PRESSURE) *= 2;
        KRATOS_CHECK_EQUAL(proxy.GetValue(PRESSURE), 2 * double(id));

        // Check EntityProxy::SetValue
        proxy.SetValue(PRESSURE, proxy.GetValue(PRESSURE) * 2);
        KRATOS_CHECK_EQUAL(proxy.GetValue(PRESSURE), 2 * double(2 * double(id)));
    }
}


KRATOS_TEST_CASE_IN_SUITE(HistoricalNodeProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    TestEntityProxy<Globals::DataLocation::NodeHistorical>(r_model_part.Nodes());
}


KRATOS_TEST_CASE_IN_SUITE(NonHistoricalNodeProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    TestEntityProxy<Globals::DataLocation::NodeNonHistorical>(r_model_part.Nodes());
}


KRATOS_TEST_CASE_IN_SUITE(ElementProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    TestEntityProxy<Globals::DataLocation::Element>(r_model_part.Elements());
}


KRATOS_TEST_CASE_IN_SUITE(ConditionProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    TestEntityProxy<Globals::DataLocation::Condition>(r_model_part.Conditions());
}


} // namespace Kratos::Testing
