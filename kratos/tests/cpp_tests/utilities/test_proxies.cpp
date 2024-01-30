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


namespace {


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
    r_model_part.CreateNewElement("Element2D3N", 2, {3, 4, 1}, p_properties);
    for (auto& r_element : r_model_part.Elements()) {
        r_element.GetValue(PRESSURE) = r_element.Id();
    }

    r_model_part.CreateNewCondition("LineCondition2D2N", 1, std::vector<ModelPart::IndexType> {1, 2}, p_properties);
    r_model_part.CreateNewCondition("LineCondition2D2N", 2, std::vector<ModelPart::IndexType> {3, 4}, p_properties);
    for (auto& r_condition : r_model_part.Conditions()) {
        r_condition.GetValue(PRESSURE) = r_condition.Id();
    }

    r_model_part.GetProcessInfo().GetValue(PRESSURE) = 2.0;
    r_model_part.GetValue(PRESSURE) = 2.0;

    return p_model;
}


/// Compile-time check whether EntityProxy::GetValue returns a mutable reference
template <class TEntityProxy,
          std::enable_if_t<std::is_same_v<decltype(std::declval<TEntityProxy>().GetValue(PRESSURE)),double&>,bool> = true>
bool IsEntityProxyMutable(TEntityProxy Proxy)
{
    return true;
}


/// Compile-time check whether EntityProxy::GetValue returns a mutable reference
template <class TEntityProxy,
          std::enable_if_t<
            std::is_same_v<decltype(std::declval<TEntityProxy>().GetValue(PRESSURE)),double>
            || std::is_same_v<decltype(std::declval<TEntityProxy>().GetValue(PRESSURE)),const double&>,bool> = true>
bool IsEntityProxyMutable(TEntityProxy Proxy)
{
    return false;
}


template <Globals::DataLocation TLocation, class TContainer>
void TestEntityProxy(TContainer& rEntities)
{
    KRATOS_EXPECT_FALSE(rEntities.empty());
    const TContainer& r_immutable_entities = rEntities;

    // Check immutable EntityProxy
    for (const auto& r_entity : r_immutable_entities) {
        const auto id = r_entity.Id();
        auto proxy = MakeProxy<TLocation>(r_entity);

        // Check const-correctness
        KRATOS_EXPECT_FALSE(IsEntityProxyMutable(proxy));

        // Check EntityProxy::HasValue
        KRATOS_EXPECT_TRUE(proxy.HasValue(PRESSURE));
        KRATOS_EXPECT_FALSE(proxy.HasValue(VELOCITY));

        // Check EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), double(id));
    }

    // Check mutable EntityProxy
    for (auto& r_entity : rEntities) {
        const auto id = r_entity.Id();
        auto proxy = MakeProxy<TLocation>(r_entity);

        // Check const-correctness
        KRATOS_EXPECT_TRUE(IsEntityProxyMutable(proxy));

        // Check EntityProxy::HasValue
        KRATOS_EXPECT_TRUE(proxy.HasValue(PRESSURE));
        KRATOS_EXPECT_FALSE(proxy.HasValue(VELOCITY));

        // Check EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), double(id));

        // Check mutable EntityProxy::GetValue
        proxy.GetValue(PRESSURE) *= 2;
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2 * double(id));

        // Check EntityProxy::SetValue
        proxy.SetValue(PRESSURE, proxy.GetValue(PRESSURE) * 2);
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2 * double(2 * double(id)));
    }
}


} // unnamed namespace


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


KRATOS_TEST_CASE_IN_SUITE(ProcessInfoProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();

    ModelPart& r_mutable_model_part = p_model->GetModelPart("root");
    const ModelPart& r_immutable_model_part = r_mutable_model_part;

    {
        auto proxy = MakeProxy<Globals::DataLocation::ProcessInfo>(r_immutable_model_part.GetProcessInfo());

        // Check immutable EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);
    }

    {
        auto proxy = MakeProxy<Globals::DataLocation::ProcessInfo>(r_mutable_model_part.GetProcessInfo());

        // Check immutable EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);

        // Check mutable EntityProxy::GetValue
        proxy.GetValue(PRESSURE) *= 3.0;
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 3.0 * 2.0);

        // Check EntityProxy::SetValue
        proxy.SetValue(PRESSURE, 2.0);
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);
    }

    {
        auto proxy = MakeProxy<Globals::DataLocation::ProcessInfo>(r_immutable_model_part);

        // Check immutable EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);
    }

    {
        auto proxy = MakeProxy<Globals::DataLocation::ProcessInfo>(r_mutable_model_part);

        // Check immutable EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);

        // Check mutable EntityProxy::GetValue
        proxy.GetValue(PRESSURE) *= 3.0;
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 3.0 * 2.0);

        // Check EntityProxy::SetValue
        proxy.SetValue(PRESSURE, 2.0);
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);
    }
}


KRATOS_TEST_CASE_IN_SUITE(ModelPartProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();

    ModelPart& r_mutable_model_part = p_model->GetModelPart("root");
    const ModelPart& r_immutable_model_part = r_mutable_model_part;

    {
        auto proxy = MakeProxy<Globals::DataLocation::ModelPart>(r_immutable_model_part);

        // Check immutable EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);
    }

    {
        auto proxy = MakeProxy<Globals::DataLocation::ModelPart>(r_mutable_model_part);

        // Check immutable EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);

        // Check mutable EntityProxy::GetValue
        proxy.GetValue(PRESSURE) *= 3.0;
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 3.0 * 2.0);

        // Check EntityProxy::SetValue
        proxy.SetValue(PRESSURE, 2.0);
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2.0);
    }
}


template <class TMutableContainerProxy, class TImmutableContainerProxy>
void TestContainerProxy(TMutableContainerProxy MutableProxies, TImmutableContainerProxy ImmutableProxies)
{
    // Check ContainerProxy::empty
    KRATOS_EXPECT_FALSE(MutableProxies.empty());
    KRATOS_EXPECT_FALSE(ImmutableProxies.empty());

    // Check ContainerProxy::size
    KRATOS_EXPECT_EQ(MutableProxies.size(), ImmutableProxies.size());

    // Make sure that the two container proxies point to the same range
    KRATOS_EXPECT_EQ(&(*MutableProxies.begin()).GetEntity(), &(*ImmutableProxies.begin()).GetEntity());

    // Check immutable ContainerProxy
    for (auto proxy : ImmutableProxies) {
        const auto id = proxy.GetEntity().Id();

        // Check const-correctness
        KRATOS_EXPECT_FALSE(IsEntityProxyMutable(proxy));

        // Check EntityProxy::HasValue
        KRATOS_EXPECT_TRUE(proxy.HasValue(PRESSURE));
        KRATOS_EXPECT_FALSE(proxy.HasValue(VELOCITY));

        // Check EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), double(id));
    }

    // Check mutable ContainerProxy
    for (auto proxy : MutableProxies) {
        const auto id = proxy.GetEntity().Id();

        // Check const-correctness
        KRATOS_EXPECT_TRUE(IsEntityProxyMutable(proxy));

        // Check EntityProxy::HasValue
        KRATOS_EXPECT_TRUE(proxy.HasValue(PRESSURE));
        KRATOS_EXPECT_FALSE(proxy.HasValue(VELOCITY));

        // Check EntityProxy::GetValue
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), double(id));

        // Check mutable EntityProxy::GetValue
        proxy.GetValue(PRESSURE) *= 2;
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2 * double(id));

        // Check EntityProxy::SetValue
        proxy.SetValue(PRESSURE, proxy.GetValue(PRESSURE) * 2);
        KRATOS_EXPECT_EQ(proxy.GetValue(PRESSURE), 2 * double(2 * double(id)));
    }
}


KRATOS_TEST_CASE_IN_SUITE(HistoricalNodeContainerProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    const ModelPart& r_immutable_model_part = r_model_part;
    TestContainerProxy(MakeProxy<Globals::DataLocation::NodeHistorical>(r_model_part),
                       MakeProxy<Globals::DataLocation::NodeHistorical>(r_immutable_model_part));
}


KRATOS_TEST_CASE_IN_SUITE(NonHistoricalNodeContainerProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    const ModelPart& r_immutable_model_part = r_model_part;
    TestContainerProxy(MakeProxy<Globals::DataLocation::NodeNonHistorical>(r_model_part),
                       MakeProxy<Globals::DataLocation::NodeNonHistorical>(r_immutable_model_part));
}


KRATOS_TEST_CASE_IN_SUITE(ElementContainerProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    const ModelPart& r_immutable_model_part = r_model_part;
    TestContainerProxy(MakeProxy<Globals::DataLocation::Element>(r_model_part),
                       MakeProxy<Globals::DataLocation::Element>(r_immutable_model_part));
}


KRATOS_TEST_CASE_IN_SUITE(ConditionContainerProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    const ModelPart& r_immutable_model_part = r_model_part;
    TestContainerProxy(MakeProxy<Globals::DataLocation::Condition>(r_model_part),
                       MakeProxy<Globals::DataLocation::Condition>(r_immutable_model_part));
}


} // namespace Kratos::Testing
