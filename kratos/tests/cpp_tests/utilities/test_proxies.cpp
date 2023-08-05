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
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    auto p_properties = Properties::Pointer(new Properties());

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        const auto id = r_node.Id();
        r_node.GetSolutionStepValue(DISPLACEMENT) = array_1d<double,3> {0.0, double(id), double(2 * id)};
        r_node.GetValue(PRESSURE) = double(3 * id);
    }

    r_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);
    r_model_part.CreateNewElement("Element2D3N", 2, {2, 3, 1}, p_properties);
    for (auto& r_element : r_model_part.Elements()) {
        r_element.GetValue(PRESSURE) = double(3 * r_element.Id());
    }

    r_model_part.CreateNewCondition("LineCondition2D2N", 1, std::vector<ModelPart::IndexType> {1, 2}, p_properties);
    r_model_part.CreateNewCondition("LineCondition2D2N", 2, std::vector<ModelPart::IndexType> {3, 4}, p_properties);
    for (auto& r_condition : r_model_part.Conditions()) {
        r_condition.GetValue(PRESSURE) = double(3 * r_condition.Id());
    }

    return p_model;
}


KRATOS_TEST_CASE_IN_SUITE(HistoricalNodeProxy, KratosCoreFastSuite)
{
    auto p_model = MakeProxyTestModel();
    ModelPart& r_model_part = p_model->GetModelPart("root");
    const ModelPart& r_const_model_part = p_model->GetModelPart("root");

    // Check immutable EntityProxy
    for (const auto& r_node : r_const_model_part.Nodes()) {
        const auto id = r_node.Id();
        auto proxy = MakeProxy<Globals::DataLocation::NodeHistorical>(r_node);

        // Check EntityProxy::GetValue
        {
            const array_1d<double,3> reference {0.0, double(id), 2 * double(id)};
            KRATOS_CHECK_VECTOR_EQUAL(r_node.GetSolutionStepValue(DISPLACEMENT), reference);
            KRATOS_CHECK_VECTOR_EQUAL(proxy.GetValue(DISPLACEMENT), reference);
        }
    }

    // Check mutable EntityProxy
    for (auto& r_node : r_model_part.Nodes()) {
        const auto id = r_node.Id();
        auto proxy = MakeProxy<Globals::DataLocation::NodeHistorical>(r_node);

        // Check EntityProxy::GetValue
        {
            const array_1d<double,3> reference {0.0, double(id), 2 * double(id)};
            KRATOS_CHECK_VECTOR_EQUAL(proxy.GetValue(DISPLACEMENT), reference);
        }

        // Check mutable EntityProxy::GetValue
        {
            proxy.GetValue(DISPLACEMENT) *= 2;
            const array_1d<double,3> reference {0.0, 2 * double(id), 4 * double(id)};
            KRATOS_CHECK_VECTOR_EQUAL(proxy.GetValue(DISPLACEMENT), reference);
        }

        // Check EntityProxy::SetValue
        {
            proxy.SetValue(DISPLACEMENT, proxy.GetValue(DISPLACEMENT) * 2);
            const array_1d<double,3> reference {0.0, 4 * double(id), 8 * double(id)};
            KRATOS_CHECK_VECTOR_EQUAL(proxy.GetValue(DISPLACEMENT), reference);
        }
    }
}


} // namespace Kratos::Testing
