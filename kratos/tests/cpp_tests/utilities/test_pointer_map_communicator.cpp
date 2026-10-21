//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ihar Antonau
//

// System includes
#include <memory>

// Project includes
#include "containers/model.h"
#include "includes/global_pointer.h"
#include "includes/variables.h"
#include "testing/testing.h"
#include "utilities/pointer_map_communicator.h"

namespace Kratos::Testing {

// SensitivityBuilder passes a temporary lambda to GetApplyProxy. The proxy must keep its own copy (seen through
// the captured pointer's use count), not a reference to the temporary, which dies at the end of that statement.
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerMapCommunicatorApplyProxyOwnsTemporaryFunctor, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    auto p_node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    p_node->SetValue(TEMPERATURE, 1.0);

    GlobalPointerMapCommunicator<Node, double> pointer_map_communicator(r_model_part.GetCommunicator().GetDataCommunicator());

    auto p_scale = std::make_shared<double>(2.0);
    auto apply_proxy = pointer_map_communicator.GetApplyProxy([p_scale](Node& rNode, const double& rValue) {
        rNode.GetValue(TEMPERATURE) += rValue * (*p_scale);
    });
    KRATOS_EXPECT_EQ(p_scale.use_count(), 2);

    GlobalPointer<Node> gp(p_node.get());
    apply_proxy.Assign(gp, 3.0);
    apply_proxy.SendAndApplyRemotely();
    KRATOS_EXPECT_DOUBLE_EQ(p_node->GetValue(TEMPERATURE), 7.0);
}

// A temporary functor is moved into the proxy, so a move-only functor (owning a std::unique_ptr) is accepted.
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerMapCommunicatorApplyProxyMovesTemporaryFunctor, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    auto p_node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    p_node->SetValue(TEMPERATURE, 1.0);

    GlobalPointerMapCommunicator<Node, double> pointer_map_communicator(r_model_part.GetCommunicator().GetDataCommunicator());

    auto apply_proxy = pointer_map_communicator.GetApplyProxy(
        [p_scale = std::make_unique<double>(2.0)](Node& rNode, const double& rValue) {
            rNode.GetValue(TEMPERATURE) += rValue * (*p_scale);
        });

    GlobalPointer<Node> gp(p_node.get());
    apply_proxy.Assign(gp, 3.0);
    apply_proxy.SendAndApplyRemotely();
    KRATOS_EXPECT_DOUBLE_EQ(p_node->GetValue(TEMPERATURE), 7.0);
}

// An lvalue functor is referenced, not copied: the use count stays at the named functor's single copy.
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerMapCommunicatorApplyProxyReferencesLvalueFunctor, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    auto p_node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    p_node->SetValue(TEMPERATURE, 1.0);

    GlobalPointerMapCommunicator<Node, double> pointer_map_communicator(r_model_part.GetCommunicator().GetDataCommunicator());

    auto p_scale = std::make_shared<double>(2.0);
    const auto apply_functor = [p_scale](Node& rNode, const double& rValue) {
        rNode.GetValue(TEMPERATURE) += rValue * (*p_scale);
    };
    auto apply_proxy = pointer_map_communicator.GetApplyProxy(apply_functor);
    KRATOS_EXPECT_EQ(p_scale.use_count(), 2);

    GlobalPointer<Node> gp(p_node.get());
    apply_proxy.Assign(gp, 3.0);
    apply_proxy.SendAndApplyRemotely();
    KRATOS_EXPECT_DOUBLE_EQ(p_node->GetValue(TEMPERATURE), 7.0);
}

} // namespace Kratos::Testing
