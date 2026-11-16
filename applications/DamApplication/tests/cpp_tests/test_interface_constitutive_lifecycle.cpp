//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers

// The hydro-mechanical uplift coupling solves the mechanical problem several
// times inside one time step. Each solve evaluates the interface constitutive
// law at the current trial state. This test verifies that those trial
// evaluations do not commit constitutive history: only FinalizeSolutionStep
// does, exactly once.

#include "dam_fast_suite.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/variables.h"

#include "custom_elements/small_displacement_interface_element.hpp"
#include "custom_constitutive/interface_element_laws/bilinear_cohesive_2D_law.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

constexpr double critical_displacement = 0.01;
constexpr double damage_threshold = 0.3;

struct JointElement
{
    ModelPart* p_model_part;
    Element::Pointer p_element;
    std::vector<Node::Pointer> nodes;
};

JointElement CreateJoint(Model& rModel,
                         const std::string& rModelPartName,
                         const int ElementId,
                         const int PropertiesId)
{
    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_JOINT_WIDTH);
    r_model_part.AddNodalSolutionStepVariable(NODAL_JOINT_AREA);
    r_model_part.GetProcessInfo()[IS_CONVERGED] = true;

    auto p_properties = r_model_part.CreateNewProperties(PropertiesId);
    p_properties->SetValue(CONSTITUTIVE_LAW, ConstitutiveLaw::Pointer(new BilinearCohesive2DLaw()));
    p_properties->SetValue(INITIAL_JOINT_WIDTH, 0.0);
    p_properties->SetValue(THICKNESS, 1.0);
    p_properties->SetValue(CRITICAL_DISPLACEMENT, critical_displacement);
    p_properties->SetValue(DAMAGE_THRESHOLD, damage_threshold);
    p_properties->SetValue(YIELD_STRESS, 1000.0);
    p_properties->SetValue(YOUNG_MODULUS, 1.0e6);
    p_properties->SetValue(FRICTION_COEFFICIENT, 0.0);

    std::vector<Node::Pointer> nodes;
    const double coordinates[4][2] = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}};
    for (int i = 0; i < 4; ++i)
    {
        Node::Pointer p_node = r_model_part.CreateNewNode(
            10 * ElementId + i, coordinates[i][0], coordinates[i][1], 0.0);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        const Matrix zero_stress = ZeroMatrix(2, 2);
        p_node->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = zero_stress;
        nodes.push_back(p_node);
    }

    Element::Pointer p_element = r_model_part.CreateNewElement(
        "SmallDisplacementInterfaceElement2D4N",
        ElementId,
        std::vector<IndexType>{nodes[0]->Id(), nodes[1]->Id(), nodes[2]->Id(), nodes[3]->Id()},
        p_properties);
    p_element->Initialize(r_model_part.GetProcessInfo());

    return {&r_model_part, p_element, nodes};
}

void SetNormalOpening(const std::vector<Node::Pointer>& rNodes, const double Opening)
{
    for (std::size_t i = 0; i < rNodes.size(); ++i)
    {
        array_1d<double, 3> displacement = ZeroVector(3);
        if (i >= 2)
        {
            displacement[1] = Opening;
        }
        rNodes[i]->FastGetSolutionStepValue(DISPLACEMENT) = displacement;
    }
}

std::vector<double> CommittedState(const JointElement& rJoint)
{
    std::vector<double> values;
    rJoint.p_element->CalculateOnIntegrationPoints(
        STATE_VARIABLE, values, rJoint.p_model_part->GetProcessInfo());
    return values;
}

void EvaluateTrial(const JointElement& rJoint)
{
    Matrix left_hand_side;
    Vector right_hand_side;
    rJoint.p_element->CalculateLocalSystem(
        left_hand_side, right_hand_side, rJoint.p_model_part->GetProcessInfo());
}

void Commit(const JointElement& rJoint)
{
    rJoint.p_element->FinalizeSolutionStep(rJoint.p_model_part->GetProcessInfo());
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(InterfaceTrialsDoNotCommitConstitutiveHistory, KratosDamFastSuite)
{
    Model model;
    JointElement repeated = CreateJoint(model, "RepeatedTrial", 1, 1);
    SetNormalOpening(repeated.nodes, 5.0e-3);

    // The committed state starts at the damage threshold.
    auto repeated_state = CommittedState(repeated);
    KRATOS_EXPECT_EQ(repeated_state.size(), 2u);
    for (double value : repeated_state)
    {
        KRATOS_EXPECT_NEAR(value, damage_threshold, 1.0e-14);
    }

    // Ten trial evaluations, as performed by ten coupling iterations, must not
    // modify the committed constitutive state.
    for (int i = 0; i < 10; ++i)
    {
        EvaluateTrial(repeated);
    }
    repeated_state = CommittedState(repeated);
    for (double value : repeated_state)
    {
        KRATOS_EXPECT_NEAR(value, damage_threshold, 1.0e-14);
    }

    // A single finalization commits the converged state.
    Commit(repeated);
    repeated_state = CommittedState(repeated);
    const double expected_equivalent_strain = 5.0e-3 / critical_displacement;
    for (double value : repeated_state)
    {
        KRATOS_EXPECT_NEAR(value, expected_equivalent_strain, 1.0e-10);
    }

    // One trial plus one finalization must give exactly the same committed
    // state: the result does not depend on the number of outer iterations.
    JointElement single = CreateJoint(model, "SingleTrial", 2, 2);
    SetNormalOpening(single.nodes, 5.0e-3);
    EvaluateTrial(single);
    Commit(single);
    const auto single_state = CommittedState(single);
    KRATOS_EXPECT_EQ(single_state.size(), repeated_state.size());
    for (std::size_t i = 0; i < single_state.size(); ++i)
    {
        KRATOS_EXPECT_NEAR(single_state[i], repeated_state[i], 1.0e-12);
    }

    // Unloading and reloading below the previous maximum must not reduce or
    // accumulate the committed state (irreversible history variable).
    SetNormalOpening(repeated.nodes, 0.0);
    EvaluateTrial(repeated);
    Commit(repeated);
    for (double value : CommittedState(repeated))
    {
        KRATOS_EXPECT_NEAR(value, expected_equivalent_strain, 1.0e-10);
    }

    SetNormalOpening(repeated.nodes, 2.0e-3);
    EvaluateTrial(repeated);
    Commit(repeated);
    for (double value : CommittedState(repeated))
    {
        KRATOS_EXPECT_NEAR(value, expected_equivalent_strain, 1.0e-10);
    }
}

} // namespace Testing
} // namespace Kratos
