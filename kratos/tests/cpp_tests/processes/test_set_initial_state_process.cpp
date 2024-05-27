//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include"tests/test_utilities/test_constitutive_law.h"
#include"tests/test_utilities/test_element.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"

/* Processes */
#include "processes/set_initial_state_process.h"

namespace Kratos::Testing
{

/**
* Checks the correct work of the set initial state process in 2D
* Test triangle
*/
KRATOS_TEST_CASE_IN_SUITE(ImposingInitialState2D, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    Properties::Pointer p_prop = this_model_part.CreateNewProperties(0);

    TestConstitutiveLaw r_clone_cl = TestConstitutiveLaw();
    p_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    auto p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    auto p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);

    std::vector<Node::Pointer> geom(4);
    geom[0] = p_node_1;
    geom[1] = p_node_2;
    geom[2] = p_node_3;
    geom[3] = p_node_4;
    auto pgeom = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{geom});

    auto p_elem = Kratos::make_intrusive<TestElement>( 0, pgeom, p_prop, TestElement::ResidualType::LINEAR );
    p_elem->Initialize(process_info);
    this_model_part.AddElement(p_elem);

    Vector initial_E = ZeroVector(3);
    initial_E(0) = 0.01;
    initial_E(1) = 0.02;
    initial_E(2) = 0.03;

    Vector initial_S = ZeroVector(3);
    initial_S(0) = 1.0e6;
    initial_S(1) = 2.0e6;
    initial_S(2) = 3.0e6;

    Matrix initial_F = ZeroMatrix(2,2);
    initial_F(0,0) = 0.001;
    initial_F(0,1) = 0.0001;
    initial_F(1,0) = initial_F(0,1);
    initial_F(1,1) = 0.002;

    // Set the initial state
    auto process = SetInitialStateProcess<2>(this_model_part, initial_E, initial_S, initial_F);
    process.ExecuteInitializeSolutionStep();

    const double tolerance = 1.0e-4;
    std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector;
    for (auto i_element = this_model_part.ElementsBegin(); i_element != this_model_part.ElementsEnd(); i_element++) {
        const auto& r_integration_points = i_element->GetGeometry().IntegrationPoints(i_element->GetIntegrationMethod());
        i_element->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_law_vector, this_model_part.GetProcessInfo());

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            auto p_initial_state = constitutive_law_vector[point_number]->pGetInitialState();
            const auto& r_imposed_F = p_initial_state->GetInitialDeformationGradientMatrix();
            const auto& r_imposed_E = p_initial_state->GetInitialStrainVector();
            const auto& r_imposed_S = p_initial_state->GetInitialStressVector();

            for (IndexType component = 0; component < 3; component++) {
                KRATOS_EXPECT_LE(r_imposed_E(component) - initial_E(component), tolerance);
                KRATOS_EXPECT_LE(r_imposed_S(component) - initial_S(component), tolerance);
            }
            KRATOS_EXPECT_LE(r_imposed_F(0,0) - initial_F(0,0), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(0,1) - initial_F(0,1), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(1,1) - initial_F(1,1), tolerance);
        }
    }
}

/**
* Checks the correct work of the set initial state process in 3D
* Test tetrahedra
*/
KRATOS_TEST_CASE_IN_SUITE(ImposingInitialState3D, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    Properties::Pointer p_prop = this_model_part.CreateNewProperties(0);

    TestConstitutiveLaw r_clone_cl = TestConstitutiveLaw();
    p_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    auto p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 1.0);
    auto p_node_2 = this_model_part.CreateNewNode(2, 0.0 , 1.0 , 0.0);
    auto p_node_3 = this_model_part.CreateNewNode(3, 0.0 , 1.0 , 1.0);
    auto p_node_4 = this_model_part.CreateNewNode(4, 1.0 , 1.0 , 0.0);

    std::vector<Node::Pointer> geom(4);
    geom[0] = p_node_1;
    geom[1] = p_node_2;
    geom[2] = p_node_3;
    geom[3] = p_node_4;
    auto pgeom = Kratos::make_shared<Tetrahedra3D4<Node>>(PointerVector<Node>{geom});

    auto p_elem = Kratos::make_intrusive<TestElement>( 0, pgeom, p_prop, TestElement::ResidualType::LINEAR );
    p_elem->Initialize(process_info);
    this_model_part.AddElement(p_elem);

    Vector initial_E = ZeroVector(6);
    initial_E(0) = 0.01;
    initial_E(1) = 0.02;
    initial_E(2) = 0.03;
    initial_E(3) = 0.04;
    initial_E(4) = 0.05;
    initial_E(5) = 0.06;

    Vector initial_S = ZeroVector(6);
    initial_S(0) = 1.0e6;
    initial_S(1) = 2.0e6;
    initial_S(2) = 3.0e6;
    initial_S(3) = 4.0e6;
    initial_S(4) = 5.0e6;
    initial_S(5) = 6.0e6;

    Matrix initial_F = ZeroMatrix(3,3);
    initial_F(0,0) = 0.001;
    initial_F(0,1) = 0.0001;
    initial_F(0,2) = 0.0002;
    initial_F(1,2) = 0.0003;
    initial_F(2,1) = initial_F(1,2);
    initial_F(2,0) = initial_F(0,2);
    initial_F(1,0) = initial_F(0,1);
    initial_F(1,1) = 0.002;
    initial_F(2,2) = 0.003;

    // Set the initial state
    auto process = SetInitialStateProcess<3>(this_model_part, initial_E, initial_S, initial_F);
    process.ExecuteInitializeSolutionStep();

    const double tolerance = 1.0e-4;
    std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector;
    for (auto i_element = this_model_part.ElementsBegin(); i_element != this_model_part.ElementsEnd(); i_element++) {
        const auto& r_integration_points = i_element->GetGeometry().IntegrationPoints(i_element->GetIntegrationMethod());
        i_element->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_law_vector, this_model_part.GetProcessInfo());

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            auto p_initial_state = constitutive_law_vector[point_number]->pGetInitialState();
            const auto& r_imposed_F = p_initial_state->GetInitialDeformationGradientMatrix();
            const auto& r_imposed_E = p_initial_state->GetInitialStrainVector();
            const auto& r_imposed_S = p_initial_state->GetInitialStressVector();

            for (IndexType component = 0; component < 6; component++) {
                KRATOS_EXPECT_LE(r_imposed_E(component) - initial_E(component), tolerance);
                KRATOS_EXPECT_LE(r_imposed_S(component) - initial_S(component), tolerance);
            }
            KRATOS_EXPECT_LE(r_imposed_F(0,0) - initial_F(0,0), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(0,1) - initial_F(0,1), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(1,1) - initial_F(1,1), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(0,1) - initial_F(0,1), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(1,0) - initial_F(1,0), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(1,2) - initial_F(1,2), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(2,1) - initial_F(2,1), tolerance);
            KRATOS_EXPECT_LE(r_imposed_F(2,2) - initial_F(2,2), tolerance);
        }
    }
}
}  // namespace Kratos::Testing.
