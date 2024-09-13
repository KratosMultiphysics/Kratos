// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#include <string>

// Project includes
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/plane_strain_stress_state.h"

#include "geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

class MockElement : public UPwSmallStrainElement<2, 4>
{
public:
    MockElement() : UPwSmallStrainElement<2, 4>() {}

    MockElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : UPwSmallStrainElement(NewId, pGeometry, pProperties, std::make_unique<PlaneStrainStressState>())
    {
    }

    //Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const override
    //{
    //    return Kratos::make_intrusive<MockElement>(NewId, pGeometry, pProperties);
    //}

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        rOutput.resize(4);

        if (rVariable == SHEAR_STIFFNESS) {
            double shear_stiffness = 5.0;
            rOutput[0]             = shear_stiffness;
            rOutput[1]             = shear_stiffness;
            rOutput[2]             = shear_stiffness;
            rOutput[3]             = shear_stiffness;
        } else if (rVariable == CONFINED_STIFFNESS) {
            double confined_stiffness = 4.0;
            rOutput[0]                = confined_stiffness;
            rOutput[1]                = confined_stiffness;
            rOutput[2]                = confined_stiffness;
            rOutput[3]                = confined_stiffness;

        } else if (rVariable == DEGREE_OF_SATURATION) {
            double degree_of_saturation = 0.5;
            rOutput[0]                  = degree_of_saturation;
            rOutput[1]                  = degree_of_saturation;
            rOutput[2]                  = degree_of_saturation;
            rOutput[3]                  = degree_of_saturation;
        } else {
            Element::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
        }
    };
};

ModelPart::ConditionType::Pointer SetUpUPwLysmerAbsorbingCondition2D2NCondition(ModelPart& rModelPart)
{
    // Set the element properties
    auto cond_prop = rModelPart.CreateNewProperties(0);

    // Create the test piping element nodes
    auto node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);

    // Create the test Lysmer condition
    std::vector<ModelPart::IndexType> cond_nodes{1, 2};
    ModelPart::ConditionType::Pointer cond =
        rModelPart.CreateNewCondition("UPwLysmerAbsorbingCondition2D2N", 1, cond_nodes, cond_prop);

    // set properties of the condition
    Vector absorbing_factors = ZeroVector(2);
    absorbing_factors[0]     = 1.0;
    absorbing_factors[1]     = 1.0;

    cond->SetValue(ABSORBING_FACTORS, absorbing_factors);
    cond->SetValue(VIRTUAL_THICKNESS, 100.0);

    // create mock neighbour element
    auto p_neighbour_prop = rModelPart.CreateNewProperties(1);
    p_neighbour_prop->SetValue(POROSITY, 1.0);
    p_neighbour_prop->SetValue(DENSITY_WATER, 1000);
    p_neighbour_prop->SetValue(DENSITY_SOLID, 2000);

    auto neighbour_node_1 = rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    auto neighbour_node_2 = rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);

    // register mock element
    Element::GeometryType::Pointer p_geometry = Kratos::make_shared<Kratos::Quadrilateral2D4<Element::NodeType>>(node_1, node_2, neighbour_node_1, neighbour_node_2);
    const auto mock_element = Kratos::make_intrusive<MockElement>(
        1, p_geometry, Element::PropertiesType::Pointer(0), p_neighbour_prop);
    rModelPart.AddElement(mock_element);

    // add neighbour element to condition
    GlobalPointersVector<Element> vector_of_neighbours;
    vector_of_neighbours.push_back(Kratos::GlobalPointer<Kratos::Element>(mock_element.get()));
    cond->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);

    // Initialize the element
    const auto& r_process_info = rModelPart.GetProcessInfo();
    cond->Initialize(r_process_info);

    return cond;
}

/// <summary>
/// Tests the calculation of the left hand side matrix of the UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateLeftHandSideUPwNormalLysmerAbsorbingCondition, KratosGeoMechanicsFastSuite)
{
    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

    // set up condition
    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition2D2NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Perfrom test, calculate left hand side
    unsigned int conditionSize       = 6;
    Matrix       rLeftHandSideMatrix = ZeroMatrix(conditionSize, conditionSize);

    p_cond->CalculateLeftHandSide(rLeftHandSideMatrix, r_process_info);

    // set expected_results
    Matrix expected_matrix = ZeroMatrix(conditionSize, conditionSize);
    expected_matrix(0, 0)  = 5.0 / 100 / 4;
    expected_matrix(0, 2)  = expected_matrix(0, 0);
    expected_matrix(2, 0)  = expected_matrix(0, 0);
    expected_matrix(2, 2)  = expected_matrix(0, 0);

    expected_matrix(1, 1) = 4.0 / 100 / 4;
    expected_matrix(1, 3) = expected_matrix(1, 1);
    expected_matrix(3, 1) = expected_matrix(1, 1);
    expected_matrix(3, 3) = expected_matrix(1, 1);

    // compare results
    for (unsigned int i = 0; i < conditionSize; ++i) {
        for (unsigned int j = 0; j < conditionSize; ++j) {
            KRATOS_EXPECT_NEAR(rLeftHandSideMatrix(i, j), expected_matrix(i, j), 1.0e-6);
        }
    }
}

/// <summary>
/// Tests the calculation of both the left hand side matrix and the right hand side vector of the UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateLocalSystemUPwNormalLysmerAbsorbingCondition, KratosGeoMechanicsFastSuite)
{
    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition2D2NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    for (auto& node : r_model_part.Nodes()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(WATER_PRESSURE);
    }

    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.0;

    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) = 3.0;
    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) = 4.0;

    // Perform test, calculate left hand side
    constexpr size_t condition_size         = 6;
    Matrix rLeftHandSideMatrix    = ZeroMatrix(condition_size, condition_size);
    Vector right_hand_side_vector = ZeroVector(condition_size);

    p_cond->CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, r_process_info);

    // set expected_results
    Matrix expected_matrix = ZeroMatrix(condition_size, condition_size);
    expected_matrix(0, 0)  = 5.0 / 100 / 4;
    expected_matrix(0, 2)  = expected_matrix(0, 0);
    expected_matrix(2, 0)  = expected_matrix(0, 0);
    expected_matrix(2, 2)  = expected_matrix(0, 0);

    expected_matrix(1, 1) = 4.0 / 100 / 4;
    expected_matrix(1, 3) = expected_matrix(1, 1);
    expected_matrix(3, 1) = expected_matrix(1, 1);
    expected_matrix(3, 3) = expected_matrix(1, 1);

    Vector expected_rhs = ZeroVector(condition_size);
    expected_rhs[0]     = -(expected_matrix(0, 0) * 1.0 + expected_matrix(0, 2) * 3.0);
    expected_rhs[1]     = -(expected_matrix(1, 1) * 2.0 + expected_matrix(1, 3) * 4.0);
    expected_rhs[2]     = -(expected_matrix(2, 0) * 1.0 + expected_matrix(2, 2) * 3.0);
    expected_rhs[3]     = -(expected_matrix(3, 1) * 2.0 + expected_matrix(3, 3) * 4.0);

    // compare results
    KRATOS_EXPECT_MATRIX_NEAR(rLeftHandSideMatrix, expected_matrix, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(right_hand_side_vector, expected_rhs, 1.0e-6);
}
} // namespace Kratos::Testing
