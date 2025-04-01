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

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

class MockElement : public UPwBaseElement
{
public:
    MockElement() : UPwBaseElement() {}

    MockElement(IndexType NewId, const GeometryType::Pointer& rGeometry, const PropertiesType::Pointer& rProperties)
        : UPwBaseElement(NewId, rGeometry, rProperties, std::make_unique<PlaneStrainStressState>())
    {
    }

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        unsigned int       number_of_integration_points = 0;
        const unsigned int n_nodes                      = this->GetGeometry().size();
        const unsigned int n_dim = this->GetGeometry().WorkingSpaceDimension();

        // Get number of integration points based on the element type
        if (n_dim == 2) {
            if (n_nodes == 3 || n_nodes == 6) {
                number_of_integration_points = 3;
            } else if (n_nodes == 4 || n_nodes == 8) {
                number_of_integration_points = 4;
            }
        } else if (n_dim == 3) {
            if (n_nodes == 4 || n_nodes == 10) {
                number_of_integration_points = 4;
            } else if (n_nodes == 8 || n_nodes == 20) {
                number_of_integration_points = 8;
            }
        }
        rOutput.resize(number_of_integration_points);

        // Set some dummy values for the variables
        if (rVariable == SHEAR_STIFFNESS) {
            double shear_stiffness = 5.0;
            rOutput = std::vector<double>(number_of_integration_points, shear_stiffness);
        } else if (rVariable == CONFINED_STIFFNESS) {
            double confined_stiffness = 4.0;
            rOutput = std::vector<double>(number_of_integration_points, confined_stiffness);
        } else if (rVariable == DEGREE_OF_SATURATION) {
            double degree_of_saturation = 0.5;
            rOutput = std::vector<double>(number_of_integration_points, degree_of_saturation);
        } else {
            UPwBaseElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
        }
    };
};

/* Copies upper triangle to bottom. This is used to make a symmetric matrix */
static void SymmetrizeMatrix(Matrix& rMatrix)
{
    const unsigned int matrix_size = rMatrix.size1();

    // Copy the upper triangle to the lower triangle.
    for (unsigned int i = 0; i < matrix_size; ++i) {
        for (unsigned int j = i + 1; j < matrix_size; ++j) {
            rMatrix(j, i) = rMatrix(i, j);
        }
    }
}

/// <summary>
/// Creates a 2D 4N quadrilateral geometry with 2 nodes from the condition and 2 new nodes.
/// </summary>
static Element::GeometryType::Pointer CreateMockGeometry2D4N(Kratos::PointerVector<Node>& rConditionNodes,
                                                             ModelPart& rModelPart)
{
    auto node_1 = rConditionNodes(0);
    auto node_2 = rConditionNodes(1);

    auto node_3 = rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto node_4 = rModelPart.CreateNewNode(4, 1.0, 1.0, 0.0);

    auto p_geometry =
        Kratos::make_shared<Kratos::Quadrilateral2D4<Element::NodeType>>(node_1, node_2, node_3, node_4);
    return p_geometry;
}

/// <summary>
/// Creates a 2D 6N triangle geometry with 3 nodes from the condition and 3 new nodes.
/// </summary>
static Element::GeometryType::Pointer CreateMockGeometry2D6N(Kratos::PointerVector<Node>& rConditionNodes,
                                                             ModelPart& rModelPart)
{
    auto node_1 = rConditionNodes(0);
    auto node_2 = rConditionNodes(1);
    auto node_4 = rConditionNodes(2);

    auto node_3 = rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    auto node_5 = rModelPart.CreateNewNode(5, 1.0, 1.0, 0.0);
    auto node_6 = rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);

    auto p_geometry = Kratos::make_shared<Kratos::Triangle2D6<Element::NodeType>>(
        node_1, node_2, node_3, node_4, node_5, node_6);
    return p_geometry;
}

/// <summary>
/// Creates a 3D 10N tetrahedra geometry with 6 nodes from the condition and 4 new nodes.
/// </summary>
static Element::GeometryType::Pointer CreateMockGeometry3D10N(Kratos::PointerVector<Node>& rConditionNodes,
                                                              ModelPart& rModelPart)
{
    auto node_1 = rConditionNodes(0);
    auto node_2 = rConditionNodes(1);
    auto node_3 = rConditionNodes(2);
    auto node_4 = rConditionNodes(3);
    auto node_5 = rConditionNodes(4);
    auto node_6 = rConditionNodes(5);

    auto node_7  = rModelPart.CreateNewNode(7, 0.0, 0.0, 1.0);
    auto node_8  = rModelPart.CreateNewNode(8, 0.0, 0.0, 0.5);
    auto node_9  = rModelPart.CreateNewNode(9, 0.5, 0.0, 0.5);
    auto node_10 = rModelPart.CreateNewNode(10, 0.0, 0.5, 0.5);

    auto p_geometry = Kratos::make_shared<Kratos::Tetrahedra3D10<Element::NodeType>>(
        node_1, node_2, node_3, node_7, node_4, node_5, node_6, node_8, node_9, node_10);

    return p_geometry;
}

/// <summary>
/// Creates a 3D 20N hexahedra geometry with 8 nodes from the condition and 12 new nodes.
/// </summary>
static Element::GeometryType::Pointer CreateMockGeometry3D20N(Kratos::PointerVector<Node>& rConditionNodes,
                                                              ModelPart& rModelPart)
{
    auto node_1 = rConditionNodes(0);
    auto node_2 = rConditionNodes(1);
    auto node_3 = rConditionNodes(2);
    auto node_4 = rConditionNodes(3);
    auto node_5 = rConditionNodes(4);
    auto node_6 = rConditionNodes(5);
    auto node_7 = rConditionNodes(6);
    auto node_8 = rConditionNodes(7);

    auto node_9  = rModelPart.CreateNewNode(9, 0.0, 0.0, 1.0);
    auto node_10 = rModelPart.CreateNewNode(10, 1.0, 0.0, 1.0);
    auto node_11 = rModelPart.CreateNewNode(11, 1.0, 1.0, 1.0);
    auto node_12 = rModelPart.CreateNewNode(12, 0.0, 1.0, 1.0);
    auto node_13 = rModelPart.CreateNewNode(13, 0.0, 0.0, 0.5);
    auto node_14 = rModelPart.CreateNewNode(14, 1.0, 0.0, 0.5);
    auto node_15 = rModelPart.CreateNewNode(15, 1.0, 1.0, 0.5);
    auto node_16 = rModelPart.CreateNewNode(16, 0.0, 1.0, 0.5);
    auto node_17 = rModelPart.CreateNewNode(17, 0.5, 0.0, 1.0);
    auto node_18 = rModelPart.CreateNewNode(18, 1.0, 0.5, 1.0);
    auto node_19 = rModelPart.CreateNewNode(19, 0.5, 1.0, 1.0);
    auto node_20 = rModelPart.CreateNewNode(20, 0.0, 0.5, 1.0);

    auto p_geometry = Kratos::make_shared<Kratos::Hexahedra3D20<Element::NodeType>>(
        node_1, node_2, node_3, node_4, node_9, node_10, node_11, node_12, node_5, node_6, node_7,
        node_8, node_13, node_14, node_15, node_16, node_17, node_18, node_19, node_20);

    return p_geometry;
}

/// <summary>
/// Sets properties of the condition and the neighbour element. And initializes the condition.
/// </summary>
static void SetPropertiesAndInitialize(ModelPart::ConditionType::Pointer pCondition,
                                       Element::GeometryType::Pointer    pNeighbourGeometry,
                                       ModelPart&                        rModelPart)
{
    // set properties of the condition
    Vector absorbing_factors = ZeroVector(2);
    absorbing_factors[0]     = 1.0;
    absorbing_factors[1]     = 1.0;
    pCondition->SetValue(ABSORBING_FACTORS, absorbing_factors);
    pCondition->SetValue(VIRTUAL_THICKNESS, 100.0);
    // add neighbour element to condition

    auto p_neighbour_prop = rModelPart.CreateNewProperties(1);
    p_neighbour_prop->SetValue(POROSITY, 0.0);
    p_neighbour_prop->SetValue(DENSITY_WATER, 1000);
    p_neighbour_prop->SetValue(DENSITY_SOLID, 2000);

    const auto p_neighbour_element =
        Kratos::make_intrusive<MockElement>(1, pNeighbourGeometry, p_neighbour_prop);
    rModelPart.AddElement(p_neighbour_element);

    GlobalPointersVector<Element> vector_of_neighbours;
    vector_of_neighbours.push_back(Kratos::GlobalPointer<Kratos::Element>(p_neighbour_element.get()));
    pCondition->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);

    p_neighbour_element->SetValue(POROSITY, 0.0);
    p_neighbour_element->SetValue(DENSITY_WATER, 1000);
    p_neighbour_element->SetValue(DENSITY_SOLID, 2000);

    // Initialize the element
    const auto& r_process_info = rModelPart.GetProcessInfo();
    pCondition->Initialize(r_process_info);
}

/// <summary>
/// Sets up a 2D 2N UPwLysmerAbsorbingCondition with a 2D 4N neighbour element
/// </summary>
static ModelPart::ConditionType::Pointer SetUpUPwLysmerAbsorbingCondition2D2NCondition(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Set the element properties
    auto cond_prop = rModelPart.CreateNewProperties(0);

    auto node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);

    // Create the test Lysmer condition
    std::vector<ModelPart::IndexType> cond_nodes{1, 2};
    ModelPart::ConditionType::Pointer p_cond =
        rModelPart.CreateNewCondition("UPwLysmerAbsorbingCondition2D2N", 1, cond_nodes, cond_prop);

    Element::GeometryType::Pointer p_neighbour_geometry =
        CreateMockGeometry2D4N(p_cond->pGetGeometry()->Points(), rModelPart);

    for (auto& node : rModelPart.Nodes()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(WATER_PRESSURE);
    }

    SetPropertiesAndInitialize(p_cond, p_neighbour_geometry, rModelPart);

    return p_cond;
}

/// <summary>
/// Sets up a 2D 3N UPwLysmerAbsorbingCondition with a 2D 6N neighbour element
/// </summary>
static ModelPart::ConditionType::Pointer SetUpUPwLysmerAbsorbingCondition2D3NCondition(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Set the element properties
    auto cond_prop = rModelPart.CreateNewProperties(0);

    auto node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto node_3 = rModelPart.CreateNewNode(3, 0.5, 0.0, 0.0);

    // Create the test Lysmer condition
    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3};
    ModelPart::ConditionType::Pointer p_cond =
        rModelPart.CreateNewCondition("UPwLysmerAbsorbingCondition2D3N", 1, cond_nodes, cond_prop);

    Element::GeometryType::Pointer p_neighbour_geometry =
        CreateMockGeometry2D6N(p_cond->pGetGeometry()->Points(), rModelPart);

    for (auto& node : rModelPart.Nodes()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(WATER_PRESSURE);
    }

    SetPropertiesAndInitialize(p_cond, p_neighbour_geometry, rModelPart);

    return p_cond;
}

/// <summary>
/// Sets up a 3D 6N UPwLysmerAbsorbingCondition with a 3D 10N neighbour element
/// </summary>
static ModelPart::ConditionType::Pointer SetUpUPwLysmerAbsorbingCondition3D6NCondition(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Set the element properties
    auto cond_prop = rModelPart.CreateNewProperties(0);

    auto node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto node_3 = rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto node_4 = rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    auto node_5 = rModelPart.CreateNewNode(5, 0.5, 0.5, 0.0);
    auto node_6 = rModelPart.CreateNewNode(6, 0.0, 0.5, 0.0);

    // Create the test Lysmer condition
    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6};
    ModelPart::ConditionType::Pointer p_cond =
        rModelPart.CreateNewCondition("UPwLysmerAbsorbingCondition3D6N", 1, cond_nodes, cond_prop);
    Element::GeometryType::Pointer p_neighbour_geometry =
        CreateMockGeometry3D10N(p_cond->pGetGeometry()->Points(), rModelPart);

    for (auto& node : rModelPart.Nodes()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
        node.AddDof(WATER_PRESSURE);
    }

    SetPropertiesAndInitialize(p_cond, p_neighbour_geometry, rModelPart);
    return p_cond;
}

/// <summary>
/// Sets up a 3D 8N UPwLysmerAbsorbingCondition with a 3D 20N neighbour element
/// </summary>
static ModelPart::ConditionType::Pointer SetUpUPwLysmerAbsorbingCondition3D8NCondition(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Set the element properties
    auto cond_prop = rModelPart.CreateNewProperties(0);

    auto node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto node_3 = rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    auto node_4 = rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    auto node_5 = rModelPart.CreateNewNode(5, 0.5, 0.0, 0.0);
    auto node_6 = rModelPart.CreateNewNode(6, 1.0, 0.5, 0.0);
    auto node_7 = rModelPart.CreateNewNode(7, 0.5, 1.0, 0.0);
    auto node_8 = rModelPart.CreateNewNode(8, 0.0, 0.5, 0.0);

    // Create the test Lysmer condition
    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6, 7, 8};
    ModelPart::ConditionType::Pointer p_cond =
        rModelPart.CreateNewCondition("UPwLysmerAbsorbingCondition3D8N", 1, cond_nodes, cond_prop);
    Element::GeometryType::Pointer p_neighbour_geometry =
        CreateMockGeometry3D20N(p_cond->pGetGeometry()->Points(), rModelPart);

    for (auto& node : rModelPart.Nodes()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
        node.AddDof(WATER_PRESSURE);
    }

    SetPropertiesAndInitialize(p_cond, p_neighbour_geometry, rModelPart);

    return p_cond;
}

/// <summary>
/// Tests the calculation of the LHS matrix of a 2D 2N UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateLeftHandSideUPwNormalLysmerAbsorbingCondition, KratosGeoMechanicsFastSuite)
{
    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition2D2NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Perform test, calculate left hand side
    constexpr size_t conditionSize       = 6;
    Matrix           rLeftHandSideMatrix = ZeroMatrix(conditionSize, conditionSize);

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
    KRATOS_EXPECT_MATRIX_NEAR(rLeftHandSideMatrix, expected_matrix, 1.0e-6);
}

/// <summary>
/// Tests the calculation of the LHS matrix and RHS vector of a 2D 2N UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateLocalSystemUPwNormalLysmerAbsorbingCondition2D2N, KratosGeoMechanicsFastSuite)
{
    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition2D2NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.0;

    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) = 3.0;
    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) = 4.0;

    // Perform test, calculate left hand side
    constexpr size_t condition_size         = 6;
    Matrix           rLeftHandSideMatrix    = ZeroMatrix(condition_size, condition_size);
    Vector           right_hand_side_vector = ZeroVector(condition_size);

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

/// <summary>
/// Tests the calculation of the LHS matrix and RHS vector of a 2D 3N UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateLocalSystemUPwNormalLysmerAbsorbingCondition2D3N, KratosGeoMechanicsFastSuite)
{
    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition2D3NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.0;

    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) = 3.0;
    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) = 4.0;

    p_cond->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT_X) = 2.0;
    p_cond->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT_Y) = 3.0;

    // Perform test, calculate left hand side
    constexpr size_t condition_size         = 9;
    Matrix           rLeftHandSideMatrix    = ZeroMatrix(condition_size, condition_size);
    Vector           right_hand_side_vector = ZeroVector(condition_size);

    p_cond->CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, r_process_info);

    // set expected_results

    const double shear_stiffness_over_virt_thickness    = 5.0 / 100;
    const double confined_stiffness_over_virt_thickness = 4.0 / 100;

    Matrix expected_matrix = ZeroMatrix(condition_size, condition_size);
    expected_matrix(0, 0)  = shear_stiffness_over_virt_thickness / 9;
    expected_matrix(0, 2)  = -shear_stiffness_over_virt_thickness / 18;
    expected_matrix(0, 4)  = expected_matrix(0, 0);
    expected_matrix(2, 0)  = expected_matrix(0, 2);
    expected_matrix(4, 0)  = expected_matrix(0, 0);

    expected_matrix(2, 2) = expected_matrix(0, 0);
    expected_matrix(2, 4) = expected_matrix(0, 0);

    expected_matrix(4, 2) = expected_matrix(0, 0);
    expected_matrix(4, 4) = shear_stiffness_over_virt_thickness * 4 / 9;

    expected_matrix(1, 1) = confined_stiffness_over_virt_thickness / 9;
    expected_matrix(1, 3) = -confined_stiffness_over_virt_thickness / 18;
    expected_matrix(1, 5) = expected_matrix(1, 1);
    expected_matrix(3, 1) = expected_matrix(1, 3);
    expected_matrix(3, 3) = expected_matrix(1, 1);

    expected_matrix(3, 5) = expected_matrix(1, 1);
    expected_matrix(5, 1) = expected_matrix(1, 1);

    expected_matrix(5, 3) = expected_matrix(1, 1);
    expected_matrix(5, 5) = confined_stiffness_over_virt_thickness * 4 / 9;

    Vector expected_rhs = ZeroVector(condition_size);
    expected_rhs[0] =
        -(expected_matrix(0, 0) * 1.0 + expected_matrix(0, 2) * 3.0 + expected_matrix(0, 4) * 2.0);
    expected_rhs[1] =
        -(expected_matrix(1, 1) * 2.0 + expected_matrix(1, 3) * 4.0 + expected_matrix(1, 5) * 3.0);
    expected_rhs[2] =
        -(expected_matrix(2, 0) * 1.0 + expected_matrix(2, 2) * 3.0 + expected_matrix(2, 4) * 2.0);
    expected_rhs[3] =
        -(expected_matrix(3, 1) * 2.0 + expected_matrix(3, 3) * 4.0 + expected_matrix(3, 5) * 3.0);
    expected_rhs[4] =
        -(expected_matrix(4, 0) * 1.0 + expected_matrix(4, 2) * 3.0 + expected_matrix(4, 4) * 2.0);
    expected_rhs[5] =
        -(expected_matrix(5, 1) * 2.0 + expected_matrix(5, 3) * 4.0 + expected_matrix(5, 5) * 3.0);

    // compare results
    KRATOS_EXPECT_MATRIX_NEAR(rLeftHandSideMatrix, expected_matrix, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(right_hand_side_vector, expected_rhs, 1.0e-6);
}

/// <summary>
/// Tests the calculation of the damping matrix of a 2D 3N UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateDampingMatrixUPwNormalLysmerAbsorbingCondition2D3N, KratosGeoMechanicsFastSuite)
{
    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition2D3NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Perform test, calculate left hand side
    constexpr size_t condition_size         = 9;
    Matrix           rLeftHandSideMatrix    = ZeroMatrix(condition_size, condition_size);
    Vector           right_hand_side_vector = ZeroVector(condition_size);

    p_cond->CalculateDampingMatrix(rLeftHandSideMatrix, r_process_info);

    // set expected_results
    const double vp = sqrt(4.0 / 2000); // sqrt(confined_stiffness / density)
    const double vs = sqrt(5.0 / 2000); // sqrt(shear_stiffness / density)

    const double perpendicular_damping = vp * 2000 * 1; // vp * density * perpendicular_damping_factor
    const double shear_damping = vs * 2000 * 1;         // vs * density * shear_damping_factor

    Matrix expected_matrix = ZeroMatrix(condition_size, condition_size);
    expected_matrix(0, 0)  = shear_damping / 9;
    expected_matrix(0, 2)  = -shear_damping / 18;
    expected_matrix(0, 4)  = expected_matrix(0, 0);
    expected_matrix(2, 0)  = expected_matrix(0, 2);
    expected_matrix(4, 0)  = expected_matrix(0, 0);

    expected_matrix(2, 2) = expected_matrix(0, 0);
    expected_matrix(2, 4) = expected_matrix(0, 0);

    expected_matrix(4, 2) = expected_matrix(0, 0);
    expected_matrix(4, 4) = shear_damping * 4 / 9;

    expected_matrix(1, 1) = perpendicular_damping / 9;
    expected_matrix(1, 3) = -perpendicular_damping / 18;
    expected_matrix(1, 5) = expected_matrix(1, 1);
    expected_matrix(3, 1) = expected_matrix(1, 3);
    expected_matrix(3, 3) = expected_matrix(1, 1);

    expected_matrix(3, 5) = expected_matrix(1, 1);
    expected_matrix(5, 1) = expected_matrix(1, 1);

    expected_matrix(5, 3) = expected_matrix(1, 1);
    expected_matrix(5, 5) = perpendicular_damping * 4 / 9;

    // compare results
    KRATOS_EXPECT_MATRIX_NEAR(rLeftHandSideMatrix, expected_matrix, 1.0e-6);
}

/// <summary>
/// Tests the calculation of the LHS matrix and RHS vector of a 3D 6N UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateLocalSystemUPwNormalLysmerAbsorbingCondition3D6N, KratosGeoMechanicsFastSuite)
{
    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition3D6NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.0;
    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z) = 3.0;

    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) = 2.0;
    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) = 3.0;
    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) = 4.0;

    p_cond->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT_X) = 3.0;
    p_cond->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT_Y) = 4.0;
    p_cond->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT_Z) = 5.0;

    p_cond->GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT_X) = 1.5;
    p_cond->GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.5;
    p_cond->GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT_Z) = 3.5;

    p_cond->GetGeometry()[4].FastGetSolutionStepValue(DISPLACEMENT_X) = 2.5;
    p_cond->GetGeometry()[4].FastGetSolutionStepValue(DISPLACEMENT_Y) = 3.5;
    p_cond->GetGeometry()[4].FastGetSolutionStepValue(DISPLACEMENT_Z) = 4.5;

    p_cond->GetGeometry()[5].FastGetSolutionStepValue(DISPLACEMENT_X) = 2.0;
    p_cond->GetGeometry()[5].FastGetSolutionStepValue(DISPLACEMENT_Y) = 3.0;
    p_cond->GetGeometry()[5].FastGetSolutionStepValue(DISPLACEMENT_Z) = 4.0;

    // Perform test, calculate left hand side
    constexpr size_t condition_size         = 24;
    Matrix           rLeftHandSideMatrix    = ZeroMatrix(condition_size, condition_size);
    Vector           right_hand_side_vector = ZeroVector(condition_size);

    p_cond->CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, r_process_info);

    // set expected_results

    const double shear_stiffness_over_virt_thickness    = 5.0 / 100;
    const double confined_stiffness_over_virt_thickness = 4.0 / 100;

    Matrix expected_matrix = ZeroMatrix(condition_size, condition_size);
    expected_matrix(0, 0)  = shear_stiffness_over_virt_thickness / 81;
    expected_matrix(0, 3)  = -shear_stiffness_over_virt_thickness / 162;
    expected_matrix(0, 6)  = expected_matrix(0, 3);
    expected_matrix(0, 9)  = -expected_matrix(0, 3);
    expected_matrix(0, 12) = -expected_matrix(0, 0);
    expected_matrix(0, 15) = -expected_matrix(0, 3);

    expected_matrix(1, 1)  = shear_stiffness_over_virt_thickness / 81;
    expected_matrix(1, 4)  = -shear_stiffness_over_virt_thickness / 162;
    expected_matrix(1, 7)  = expected_matrix(1, 4);
    expected_matrix(1, 10) = -expected_matrix(1, 4);
    expected_matrix(1, 13) = -expected_matrix(1, 1);
    expected_matrix(1, 16) = -expected_matrix(1, 4);

    expected_matrix(2, 2)  = confined_stiffness_over_virt_thickness / 81;
    expected_matrix(2, 5)  = -confined_stiffness_over_virt_thickness / 162;
    expected_matrix(2, 8)  = expected_matrix(2, 5);
    expected_matrix(2, 11) = -expected_matrix(2, 5);
    expected_matrix(2, 14) = -expected_matrix(2, 2);
    expected_matrix(2, 17) = -expected_matrix(2, 5);

    expected_matrix(3, 3)  = expected_matrix(0, 0);
    expected_matrix(3, 6)  = expected_matrix(0, 3);
    expected_matrix(3, 9)  = -expected_matrix(0, 3);
    expected_matrix(3, 12) = -expected_matrix(0, 3);
    expected_matrix(3, 15) = -expected_matrix(0, 0);

    expected_matrix(4, 4)  = expected_matrix(1, 1);
    expected_matrix(4, 7)  = expected_matrix(1, 4);
    expected_matrix(4, 10) = -expected_matrix(1, 4);
    expected_matrix(4, 13) = -expected_matrix(1, 4);
    expected_matrix(4, 16) = -expected_matrix(1, 1);

    expected_matrix(5, 5)  = expected_matrix(2, 2);
    expected_matrix(5, 8)  = expected_matrix(2, 5);
    expected_matrix(5, 11) = -expected_matrix(2, 5);
    expected_matrix(5, 14) = -expected_matrix(2, 5);
    expected_matrix(5, 17) = -expected_matrix(2, 2);

    expected_matrix(6, 6)  = expected_matrix(0, 0);
    expected_matrix(6, 9)  = -expected_matrix(0, 0);
    expected_matrix(6, 12) = -expected_matrix(0, 3);
    expected_matrix(6, 15) = -expected_matrix(0, 3);

    expected_matrix(7, 7)  = expected_matrix(1, 1);
    expected_matrix(7, 10) = -expected_matrix(1, 1);
    expected_matrix(7, 13) = -expected_matrix(1, 4);
    expected_matrix(7, 16) = -expected_matrix(1, 4);

    expected_matrix(8, 8)  = expected_matrix(2, 2);
    expected_matrix(8, 11) = -expected_matrix(2, 2);
    expected_matrix(8, 14) = -expected_matrix(2, 5);
    expected_matrix(8, 17) = -expected_matrix(2, 5);

    expected_matrix(9, 9)  = shear_stiffness_over_virt_thickness * 11 / 162;
    expected_matrix(9, 12) = shear_stiffness_over_virt_thickness * 4 / 81;
    expected_matrix(9, 15) = expected_matrix(9, 12);

    expected_matrix(10, 10) = shear_stiffness_over_virt_thickness * 11 / 162;
    expected_matrix(10, 13) = shear_stiffness_over_virt_thickness * 4 / 81;
    expected_matrix(10, 16) = expected_matrix(10, 13);

    expected_matrix(11, 11) = confined_stiffness_over_virt_thickness * 11 / 162;
    expected_matrix(11, 14) = confined_stiffness_over_virt_thickness * 4 / 81;
    expected_matrix(11, 17) = expected_matrix(11, 14);

    expected_matrix(12, 12) = expected_matrix(9, 9);
    expected_matrix(12, 15) = expected_matrix(9, 12);

    expected_matrix(13, 13) = expected_matrix(10, 10);
    expected_matrix(13, 16) = expected_matrix(10, 13);

    expected_matrix(14, 14) = expected_matrix(11, 11);
    expected_matrix(14, 17) = expected_matrix(11, 14);

    expected_matrix(15, 15) = expected_matrix(9, 9);

    expected_matrix(16, 16) = expected_matrix(10, 10);

    expected_matrix(17, 17) = expected_matrix(11, 11);

    SymmetrizeMatrix(expected_matrix);

    Vector expected_rhs = ZeroVector(condition_size);
    expected_rhs[0] =
        -(expected_matrix(0, 0) * 1.0 + expected_matrix(0, 3) * 2.0 + expected_matrix(0, 6) * 3.0 +
          expected_matrix(0, 9) * 1.5 + expected_matrix(0, 12) * 2.5 + expected_matrix(0, 15) * 2.0);
    expected_rhs[1] =
        -(expected_matrix(1, 1) * 2.0 + expected_matrix(1, 4) * 3.0 + expected_matrix(1, 7) * 4.0 +
          expected_matrix(1, 10) * 2.5 + expected_matrix(1, 13) * 3.5 + expected_matrix(1, 16) * 3.0);
    expected_rhs[2] =
        -(expected_matrix(2, 2) * 3.0 + expected_matrix(2, 5) * 4.0 + expected_matrix(2, 8) * 5.0 +
          expected_matrix(2, 11) * 3.5 + expected_matrix(2, 14) * 4.5 + expected_matrix(2, 17) * 4.0);

    expected_rhs[3] =
        -(expected_matrix(3, 0) * 1.0 + expected_matrix(3, 3) * 2.0 + expected_matrix(3, 6) * 3.0 +
          expected_matrix(3, 9) * 1.5 + expected_matrix(3, 12) * 2.5 + expected_matrix(3, 15) * 2.0);
    expected_rhs[4] =
        -(expected_matrix(4, 1) * 2.0 + expected_matrix(4, 4) * 3.0 + expected_matrix(4, 7) * 4.0 +
          expected_matrix(4, 10) * 2.5 + expected_matrix(4, 13) * 3.5 + expected_matrix(4, 16) * 3.0);
    expected_rhs[5] =
        -(expected_matrix(5, 2) * 3.0 + expected_matrix(5, 5) * 4.0 + expected_matrix(5, 8) * 5.0 +
          expected_matrix(5, 11) * 3.5 + expected_matrix(5, 14) * 4.5 + expected_matrix(5, 17) * 4.0);

    expected_rhs[6] =
        -(expected_matrix(6, 0) * 1.0 + expected_matrix(6, 3) * 2.0 + expected_matrix(6, 6) * 3.0 +
          expected_matrix(6, 9) * 1.5 + expected_matrix(6, 12) * 2.5 + expected_matrix(6, 15) * 2.0);
    expected_rhs[7] =
        -(expected_matrix(7, 1) * 2.0 + expected_matrix(7, 4) * 3.0 + expected_matrix(7, 7) * 4.0 +
          expected_matrix(7, 10) * 2.5 + expected_matrix(7, 13) * 3.5 + expected_matrix(7, 16) * 3.0);
    expected_rhs[8] =
        -(expected_matrix(8, 2) * 3.0 + expected_matrix(8, 5) * 4.0 + expected_matrix(8, 8) * 5.0 +
          expected_matrix(8, 11) * 3.5 + expected_matrix(8, 14) * 4.5 + expected_matrix(8, 17) * 4.0);

    expected_rhs[9] =
        -(expected_matrix(9, 0) * 1.0 + expected_matrix(9, 3) * 2.0 + expected_matrix(9, 6) * 3.0 +
          expected_matrix(9, 9) * 1.5 + expected_matrix(9, 12) * 2.5 + expected_matrix(9, 15) * 2.0);
    expected_rhs[10] =
        -(expected_matrix(10, 1) * 2.0 + expected_matrix(10, 4) * 3.0 + expected_matrix(10, 7) * 4.0 +
          expected_matrix(10, 10) * 2.5 + expected_matrix(10, 13) * 3.5 + expected_matrix(10, 16) * 3.0);
    expected_rhs[11] =
        -(expected_matrix(11, 2) * 3.0 + expected_matrix(11, 5) * 4.0 + expected_matrix(11, 8) * 5.0 +
          expected_matrix(11, 11) * 3.5 + expected_matrix(11, 14) * 4.5 + expected_matrix(11, 17) * 4.0);

    expected_rhs[12] =
        -(expected_matrix(12, 0) * 1.0 + expected_matrix(12, 3) * 2.0 + expected_matrix(12, 6) * 3.0 +
          expected_matrix(12, 9) * 1.5 + expected_matrix(12, 12) * 2.5 + expected_matrix(12, 15) * 2.0);
    expected_rhs[13] =
        -(expected_matrix(13, 1) * 2.0 + expected_matrix(13, 4) * 3.0 + expected_matrix(13, 7) * 4.0 +
          expected_matrix(13, 10) * 2.5 + expected_matrix(13, 13) * 3.5 + expected_matrix(13, 16) * 3.0);
    expected_rhs[14] =
        -(expected_matrix(14, 2) * 3.0 + expected_matrix(14, 5) * 4.0 + expected_matrix(14, 8) * 5.0 +
          expected_matrix(14, 11) * 3.5 + expected_matrix(14, 14) * 4.5 + expected_matrix(14, 17) * 4.0);

    expected_rhs[15] =
        -(expected_matrix(15, 0) * 1.0 + expected_matrix(15, 3) * 2.0 + expected_matrix(15, 6) * 3.0 +
          expected_matrix(15, 9) * 1.5 + expected_matrix(15, 12) * 2.5 + expected_matrix(15, 15) * 2.0);
    expected_rhs[16] =
        -(expected_matrix(16, 1) * 2.0 + expected_matrix(16, 4) * 3.0 + expected_matrix(16, 7) * 4.0 +
          expected_matrix(16, 10) * 2.5 + expected_matrix(16, 13) * 3.5 + expected_matrix(16, 16) * 3.0);
    expected_rhs[17] =
        -(expected_matrix(17, 2) * 3.0 + expected_matrix(17, 5) * 4.0 + expected_matrix(17, 8) * 5.0 +
          expected_matrix(17, 11) * 3.5 + expected_matrix(17, 14) * 4.5 + expected_matrix(17, 17) * 4.0);

    // compare results
    KRATOS_EXPECT_MATRIX_NEAR(rLeftHandSideMatrix, expected_matrix, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(right_hand_side_vector, expected_rhs, 1.0e-6);

    std::cout << "calculated rhs (2): " << right_hand_side_vector(2) << std::endl;
}

/// <summary>
/// Tests the calculation of the damping matrix of a 3D 6N UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateDampingMatrixUPwNormalLysmerAbsorbingCondition3D6N, KratosGeoMechanicsFastSuite)
{
    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition3D6NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Perform test, calculate left hand side
    constexpr size_t condition_size         = 24;
    Matrix           damping_matrix         = ZeroMatrix(condition_size, condition_size);
    Vector           right_hand_side_vector = ZeroVector(condition_size);

    p_cond->CalculateDampingMatrix(damping_matrix, r_process_info);

    // set expected_results
    const double vp = sqrt(4.0 / 2000); // sqrt(confined_stiffness / density)
    const double vs = sqrt(5.0 / 2000); // sqrt(shear_stiffness / density)

    const double perpendicular_damping = vp * 2000 * 1; // vp * density * perpendicular_damping_factor
    const double shear_damping = vs * 2000 * 1;         // vs * density * shear_damping_factor

    Matrix expected_matrix = ZeroMatrix(condition_size, condition_size);
    expected_matrix(0, 0)  = shear_damping / 81;
    expected_matrix(0, 3)  = -shear_damping / 162;
    expected_matrix(0, 6)  = expected_matrix(0, 3);
    expected_matrix(0, 9)  = -expected_matrix(0, 3);
    expected_matrix(0, 12) = -expected_matrix(0, 0);
    expected_matrix(0, 15) = -expected_matrix(0, 3);

    expected_matrix(1, 1)  = shear_damping / 81;
    expected_matrix(1, 4)  = -shear_damping / 162;
    expected_matrix(1, 7)  = expected_matrix(1, 4);
    expected_matrix(1, 10) = -expected_matrix(1, 4);
    expected_matrix(1, 13) = -expected_matrix(1, 1);
    expected_matrix(1, 16) = -expected_matrix(1, 4);

    expected_matrix(2, 2)  = perpendicular_damping / 81;
    expected_matrix(2, 5)  = -perpendicular_damping / 162;
    expected_matrix(2, 8)  = expected_matrix(2, 5);
    expected_matrix(2, 11) = -expected_matrix(2, 5);
    expected_matrix(2, 14) = -expected_matrix(2, 2);
    expected_matrix(2, 17) = -expected_matrix(2, 5);

    expected_matrix(3, 3)  = expected_matrix(0, 0);
    expected_matrix(3, 6)  = expected_matrix(0, 3);
    expected_matrix(3, 9)  = -expected_matrix(0, 3);
    expected_matrix(3, 12) = -expected_matrix(0, 3);
    expected_matrix(3, 15) = -expected_matrix(0, 0);

    expected_matrix(4, 4)  = expected_matrix(1, 1);
    expected_matrix(4, 7)  = expected_matrix(1, 4);
    expected_matrix(4, 10) = -expected_matrix(1, 4);
    expected_matrix(4, 13) = -expected_matrix(1, 4);
    expected_matrix(4, 16) = -expected_matrix(1, 1);

    expected_matrix(5, 5)  = expected_matrix(2, 2);
    expected_matrix(5, 8)  = expected_matrix(2, 5);
    expected_matrix(5, 11) = -expected_matrix(2, 5);
    expected_matrix(5, 14) = -expected_matrix(2, 5);
    expected_matrix(5, 17) = -expected_matrix(2, 2);

    expected_matrix(6, 6)  = expected_matrix(0, 0);
    expected_matrix(6, 9)  = -expected_matrix(0, 0);
    expected_matrix(6, 12) = -expected_matrix(0, 3);
    expected_matrix(6, 15) = -expected_matrix(0, 3);

    expected_matrix(7, 7)  = expected_matrix(1, 1);
    expected_matrix(7, 10) = -expected_matrix(1, 1);
    expected_matrix(7, 13) = -expected_matrix(1, 4);
    expected_matrix(7, 16) = -expected_matrix(1, 4);

    expected_matrix(8, 8)  = expected_matrix(2, 2);
    expected_matrix(8, 11) = -expected_matrix(2, 2);
    expected_matrix(8, 14) = -expected_matrix(2, 5);
    expected_matrix(8, 17) = -expected_matrix(2, 5);

    expected_matrix(9, 9)  = shear_damping * 11 / 162;
    expected_matrix(9, 12) = shear_damping * 4 / 81;
    expected_matrix(9, 15) = expected_matrix(9, 12);

    expected_matrix(10, 10) = shear_damping * 11 / 162;
    expected_matrix(10, 13) = shear_damping * 4 / 81;
    expected_matrix(10, 16) = expected_matrix(10, 13);

    expected_matrix(11, 11) = perpendicular_damping * 11 / 162;
    expected_matrix(11, 14) = perpendicular_damping * 4 / 81;
    expected_matrix(11, 17) = expected_matrix(11, 14);

    expected_matrix(12, 12) = expected_matrix(9, 9);
    expected_matrix(12, 15) = expected_matrix(9, 12);

    expected_matrix(13, 13) = expected_matrix(10, 10);
    expected_matrix(13, 16) = expected_matrix(10, 13);

    expected_matrix(14, 14) = expected_matrix(11, 11);
    expected_matrix(14, 17) = expected_matrix(11, 14);

    expected_matrix(15, 15) = expected_matrix(9, 9);

    expected_matrix(16, 16) = expected_matrix(10, 10);

    expected_matrix(17, 17) = expected_matrix(11, 11);

    SymmetrizeMatrix(expected_matrix);

    // compare results
    KRATOS_EXPECT_MATRIX_NEAR(damping_matrix, expected_matrix, 1.0e-6);
}

/// <summary>
/// Tests the calculation of the LHS matrix and RHS vector of a 3D 8N UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateLocalSystemUPwNormalLysmerAbsorbingCondition3D8N, KratosGeoMechanicsFastSuite)
{
    GTEST_SKIP() << "This test is skipped as it depends on #13269" << std::endl;

    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition3D8NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.0;
    p_cond->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z) = 3.0;

    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) = 2.0;
    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) = 3.0;
    p_cond->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) = 4.0;

    p_cond->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT_X) = 3.0;
    p_cond->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT_Y) = 4.0;
    p_cond->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT_Z) = 5.0;

    p_cond->GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT_X) = 4.0;
    p_cond->GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT_Y) = 5.0;
    p_cond->GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT_Z) = 6.0;

    p_cond->GetGeometry()[4].FastGetSolutionStepValue(DISPLACEMENT_X) = 1.5;
    p_cond->GetGeometry()[4].FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.5;
    p_cond->GetGeometry()[4].FastGetSolutionStepValue(DISPLACEMENT_Z) = 3.5;

    p_cond->GetGeometry()[5].FastGetSolutionStepValue(DISPLACEMENT_X) = 2.5;
    p_cond->GetGeometry()[5].FastGetSolutionStepValue(DISPLACEMENT_Y) = 3.5;
    p_cond->GetGeometry()[5].FastGetSolutionStepValue(DISPLACEMENT_Z) = 4.5;

    p_cond->GetGeometry()[6].FastGetSolutionStepValue(DISPLACEMENT_X) = 3.5;
    p_cond->GetGeometry()[6].FastGetSolutionStepValue(DISPLACEMENT_Y) = 4.5;
    p_cond->GetGeometry()[6].FastGetSolutionStepValue(DISPLACEMENT_Z) = 5.5;

    p_cond->GetGeometry()[7].FastGetSolutionStepValue(DISPLACEMENT_X) = 2.5;
    p_cond->GetGeometry()[7].FastGetSolutionStepValue(DISPLACEMENT_Y) = 3.5;
    p_cond->GetGeometry()[7].FastGetSolutionStepValue(DISPLACEMENT_Z) = 4.5;

    // Perform test, calculate left hand side
    constexpr size_t condition_size         = 32;
    Matrix           rLeftHandSideMatrix    = ZeroMatrix(condition_size, condition_size);
    Vector           right_hand_side_vector = ZeroVector(condition_size);

    p_cond->CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, r_process_info);

    // set expected_results
    const double shear_stiffness_over_virt_thickness    = 5.0 / 100;
    const double confined_stiffness_over_virt_thickness = 4.0 / 100;

    Matrix expected_matrix = ZeroMatrix(condition_size, condition_size);
    expected_matrix(0, 0)  = shear_stiffness_over_virt_thickness / 30;
    expected_matrix(0, 3)  = shear_stiffness_over_virt_thickness / 90;
    expected_matrix(0, 6)  = shear_stiffness_over_virt_thickness / 60;
    expected_matrix(0, 9)  = expected_matrix(0, 3);
    expected_matrix(0, 12) = -shear_stiffness_over_virt_thickness / 30;
    expected_matrix(0, 15) = -shear_stiffness_over_virt_thickness * 2 / 45;
    expected_matrix(0, 18) = expected_matrix(0, 15);
    expected_matrix(0, 21) = expected_matrix(0, 12);

    expected_matrix(1, 1)  = shear_stiffness_over_virt_thickness / 30;
    expected_matrix(1, 4)  = shear_stiffness_over_virt_thickness / 90;
    expected_matrix(1, 7)  = shear_stiffness_over_virt_thickness / 60;
    expected_matrix(1, 10) = expected_matrix(1, 4);
    expected_matrix(1, 13) = -shear_stiffness_over_virt_thickness / 30;
    expected_matrix(1, 16) = -shear_stiffness_over_virt_thickness * 2 / 45;
    expected_matrix(1, 19) = expected_matrix(1, 16);
    expected_matrix(1, 22) = expected_matrix(1, 13);

    expected_matrix(2, 2)  = confined_stiffness_over_virt_thickness / 30;
    expected_matrix(2, 5)  = confined_stiffness_over_virt_thickness / 90;
    expected_matrix(2, 8)  = confined_stiffness_over_virt_thickness / 60;
    expected_matrix(2, 11) = expected_matrix(2, 5);
    expected_matrix(2, 14) = -confined_stiffness_over_virt_thickness / 30;
    expected_matrix(2, 17) = -confined_stiffness_over_virt_thickness * 2 / 45;
    expected_matrix(2, 20) = expected_matrix(2, 17);
    expected_matrix(2, 23) = expected_matrix(2, 14);

    expected_matrix(3, 3)  = expected_matrix(0, 0);
    expected_matrix(3, 6)  = expected_matrix(0, 3);
    expected_matrix(3, 9)  = expected_matrix(0, 6);
    expected_matrix(3, 12) = expected_matrix(0, 12);
    expected_matrix(3, 15) = expected_matrix(0, 12);
    expected_matrix(3, 18) = expected_matrix(0, 15);
    expected_matrix(3, 21) = expected_matrix(0, 18);

    expected_matrix(4, 4)  = expected_matrix(1, 1);
    expected_matrix(4, 7)  = expected_matrix(1, 4);
    expected_matrix(4, 10) = expected_matrix(1, 7);
    expected_matrix(4, 13) = expected_matrix(1, 13);
    expected_matrix(4, 16) = expected_matrix(1, 13);
    expected_matrix(4, 19) = expected_matrix(1, 16);
    expected_matrix(4, 22) = expected_matrix(1, 19);

    expected_matrix(5, 5)  = expected_matrix(2, 2);
    expected_matrix(5, 8)  = expected_matrix(2, 5);
    expected_matrix(5, 11) = expected_matrix(2, 8);
    expected_matrix(5, 14) = expected_matrix(2, 14);
    expected_matrix(5, 17) = expected_matrix(2, 14);
    expected_matrix(5, 20) = expected_matrix(2, 17);
    expected_matrix(5, 23) = expected_matrix(2, 20);

    expected_matrix(6, 6)  = expected_matrix(0, 0);
    expected_matrix(6, 9)  = expected_matrix(0, 3);
    expected_matrix(6, 12) = expected_matrix(0, 15);
    expected_matrix(6, 15) = expected_matrix(0, 12);
    expected_matrix(6, 18) = expected_matrix(0, 12);
    expected_matrix(6, 21) = expected_matrix(0, 18);

    expected_matrix(7, 7)  = expected_matrix(1, 1);
    expected_matrix(7, 10) = expected_matrix(1, 4);
    expected_matrix(7, 13) = expected_matrix(1, 16);
    expected_matrix(7, 16) = expected_matrix(1, 13);
    expected_matrix(7, 19) = expected_matrix(1, 13);
    expected_matrix(7, 22) = expected_matrix(1, 19);

    expected_matrix(8, 8)  = expected_matrix(2, 2);
    expected_matrix(8, 11) = expected_matrix(2, 5);
    expected_matrix(8, 14) = expected_matrix(2, 17);
    expected_matrix(8, 17) = expected_matrix(2, 14);
    expected_matrix(8, 20) = expected_matrix(2, 14);
    expected_matrix(8, 23) = expected_matrix(2, 20);

    expected_matrix(9, 9)  = expected_matrix(0, 0);
    expected_matrix(9, 12) = expected_matrix(0, 15);
    expected_matrix(9, 15) = expected_matrix(0, 15);
    expected_matrix(9, 18) = expected_matrix(0, 12);
    expected_matrix(9, 21) = expected_matrix(0, 12);

    expected_matrix(10, 10) = expected_matrix(1, 1);
    expected_matrix(10, 13) = expected_matrix(1, 16);
    expected_matrix(10, 16) = expected_matrix(1, 16);
    expected_matrix(10, 19) = expected_matrix(1, 13);
    expected_matrix(10, 22) = expected_matrix(1, 13);

    expected_matrix(11, 11) = expected_matrix(2, 2);
    expected_matrix(11, 14) = expected_matrix(2, 17);
    expected_matrix(11, 17) = expected_matrix(2, 17);
    expected_matrix(11, 20) = expected_matrix(2, 14);
    expected_matrix(11, 23) = expected_matrix(2, 14);

    expected_matrix(12, 12) = shear_stiffness_over_virt_thickness * 8 / 45;
    expected_matrix(12, 15) = shear_stiffness_over_virt_thickness / 9;
    expected_matrix(12, 18) = shear_stiffness_over_virt_thickness * 4 / 45;
    expected_matrix(12, 21) = expected_matrix(12, 15);

    expected_matrix(13, 13) = shear_stiffness_over_virt_thickness * 8 / 45;
    expected_matrix(13, 16) = shear_stiffness_over_virt_thickness / 9;
    expected_matrix(13, 19) = shear_stiffness_over_virt_thickness * 4 / 45;
    expected_matrix(13, 22) = expected_matrix(13, 16);

    expected_matrix(14, 14) = confined_stiffness_over_virt_thickness * 8 / 45;
    expected_matrix(14, 17) = confined_stiffness_over_virt_thickness / 9;
    expected_matrix(14, 20) = confined_stiffness_over_virt_thickness * 4 / 45;
    expected_matrix(14, 23) = expected_matrix(14, 17);

    expected_matrix(15, 15) = expected_matrix(12, 12);
    expected_matrix(15, 18) = expected_matrix(12, 15);
    expected_matrix(15, 21) = expected_matrix(12, 18);

    expected_matrix(16, 16) = expected_matrix(13, 13);
    expected_matrix(16, 19) = expected_matrix(13, 16);
    expected_matrix(16, 22) = expected_matrix(13, 19);

    expected_matrix(17, 17) = expected_matrix(14, 14);
    expected_matrix(17, 20) = expected_matrix(14, 17);
    expected_matrix(17, 23) = expected_matrix(14, 20);

    expected_matrix(18, 18) = expected_matrix(12, 12);
    expected_matrix(18, 21) = expected_matrix(12, 15);

    expected_matrix(19, 19) = expected_matrix(13, 13);
    expected_matrix(19, 22) = expected_matrix(13, 16);

    expected_matrix(20, 20) = expected_matrix(14, 14);
    expected_matrix(20, 23) = expected_matrix(14, 17);

    expected_matrix(21, 21) = expected_matrix(12, 12);

    expected_matrix(22, 22) = expected_matrix(13, 13);

    expected_matrix(23, 23) = expected_matrix(14, 14);

    SymmetrizeMatrix(expected_matrix);

    Vector expected_rhs = ZeroVector(condition_size);
    expected_rhs[0] =
        -(expected_matrix(0, 0) * 1.0 + expected_matrix(0, 3) * 2.0 + expected_matrix(0, 6) * 3.0 +
          expected_matrix(0, 9) * 4.0 + expected_matrix(0, 12) * 1.5 + expected_matrix(0, 15) * 2.5 +
          expected_matrix(0, 18) * 3.5 + expected_matrix(0, 21) * 2.5);
    expected_rhs[1] =
        -(expected_matrix(1, 1) * 2.0 + expected_matrix(1, 4) * 3.0 + expected_matrix(1, 7) * 4.0 +
          expected_matrix(1, 10) * 5.0 + expected_matrix(1, 13) * 2.5 + expected_matrix(1, 16) * 3.5 +
          expected_matrix(1, 19) * 4.5 + expected_matrix(1, 22) * 3.5);
    expected_rhs[2] =
        -(expected_matrix(2, 2) * 3.0 + expected_matrix(2, 5) * 4.0 + expected_matrix(2, 8) * 5.0 +
          expected_matrix(2, 11) * 6.0 + expected_matrix(2, 14) * 3.5 + expected_matrix(2, 17) * 4.5 +
          expected_matrix(2, 20) * 5.5 + expected_matrix(2, 23) * 4.5);
    expected_rhs[3] =
        -(expected_matrix(3, 0) * 1.0 + expected_matrix(3, 3) * 2.0 + expected_matrix(3, 6) * 3.0 +
          expected_matrix(3, 9) * 4.0 + expected_matrix(3, 12) * 1.5 + expected_matrix(3, 15) * 2.5 +
          expected_matrix(3, 18) * 3.5 + expected_matrix(3, 21) * 2.5);
    expected_rhs[4] =
        -(expected_matrix(4, 1) * 2.0 + expected_matrix(4, 4) * 3.0 + expected_matrix(4, 7) * 4.0 +
          expected_matrix(4, 10) * 5.0 + expected_matrix(4, 13) * 2.5 + expected_matrix(4, 16) * 3.5 +
          expected_matrix(4, 19) * 4.5 + expected_matrix(4, 22) * 3.5);
    expected_rhs[5] =
        -(expected_matrix(5, 2) * 3.0 + expected_matrix(5, 5) * 4.0 + expected_matrix(5, 8) * 5.0 +
          expected_matrix(5, 11) * 6.0 + expected_matrix(5, 14) * 3.5 + expected_matrix(5, 17) * 4.5 +
          expected_matrix(5, 20) * 5.5 + expected_matrix(5, 23) * 4.5);
    expected_rhs[6] =
        -(expected_matrix(6, 0) * 1.0 + expected_matrix(6, 3) * 2.0 + expected_matrix(6, 6) * 3.0 +
          expected_matrix(6, 9) * 4.0 + expected_matrix(6, 12) * 1.5 + expected_matrix(6, 15) * 2.5 +
          expected_matrix(6, 18) * 3.5 + expected_matrix(6, 21) * 2.5);
    expected_rhs[7] =
        -(expected_matrix(7, 1) * 2.0 + expected_matrix(7, 4) * 3.0 + expected_matrix(7, 7) * 4.0 +
          expected_matrix(7, 10) * 5.0 + expected_matrix(7, 13) * 2.5 + expected_matrix(7, 16) * 3.5 +
          expected_matrix(7, 19) * 4.5 + expected_matrix(7, 22) * 3.5);
    expected_rhs[8] =
        -(expected_matrix(8, 2) * 3.0 + expected_matrix(8, 5) * 4.0 + expected_matrix(8, 8) * 5.0 +
          expected_matrix(8, 11) * 6.0 + expected_matrix(8, 14) * 3.5 + expected_matrix(8, 17) * 4.5 +
          expected_matrix(8, 20) * 5.5 + expected_matrix(8, 23) * 4.5);
    expected_rhs[9] =
        -(expected_matrix(9, 0) * 1.0 + expected_matrix(9, 3) * 2.0 + expected_matrix(9, 6) * 3.0 +
          expected_matrix(9, 9) * 4.0 + expected_matrix(9, 12) * 1.5 + expected_matrix(9, 15) * 2.5 +
          expected_matrix(9, 18) * 3.5 + expected_matrix(9, 21) * 2.5);
    expected_rhs[10] =
        -(expected_matrix(10, 1) * 2.0 + expected_matrix(10, 4) * 3.0 + expected_matrix(10, 7) * 4.0 +
          expected_matrix(10, 10) * 5.0 + expected_matrix(10, 13) * 2.5 + expected_matrix(10, 16) * 3.5 +
          expected_matrix(10, 19) * 4.5 + expected_matrix(10, 22) * 3.5);
    expected_rhs[11] =
        -(expected_matrix(11, 2) * 3.0 + expected_matrix(11, 5) * 4.0 + expected_matrix(11, 8) * 5.0 +
          expected_matrix(11, 11) * 6.0 + expected_matrix(11, 14) * 3.5 + expected_matrix(11, 17) * 4.5 +
          expected_matrix(11, 20) * 5.5 + expected_matrix(11, 23) * 4.5);
    expected_rhs[12] =
        -(expected_matrix(12, 0) * 1.0 + expected_matrix(12, 3) * 2.0 + expected_matrix(12, 6) * 3.0 +
          expected_matrix(12, 9) * 4.0 + expected_matrix(12, 12) * 1.5 + expected_matrix(12, 15) * 2.5 +
          expected_matrix(12, 18) * 3.5 + expected_matrix(12, 21) * 2.5);
    expected_rhs[13] =
        -(expected_matrix(13, 1) * 2.0 + expected_matrix(13, 4) * 3.0 + expected_matrix(13, 7) * 4.0 +
          expected_matrix(13, 10) * 5.0 + expected_matrix(13, 13) * 2.5 + expected_matrix(13, 16) * 3.5 +
          expected_matrix(13, 19) * 4.5 + expected_matrix(13, 22) * 3.5);
    expected_rhs[14] =
        -(expected_matrix(14, 2) * 3.0 + expected_matrix(14, 5) * 4.0 + expected_matrix(14, 8) * 5.0 +
          expected_matrix(14, 11) * 6.0 + expected_matrix(14, 14) * 3.5 + expected_matrix(14, 17) * 4.5 +
          expected_matrix(14, 20) * 5.5 + expected_matrix(14, 23) * 4.5);
    expected_rhs[15] =
        -(expected_matrix(15, 0) * 1.0 + expected_matrix(15, 3) * 2.0 + expected_matrix(15, 6) * 3.0 +
          expected_matrix(15, 9) * 4.0 + expected_matrix(15, 12) * 1.5 + expected_matrix(15, 15) * 2.5 +
          expected_matrix(15, 18) * 3.5 + expected_matrix(15, 21) * 2.5);
    expected_rhs[16] =
        -(expected_matrix(16, 1) * 2.0 + expected_matrix(16, 4) * 3.0 + expected_matrix(16, 7) * 4.0 +
          expected_matrix(16, 10) * 5.0 + expected_matrix(16, 13) * 2.5 + expected_matrix(16, 16) * 3.5 +
          expected_matrix(16, 19) * 4.5 + expected_matrix(16, 22) * 3.5);
    expected_rhs[17] =
        -(expected_matrix(17, 2) * 3.0 + expected_matrix(17, 5) * 4.0 + expected_matrix(17, 8) * 5.0 +
          expected_matrix(17, 11) * 6.0 + expected_matrix(17, 14) * 3.5 + expected_matrix(17, 17) * 4.5 +
          expected_matrix(17, 20) * 5.5 + expected_matrix(17, 23) * 4.5);
    expected_rhs[18] =
        -(expected_matrix(18, 0) * 1.0 + expected_matrix(18, 3) * 2.0 + expected_matrix(18, 6) * 3.0 +
          expected_matrix(18, 9) * 4.0 + expected_matrix(18, 12) * 1.5 + expected_matrix(18, 15) * 2.5 +
          expected_matrix(18, 18) * 3.5 + expected_matrix(18, 21) * 2.5);
    expected_rhs[19] =
        -(expected_matrix(19, 1) * 2.0 + expected_matrix(19, 4) * 3.0 + expected_matrix(19, 7) * 4.0 +
          expected_matrix(19, 10) * 5.0 + expected_matrix(19, 13) * 2.5 + expected_matrix(19, 16) * 3.5 +
          expected_matrix(19, 19) * 4.5 + expected_matrix(19, 22) * 3.5);
    expected_rhs[20] =
        -(expected_matrix(20, 2) * 3.0 + expected_matrix(20, 5) * 4.0 + expected_matrix(20, 8) * 5.0 +
          expected_matrix(20, 11) * 6.0 + expected_matrix(20, 14) * 3.5 + expected_matrix(20, 17) * 4.5 +
          expected_matrix(20, 20) * 5.5 + expected_matrix(20, 23) * 4.5);
    expected_rhs[21] =
        -(expected_matrix(21, 0) * 1.0 + expected_matrix(21, 3) * 2.0 + expected_matrix(21, 6) * 3.0 +
          expected_matrix(21, 9) * 4.0 + expected_matrix(21, 12) * 1.5 + expected_matrix(21, 15) * 2.5 +
          expected_matrix(21, 18) * 3.5 + expected_matrix(21, 21) * 2.5);
    expected_rhs[22] =
        -(expected_matrix(22, 1) * 2.0 + expected_matrix(22, 4) * 3.0 + expected_matrix(22, 7) * 4.0 +
          expected_matrix(22, 10) * 5.0 + expected_matrix(22, 13) * 2.5 + expected_matrix(22, 16) * 3.5 +
          expected_matrix(22, 19) * 4.5 + expected_matrix(22, 22) * 3.5);
    expected_rhs[23] =
        -(expected_matrix(23, 2) * 3.0 + expected_matrix(23, 5) * 4.0 + expected_matrix(23, 8) * 5.0 +
          expected_matrix(23, 11) * 6.0 + expected_matrix(23, 14) * 3.5 + expected_matrix(23, 17) * 4.5 +
          expected_matrix(23, 20) * 5.5 + expected_matrix(23, 23) * 4.5);

    // compare results
    KRATOS_EXPECT_MATRIX_NEAR(rLeftHandSideMatrix, expected_matrix, 1.0e-6);
    KRATOS_EXPECT_VECTOR_NEAR(right_hand_side_vector, expected_rhs, 1.0e-6);
}

/// <summary>
/// Tests the calculation of the damping matrix of a 3D 8N UPwNormalLysmerAbsorbingCondition
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(CalculateDampingMatrixUPwNormalLysmerAbsorbingCondition3D8N, KratosGeoMechanicsFastSuite)
{
    GTEST_SKIP() << "This test is skipped as it depends on #13269" << std::endl;

    // initialize modelpart
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    ModelPart::ConditionType::Pointer p_cond = SetUpUPwLysmerAbsorbingCondition3D8NCondition(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Perform test, calculate left hand side
    constexpr size_t condition_size      = 32;
    Matrix           rLeftHandSideMatrix = ZeroMatrix(condition_size, condition_size);

    p_cond->CalculateDampingMatrix(rLeftHandSideMatrix, r_process_info);

    // set expected_results
    const double vp = sqrt(4.0 / 2000); // sqrt(confined_stiffness / density)
    const double vs = sqrt(5.0 / 2000); // sqrt(shear_stiffness / density)

    const double perpendicular_damping = vp * 2000 * 1; // vp * density * perpendicular_damping_factor
    const double shear_damping = vs * 2000 * 1;         // vs * density * shear_damping_factor

    Matrix expected_matrix = ZeroMatrix(condition_size, condition_size);
    expected_matrix(0, 0)  = shear_damping / 30;
    expected_matrix(0, 3)  = shear_damping / 90;
    expected_matrix(0, 6)  = shear_damping / 60;
    expected_matrix(0, 9)  = expected_matrix(0, 3);
    expected_matrix(0, 12) = -shear_damping / 30;
    expected_matrix(0, 15) = -shear_damping * 2 / 45;
    expected_matrix(0, 18) = expected_matrix(0, 15);
    expected_matrix(0, 21) = expected_matrix(0, 12);

    expected_matrix(1, 1)  = shear_damping / 30;
    expected_matrix(1, 4)  = shear_damping / 90;
    expected_matrix(1, 7)  = shear_damping / 60;
    expected_matrix(1, 10) = expected_matrix(1, 4);
    expected_matrix(1, 13) = -shear_damping / 30;
    expected_matrix(1, 16) = -shear_damping * 2 / 45;
    expected_matrix(1, 19) = expected_matrix(1, 16);
    expected_matrix(1, 22) = expected_matrix(1, 13);

    expected_matrix(2, 2)  = perpendicular_damping / 30;
    expected_matrix(2, 5)  = perpendicular_damping / 90;
    expected_matrix(2, 8)  = perpendicular_damping / 60;
    expected_matrix(2, 11) = expected_matrix(2, 5);
    expected_matrix(2, 14) = -perpendicular_damping / 30;
    expected_matrix(2, 17) = -perpendicular_damping * 2 / 45;
    expected_matrix(2, 20) = expected_matrix(2, 17);
    expected_matrix(2, 23) = expected_matrix(2, 14);

    expected_matrix(3, 3)  = expected_matrix(0, 0);
    expected_matrix(3, 6)  = expected_matrix(0, 3);
    expected_matrix(3, 9)  = expected_matrix(0, 6);
    expected_matrix(3, 12) = expected_matrix(0, 12);
    expected_matrix(3, 15) = expected_matrix(0, 12);
    expected_matrix(3, 18) = expected_matrix(0, 15);
    expected_matrix(3, 21) = expected_matrix(0, 18);

    expected_matrix(4, 4)  = expected_matrix(1, 1);
    expected_matrix(4, 7)  = expected_matrix(1, 4);
    expected_matrix(4, 10) = expected_matrix(1, 7);
    expected_matrix(4, 13) = expected_matrix(1, 13);
    expected_matrix(4, 16) = expected_matrix(1, 13);
    expected_matrix(4, 19) = expected_matrix(1, 16);
    expected_matrix(4, 22) = expected_matrix(1, 19);

    expected_matrix(5, 5)  = expected_matrix(2, 2);
    expected_matrix(5, 8)  = expected_matrix(2, 5);
    expected_matrix(5, 11) = expected_matrix(2, 8);
    expected_matrix(5, 14) = expected_matrix(2, 14);
    expected_matrix(5, 17) = expected_matrix(2, 14);
    expected_matrix(5, 20) = expected_matrix(2, 17);
    expected_matrix(5, 23) = expected_matrix(2, 20);

    expected_matrix(6, 6)  = expected_matrix(0, 0);
    expected_matrix(6, 9)  = expected_matrix(0, 3);
    expected_matrix(6, 12) = expected_matrix(0, 15);
    expected_matrix(6, 15) = expected_matrix(0, 12);
    expected_matrix(6, 18) = expected_matrix(0, 12);
    expected_matrix(6, 21) = expected_matrix(0, 18);

    expected_matrix(7, 7)  = expected_matrix(1, 1);
    expected_matrix(7, 10) = expected_matrix(1, 4);
    expected_matrix(7, 13) = expected_matrix(1, 16);
    expected_matrix(7, 16) = expected_matrix(1, 13);
    expected_matrix(7, 19) = expected_matrix(1, 13);
    expected_matrix(7, 22) = expected_matrix(1, 19);

    expected_matrix(8, 8)  = expected_matrix(2, 2);
    expected_matrix(8, 11) = expected_matrix(2, 5);
    expected_matrix(8, 14) = expected_matrix(2, 17);
    expected_matrix(8, 17) = expected_matrix(2, 14);
    expected_matrix(8, 20) = expected_matrix(2, 14);
    expected_matrix(8, 23) = expected_matrix(2, 20);

    expected_matrix(9, 9)  = expected_matrix(0, 0);
    expected_matrix(9, 12) = expected_matrix(0, 15);
    expected_matrix(9, 15) = expected_matrix(0, 15);
    expected_matrix(9, 18) = expected_matrix(0, 12);
    expected_matrix(9, 21) = expected_matrix(0, 12);

    expected_matrix(10, 10) = expected_matrix(1, 1);
    expected_matrix(10, 13) = expected_matrix(1, 16);
    expected_matrix(10, 16) = expected_matrix(1, 16);
    expected_matrix(10, 19) = expected_matrix(1, 13);
    expected_matrix(10, 22) = expected_matrix(1, 13);

    expected_matrix(11, 11) = expected_matrix(2, 2);
    expected_matrix(11, 14) = expected_matrix(2, 17);
    expected_matrix(11, 17) = expected_matrix(2, 17);
    expected_matrix(11, 20) = expected_matrix(2, 14);
    expected_matrix(11, 23) = expected_matrix(2, 14);

    expected_matrix(12, 12) = shear_damping * 8 / 45;
    expected_matrix(12, 15) = shear_damping / 9;
    expected_matrix(12, 18) = shear_damping * 4 / 45;
    expected_matrix(12, 21) = expected_matrix(12, 15);

    expected_matrix(13, 13) = shear_damping * 8 / 45;
    expected_matrix(13, 16) = shear_damping / 9;
    expected_matrix(13, 19) = shear_damping * 4 / 45;
    expected_matrix(13, 22) = expected_matrix(13, 16);

    expected_matrix(14, 14) = perpendicular_damping * 8 / 45;
    expected_matrix(14, 17) = perpendicular_damping / 9;
    expected_matrix(14, 20) = perpendicular_damping * 4 / 45;
    expected_matrix(14, 23) = expected_matrix(14, 17);

    expected_matrix(15, 15) = expected_matrix(12, 12);
    expected_matrix(15, 18) = expected_matrix(12, 15);
    expected_matrix(15, 21) = expected_matrix(12, 18);

    expected_matrix(16, 16) = expected_matrix(13, 13);
    expected_matrix(16, 19) = expected_matrix(13, 16);
    expected_matrix(16, 22) = expected_matrix(13, 19);

    expected_matrix(17, 17) = expected_matrix(14, 14);
    expected_matrix(17, 20) = expected_matrix(14, 17);
    expected_matrix(17, 23) = expected_matrix(14, 20);

    expected_matrix(18, 18) = expected_matrix(12, 12);
    expected_matrix(18, 21) = expected_matrix(12, 15);

    expected_matrix(19, 19) = expected_matrix(13, 13);
    expected_matrix(19, 22) = expected_matrix(13, 16);

    expected_matrix(20, 20) = expected_matrix(14, 14);
    expected_matrix(20, 23) = expected_matrix(14, 17);

    expected_matrix(21, 21) = expected_matrix(12, 12);

    expected_matrix(22, 22) = expected_matrix(13, 13);

    expected_matrix(23, 23) = expected_matrix(14, 14);

    SymmetrizeMatrix(expected_matrix);

    // compare results
    KRATOS_EXPECT_MATRIX_NEAR(rLeftHandSideMatrix, expected_matrix, 1.0e-6);
}

} // namespace Kratos::Testing
