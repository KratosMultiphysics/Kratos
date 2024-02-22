// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Richard Faasse
//

#include <boost/numeric/ublas/assignment.hpp>
#include <boost/range/algorithm/copy.hpp>

#include "containers/model.h"
#include "custom_utilities/dof_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/element.h"
#include "testing/testing.h"

namespace
{

using namespace Kratos;

using ConstVariableRefs = std::vector<std::reference_wrapper<const Variable<double>>>;

ModelPart& CreateTestModelPart(Model& rModel, const ConstVariableRefs& rVariables)
{
    const auto buffer_size = std::size_t{2};
    auto&      r_result    = rModel.CreateModelPart("Dummy", buffer_size);
    for (const auto& r_variable : rVariables) {
        r_result.AddNodalSolutionStepVariable(r_variable);
    }
    return r_result;
}

ModelPart& CreateTestModelPart(Model& rModel, const Variable<double>& rVariable)
{
    const auto variables = ConstVariableRefs{std::cref(rVariable)};
    return CreateTestModelPart(rModel, variables);
}

Dof<double> MakeDofWithEquationId(Dof<double>::EquationIdType EquationId)
{
    Dof<double> result;
    result.SetEquationId(EquationId);
    return result;
}

[[maybe_unused]] intrusive_ptr<Node> AddNodeWithDofs(
    ModelPart& rModelPart, ModelPart::IndexType NodeId, double x, double y, double z, const ConstVariableRefs& rVariables)
{
    auto p_result = rModelPart.CreateNewNode(NodeId, x, y, z);
    for (const auto& r_variable : rVariables) {
        p_result->AddDof(r_variable.get());
    }
    return p_result;
}

[[maybe_unused]] intrusive_ptr<Node> AddNodeWithDof(
    ModelPart& rModelPart, ModelPart::IndexType NodeId, double x, double y, double z, const Variable<double>& rVariable)
{
    const auto variables = ConstVariableRefs{std::cref(rVariable)};
    return AddNodeWithDofs(rModelPart, NodeId, x, y, z, variables);
}

void AddThreeNodesWithDofs(ModelPart& rModelPart, const ConstVariableRefs& rVariables)
{
    AddNodeWithDofs(rModelPart, 1, 0.0, 0.0, 0.0, rVariables);
    AddNodeWithDofs(rModelPart, 2, 1.0, 0.0, 0.0, rVariables);
    AddNodeWithDofs(rModelPart, 3, 0.0, 1.0, 0.0, rVariables);
}

void AddThreeNodesWithDofs(ModelPart& rModelPart, const Variable<double>& rVariable)
{
    const auto variables = ConstVariableRefs{std::cref(rVariable)};
    AddThreeNodesWithDofs(rModelPart, variables);
}

using NodeIndexAndValue = std::pair<std::size_t, double>;

void SetNodalValues(ModelPart&                            rModelPart,
                    const std::vector<NodeIndexAndValue>& rNodalValues,
                    const Variable<double>&               rVariable,
                    std::size_t                           BufferIndex)
{
    for (const auto& [index, value] : rNodalValues) {
        auto& r_node                                            = rModelPart.GetNode(index);
        r_node.FastGetSolutionStepValue(rVariable, BufferIndex) = value;
    }
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ExtractingEquationIdsFromEmptyDofListReturnsEmptyList, KratosGeoMechanicsFastSuite)
{
    std::vector<Dof<double>*> dofs;

    KRATOS_EXPECT_TRUE(Geo::DofUtilities::ExtractEquationIdsFrom(dofs).empty())
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingEquationIdsFromDofsYieldsAssociatedIds, KratosGeoMechanicsFastSuite)
{
    auto       dof1 = MakeDofWithEquationId(22);
    auto       dof2 = MakeDofWithEquationId(20);
    auto       dof3 = MakeDofWithEquationId(21);
    const auto dofs = std::vector<Dof<double>*>{&dof1, &dof2, &dof3};

    const auto expected_equation_ids = std::vector<std::size_t>{22, 20, 21};
    KRATOS_EXPECT_EQ(Geo::DofUtilities::ExtractEquationIdsFrom(dofs), expected_equation_ids);
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingDofsFromEmptyNodeCollectionYieldsEmptyVector, KratosGeoMechanicsFastSuite)
{
    const auto geometry_id    = Geometry<Node>::IndexType{0};
    const auto empty_geometry = Geometry<Node>{geometry_id};

    KRATOS_EXPECT_TRUE(Geo::DofUtilities::ExtractDofsFromNodes(empty_geometry, WATER_PRESSURE).empty())
}

KRATOS_TEST_CASE_IN_SUITE(ExpectThrowWhenExtractingNonExistingDofsFromNodes, KratosGeoMechanicsFastSuite)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Dummy");

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

    const auto node_ids  = std::vector<ModelPart::IndexType>{1, 2, 3};
    const auto p_element = r_model_part.CreateNewElement("UPwSmallStrainElement2D3N", 1, node_ids,
                                                         r_model_part.CreateNewProperties(0));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        Geo::DofUtilities::ExtractDofsFromNodes(p_element->GetGeometry(), DISPLACEMENT_X),
        "Non-existent DOF in node #1 for variable : DISPLACEMENT_X")
}

KRATOS_TEST_CASE_IN_SUITE(VariableTypeAndNodeIDsMustMatchWhenExtractingDofsFromNodes, KratosGeoMechanicsFastSuite)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Dummy");
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);

    AddNodeWithDof(r_model_part, 1, 0.0, 0.0, 0.0, DISPLACEMENT_X);
    AddNodeWithDof(r_model_part, 2, 1.0, 0.0, 0.0, DISPLACEMENT_X);
    AddNodeWithDof(r_model_part, 3, 0.0, 1.0, 0.0, DISPLACEMENT_X);

    const auto node_ids  = std::vector<ModelPart::IndexType>{1, 2, 3};
    const auto p_element = r_model_part.CreateNewElement("UPwSmallStrainElement2D3N", 1, node_ids,
                                                         r_model_part.CreateNewProperties(0));

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(p_element->GetGeometry(), DISPLACEMENT_X);

    KRATOS_EXPECT_EQ(dofs.size(), node_ids.size());
    KRATOS_EXPECT_TRUE(
        std::all_of(dofs.begin(), dofs.end(), [](const auto p_dof) { return p_dof != nullptr; }))
    KRATOS_EXPECT_TRUE(std::all_of(dofs.begin(), dofs.end(), [](const auto p_dof) {
        return p_dof->GetVariable() == DISPLACEMENT_X;
    }))
    KRATOS_EXPECT_EQ(dofs[0]->GetId(), 1);
    KRATOS_EXPECT_EQ(dofs[1]->GetId(), 2);
    KRATOS_EXPECT_EQ(dofs[2]->GetId(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingValuesFromDofsYieldsNodalValues, KratosGeoMechanicsFastSuite)
{
    auto        model        = Model{};
    const auto& r_variable   = DISPLACEMENT_X;
    auto&       r_model_part = CreateTestModelPart(model, r_variable);
    AddThreeNodesWithDofs(r_model_part, r_variable);

    const auto current_buffer_index = std::size_t{0};
    const auto current_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 1.0), std::make_pair(2, 2.0), std::make_pair(3, 3.0)};
    SetNodalValues(r_model_part, current_values, r_variable, current_buffer_index);

    const auto previous_buffer_index = std::size_t{1};
    const auto previous_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 4.0), std::make_pair(2, 5.0), std::make_pair(3, 6.0)};
    SetNodalValues(r_model_part, previous_values, r_variable, previous_buffer_index);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    auto expected_values = Vector(3);
    expected_values <<= 1.0, 2.0, 3.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSolutionStepValues(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    boost::range::copy(std::vector<double>{4.0, 5.0, 6.0}, expected_values.begin());
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSolutionStepValues(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingFirstDerivativeValuesFromDofsYieldsNodalFirstDerivativeValues, KratosGeoMechanicsFastSuite)
{
    const auto& r_variable              = DISPLACEMENT_X;
    const auto& r_first_time_derivative = VELOCITY_X;
    const auto all_variables = ConstVariableRefs{std::cref(r_variable), std::cref(r_first_time_derivative)};

    auto  model        = Model{};
    auto& r_model_part = CreateTestModelPart(model, all_variables);
    AddThreeNodesWithDofs(r_model_part, all_variables);

    const auto current_buffer_index = std::size_t{0};
    const auto current_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 1.0), std::make_pair(2, 2.0), std::make_pair(3, 3.0)};
    SetNodalValues(r_model_part, current_values, r_first_time_derivative, current_buffer_index);

    const auto previous_buffer_index = std::size_t{1};
    const auto previous_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 4.0), std::make_pair(2, 5.0), std::make_pair(3, 6.0)};
    SetNodalValues(r_model_part, previous_values, r_first_time_derivative, previous_buffer_index);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    auto expected_values = Vector(3);
    expected_values <<= 1.0, 2.0, 3.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractFirstTimeDerivatives(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    boost::range::copy(std::vector<double>{4.0, 5.0, 6.0}, expected_values.begin());
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractFirstTimeDerivatives(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingSecondDerivativeValuesFromDofsYieldsNodalSecondDerivativeValues,
                          KratosGeoMechanicsFastSuite)
{
    const auto& r_variable               = DISPLACEMENT_X;
    const auto& r_first_time_derivative  = VELOCITY_X;
    const auto& r_second_time_derivative = ACCELERATION_X;
    const auto  all_variables            = ConstVariableRefs{
        std::cref(r_variable), std::cref(r_first_time_derivative), std::cref(r_second_time_derivative)};

    auto  model        = Model{};
    auto& r_model_part = CreateTestModelPart(model, all_variables);
    AddThreeNodesWithDofs(r_model_part, all_variables);

    const auto current_buffer_index = std::size_t{0};
    const auto current_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 1.0), std::make_pair(2, 2.0), std::make_pair(3, 3.0)};
    SetNodalValues(r_model_part, current_values, r_second_time_derivative, current_buffer_index);

    const auto previous_buffer_index = std::size_t{1};
    const auto previous_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 4.0), std::make_pair(2, 5.0), std::make_pair(3, 6.0)};
    SetNodalValues(r_model_part, previous_values, r_second_time_derivative, previous_buffer_index);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    auto expected_values = Vector(3);
    expected_values <<= 1.0, 2.0, 3.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSecondTimeDerivatives(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    boost::range::copy(std::vector<double>{4.0, 5.0, 6.0}, expected_values.begin());
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSecondTimeDerivatives(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingValuesFromUPwDofsNotBeingPwYieldsNodalValues, KratosGeoMechanicsFastSuite)
{
    auto        model        = Model{};
    const auto& r_variable   = DISPLACEMENT_X;
    auto&       r_model_part = CreateTestModelPart(model, r_variable);
    AddThreeNodesWithDofs(r_model_part, r_variable);

    const auto current_buffer_index = std::size_t{0};
    const auto current_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 1.0), std::make_pair(2, 2.0), std::make_pair(3, 3.0)};
    SetNodalValues(r_model_part, current_values, r_variable, current_buffer_index);

    const auto previous_buffer_index = std::size_t{1};
    const auto previous_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 4.0), std::make_pair(2, 5.0), std::make_pair(3, 6.0)};
    SetNodalValues(r_model_part, previous_values, r_variable, previous_buffer_index);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    auto expected_values = Vector(3);
    expected_values <<= 1.0, 2.0, 3.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSolutionStepValuesOfUPwDofs(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    boost::range::copy(std::vector<double>{4.0, 5.0, 6.0}, expected_values.begin());
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSolutionStepValuesOfUPwDofs(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingPwValuesFromUPwDofsAlwaysYieldsZeroes, KratosGeoMechanicsFastSuite)
{
    auto        model        = Model{};
    const auto& r_variable   = WATER_PRESSURE;
    auto&       r_model_part = CreateTestModelPart(model, r_variable);
    AddThreeNodesWithDofs(r_model_part, r_variable);

    const auto current_buffer_index = std::size_t{0};
    const auto current_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 1.0), std::make_pair(2, 2.0), std::make_pair(3, 3.0)};
    SetNodalValues(r_model_part, current_values, r_variable, current_buffer_index);

    const auto previous_buffer_index = std::size_t{1};
    const auto previous_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 4.0), std::make_pair(2, 5.0), std::make_pair(3, 6.0)};
    SetNodalValues(r_model_part, previous_values, r_variable, previous_buffer_index);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    auto expected_values = Vector(3);
    expected_values <<= 0.0, 0.0, 0.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSolutionStepValuesOfUPwDofs(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSolutionStepValuesOfUPwDofs(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingFirstDerivativesFromUPwDofsNotBeingPwYieldsNodalFirstDerivatives,
                          KratosGeoMechanicsFastSuite)
{
    const auto& r_variable              = DISPLACEMENT_X;
    const auto& r_first_time_derivative = VELOCITY_X;
    const auto all_variables = ConstVariableRefs{std::cref(r_variable), std::cref(r_first_time_derivative)};

    auto  model        = Model{};
    auto& r_model_part = CreateTestModelPart(model, all_variables);
    AddThreeNodesWithDofs(r_model_part, all_variables);

    const auto current_buffer_index = std::size_t{0};
    const auto current_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 1.0), std::make_pair(2, 2.0), std::make_pair(3, 3.0)};
    SetNodalValues(r_model_part, current_values, r_first_time_derivative, current_buffer_index);

    const auto previous_buffer_index = std::size_t{1};
    const auto previous_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 4.0), std::make_pair(2, 5.0), std::make_pair(3, 6.0)};
    SetNodalValues(r_model_part, previous_values, r_first_time_derivative, previous_buffer_index);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    auto expected_values = Vector(3);
    expected_values <<= 1.0, 2.0, 3.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractFirstTimeDerivativesOfUPwDofs(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    boost::range::copy(std::vector<double>{4.0, 5.0, 6.0}, expected_values.begin());
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractFirstTimeDerivativesOfUPwDofs(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingDPwDtValuesFromUPwDofsAlwaysYieldsZeroes, KratosGeoMechanicsFastSuite)
{
    const auto& r_variable              = WATER_PRESSURE;
    const auto& r_first_time_derivative = DT_WATER_PRESSURE;
    const auto all_variables = ConstVariableRefs{std::cref(r_variable), std::cref(r_first_time_derivative)};

    auto  model        = Model{};
    auto& r_model_part = CreateTestModelPart(model, all_variables);
    AddThreeNodesWithDofs(r_model_part, all_variables);

    const auto current_buffer_index = std::size_t{0};
    const auto current_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 1.0), std::make_pair(2, 2.0), std::make_pair(3, 3.0)};
    SetNodalValues(r_model_part, current_values, r_first_time_derivative, current_buffer_index);

    const auto previous_buffer_index = std::size_t{1};
    const auto previous_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 4.0), std::make_pair(2, 5.0), std::make_pair(3, 6.0)};
    SetNodalValues(r_model_part, previous_values, r_first_time_derivative, previous_buffer_index);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    auto expected_values = Vector(3);
    expected_values <<= 0.0, 0.0, 0.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractFirstTimeDerivativesOfUPwDofs(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractFirstTimeDerivativesOfUPwDofs(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingSecondDerivativesFromUPwDofsNotBeingPwYieldsNodalSecondDerivatives,
                          KratosGeoMechanicsFastSuite)
{
    const auto& r_variable               = DISPLACEMENT_X;
    const auto& r_first_time_derivative  = VELOCITY_X;
    const auto& r_second_time_derivative = ACCELERATION_X;
    const auto  all_variables            = ConstVariableRefs{
        std::cref(r_variable), std::cref(r_first_time_derivative), std::cref(r_second_time_derivative)};

    auto  model        = Model{};
    auto& r_model_part = CreateTestModelPart(model, all_variables);
    AddThreeNodesWithDofs(r_model_part, all_variables);

    const auto current_buffer_index = std::size_t{0};
    const auto current_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 1.0), std::make_pair(2, 2.0), std::make_pair(3, 3.0)};
    SetNodalValues(r_model_part, current_values, r_second_time_derivative, current_buffer_index);

    const auto previous_buffer_index = std::size_t{1};
    const auto previous_values       = std::vector<NodeIndexAndValue>{
        std::make_pair(1, 4.0), std::make_pair(2, 5.0), std::make_pair(3, 6.0)};
    SetNodalValues(r_model_part, previous_values, r_second_time_derivative, previous_buffer_index);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    auto expected_values = Vector(3);
    expected_values <<= 1.0, 2.0, 3.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSecondTimeDerivativesOfUPwDofs(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    boost::range::copy(std::vector<double>{4.0, 5.0, 6.0}, expected_values.begin());
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSecondTimeDerivativesOfUPwDofs(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingD2PwDt2ValuesFromUPwDofsAlwaysYieldsZeroes, KratosGeoMechanicsFastSuite)
{
    const auto& r_variable              = WATER_PRESSURE;
    const auto& r_first_time_derivative = DT_WATER_PRESSURE;
    // There is no variable defined for the second time derivative of the water pressure
    const auto all_variables = ConstVariableRefs{std::cref(r_variable), std::cref(r_first_time_derivative)};

    auto  model        = Model{};
    auto& r_model_part = CreateTestModelPart(model, all_variables);
    AddThreeNodesWithDofs(r_model_part, all_variables);

    const auto dofs = Geo::DofUtilities::ExtractDofsFromNodes(r_model_part.Nodes(), r_variable);

    const auto current_buffer_index = std::size_t{0};
    auto       expected_values      = Vector(3);
    expected_values <<= 0.0, 0.0, 0.0;
    const auto abs_tolerance = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSecondTimeDerivativesOfUPwDofs(dofs, current_buffer_index),
                              expected_values, abs_tolerance)

    const auto previous_buffer_index = std::size_t{1};
    KRATOS_EXPECT_VECTOR_NEAR(Geo::DofUtilities::ExtractSecondTimeDerivativesOfUPwDofs(dofs, previous_buffer_index),
                              expected_values, abs_tolerance)
}

} // namespace Kratos::Testing
