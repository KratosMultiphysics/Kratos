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
#include "includes/element.h"
#include "testing/testing.h"

namespace
{

using namespace Kratos;

ModelPart& CreateTestModelPart(Model& rModel, const Variable<double>& rVariable)
{
    const auto buffer_size = std::size_t{2};
    auto&      r_result    = rModel.CreateModelPart("Dummy", buffer_size);
    r_result.AddNodalSolutionStepVariable(rVariable);
    return r_result;
}

Dof<double> MakeDofWithEquationId(Dof<double>::EquationIdType EquationId)
{
    Dof<double> result;
    result.SetEquationId(EquationId);
    return result;
}

intrusive_ptr<Node> AddNodeWithDof(ModelPart&              rModelPart,
                                   ModelPart::IndexType    NodeId,
                                   double                  x,
                                   double                  y,
                                   double                  z,
                                   const Variable<double>& rDofVariable)
{
    auto p_result = rModelPart.CreateNewNode(NodeId, x, y, z);
    p_result->AddDof(rDofVariable);
    return p_result;
}

void AddThreeNodesWithDofs(ModelPart& rModelPart, const Variable<double>& rVariable)
{
    AddNodeWithDof(rModelPart, 1, 0.0, 0.0, 0.0, rVariable);
    AddNodeWithDof(rModelPart, 2, 1.0, 0.0, 0.0, rVariable);
    AddNodeWithDof(rModelPart, 3, 0.0, 1.0, 0.0, rVariable);
}

using NodeIndexAndValue = std::pair<std::size_t, double>;

void SetNodalValues(ModelPart&                            rModelPart,
                    const std::vector<NodeIndexAndValue>& rNodalValues,
                    const Variable<double>&               rDofVariable,
                    std::size_t                           BufferIndex)
{
    for (const auto& [index, value] : rNodalValues) {
        auto& r_node                                               = rModelPart.GetNode(index);
        r_node.FastGetSolutionStepValue(rDofVariable, BufferIndex) = value;
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

KRATOS_TEST_CASE_IN_SUITE(ExtractingValuesFromDofsGetsNodalValues, KratosGeoMechanicsFastSuite)
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

} // namespace Kratos::Testing
