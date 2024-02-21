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

#include "containers/model.h"
#include "custom_utilities/dof_utilities.h"
#include "includes/element.h"
#include "testing/testing.h"

namespace
{

using namespace Kratos;

Dof<double> MakeDofWithEquationId(Dof<double>::EquationIdType EquationId)
{
    Dof<double> result;
    result.SetEquationId(EquationId);
    return result;
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

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(Geo::DofUtilities::ExtractDofsFromNodes(p_element->GetGeometry(), DISPLACEMENT_X),
                                      "Non-existent DOF in node #1 for variable : DISPLACEMENT_X")
}

} // namespace Kratos::Testing
