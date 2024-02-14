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

#include "custom_utilities/element_utilities.hpp"
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
    Element::DofsVectorType dofs;

    KRATOS_EXPECT_TRUE(GeoElementUtilities::ExtractEquationIdsFrom(dofs).empty())
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingEquationIdsFromDofsYieldsAssociatedIds, KratosGeoMechanicsFastSuite)
{
    auto       dof1 = MakeDofWithEquationId(22);
    auto       dof2 = MakeDofWithEquationId(20);
    auto       dof3 = MakeDofWithEquationId(21);
    const auto dofs = Element::DofsVectorType{&dof1, &dof2, &dof3};

    const auto expected_equation_ids = Element::EquationIdVectorType{22, 20, 21};
    KRATOS_EXPECT_EQ(GeoElementUtilities::ExtractEquationIdsFrom(dofs), expected_equation_ids);
}

} // namespace Kratos::Testing
