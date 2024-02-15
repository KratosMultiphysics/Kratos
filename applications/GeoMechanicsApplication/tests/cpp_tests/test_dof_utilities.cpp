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

    KRATOS_EXPECT_TRUE(ExtractEquationIdsFrom(dofs).empty())
}

KRATOS_TEST_CASE_IN_SUITE(ExtractingEquationIdsFromDofsYieldsAssociatedIds, KratosGeoMechanicsFastSuite)
{
    auto       dof1 = MakeDofWithEquationId(22);
    auto       dof2 = MakeDofWithEquationId(20);
    auto       dof3 = MakeDofWithEquationId(21);
    const auto dofs = std::vector<Dof<double>*>{&dof1, &dof2, &dof3};

    const auto expected_equation_ids = std::vector<std::size_t>{22, 20, 21};
    KRATOS_EXPECT_EQ(ExtractEquationIdsFrom(dofs), expected_equation_ids);
}

} // namespace Kratos::Testing
