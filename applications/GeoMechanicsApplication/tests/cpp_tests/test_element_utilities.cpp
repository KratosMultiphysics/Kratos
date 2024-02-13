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

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ExtractingEquationIdsFromEmptyDofListReturnsEmptyList, KratosGeoMechanicsFastSuite)
{
    Element::DofsVectorType dofs;

    KRATOS_EXPECT_TRUE(GeoElementUtilities::ExtractEquationIdsFrom(dofs).empty());
}

} // namespace Kratos::Testing
