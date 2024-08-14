// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Anne van de Graaf
//

#include "../geo_mechanics_fast_suite.h"
#include "custom_geometries/line_interface_geometry.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestConstructInterfaceGeometry, KratosGeoMechanicsFastSuiteWithoutKernel)
{
       LineInterfaceGeometry geometry;
}

} // namespace Kratos::Testing