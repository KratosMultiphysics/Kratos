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

#include "custom_processes/apply_seepage_boundary_process.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplySeepageBoundaryProcess_CanBeDefaultConstructed, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    ApplySeepageBoundaryProcess process;
}

} // namespace Kratos::Testing