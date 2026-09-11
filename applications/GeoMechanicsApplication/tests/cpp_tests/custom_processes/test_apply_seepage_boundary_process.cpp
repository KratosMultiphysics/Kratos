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

#include "containers/model.h"
#include "custom_processes/apply_seepage_boundary_process.h"
#include "includes/kratos_parameters.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplySeepageBoundaryProcess_InfoReturnsClassName, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(ApplySeepageBoundaryProcess{}.Info(), "ApplySeepageBoundaryProcess"s);
}

KRATOS_TEST_CASE_IN_SUITE(ApplySeepageBoundaryProcess_ConstructorThrowsWhenNoModelPartNameIsProvided,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto model = Model{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplySeepageBoundaryProcess{model, Parameters{}}),
        "ApplySeepageBoundaryProcess: \"model_part_name\" is required in the process settings");
}

} // namespace Kratos::Testing