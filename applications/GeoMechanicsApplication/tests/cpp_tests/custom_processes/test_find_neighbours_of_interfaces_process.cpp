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
//

#include "custom_processes/find_neighbours_of_interfaces_process.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(FindNeighboursOfInterfacesProcess_IsAProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto process = FindNeighboursOfInterfacesProcess{};
    EXPECT_NE(dynamic_cast<const Process*>(&process), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(FindNeighboursOfInterfacesProcess_FindsNoNeighboursWhenAdjacentModelPartIsEmpty,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto  model = Model{};
    auto& r_model_part_with_interface_element =
        ModelSetupUtilities::CreateModelPartWithASingle3D6NInterfaceElement(model);
    model.CreateModelPart("EmptyModelPart");
    auto process = FindNeighboursOfInterfacesProcess{};

    // Act
    process.ExecuteInitialize();

    // Assert
    EXPECT_EQ(r_model_part_with_interface_element.Elements().front().GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

} // namespace Kratos::Testing
