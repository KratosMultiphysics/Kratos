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
//

#include "custom_utilities/neighbouring_entity_finder.h"
#include "test_setup_utilities/element_setup_utilities.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NeighbouringEntityFinder_ReturnsEmptyListWhenNoNeighbouringElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part_for_neighbouring_elements = model.CreateModelPart("empty");

    NeighbouringEntityFinder finder;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);
    finder.InitializeBoundaryMaps(r_model_part_for_entities_for_finding.Elements());

    auto generate_generic_boundaries = [](const auto& rGeometry) {
        return rGeometry.GenerateBoundariesEntities();
    };

    finder.FindConditionNeighboursBasedOnBoundaryType(
        generate_generic_boundaries, r_model_part_for_neighbouring_elements.Elements());

    EXPECT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

} // namespace Kratos::Testing