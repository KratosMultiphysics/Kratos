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
#include "containers/variables_list.h"
#include "custom_conditions/general_U_Pw_diff_order_condition.h"
#include "custom_utilities/registration_utilities.hpp"
#include "geo_aliases.h"
#include "geo_mechanics_application.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "includes/stream_serializer.h"
#include "test_setup_utilities/element_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <string>

using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(LineLoad2DDiffOrderCondition_CanBeSavedAndLoaded, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto registration = ScopedSerializerRegistration(
        std::make_pair("LineLoad2DDiffOrderCondition"s, LineLoad2DDiffOrderCondition{}),
        std::make_pair("Line2D3"s, Line2D3<Node>{PointerVector<Node>(3)}),
        std::make_pair("Line2D2"s, Line2D2<Node>{PointerVector<Node>(2)}));

    const auto solution_step_variables =
        Geo::ConstVariableDataRefs{std::cref(WATER_PRESSURE), std::cref(DISPLACEMENT)};
    const auto degrees_of_freedom =
        Geo::ConstVariableRefs{std::cref(WATER_PRESSURE), std::cref(DISPLACEMENT_X),
                               std::cref(DISPLACEMENT_Y), std::cref(DISPLACEMENT_Z)};
    auto p_condition = ElementSetupUtilities::Create2D3NLineCondition();
    ElementSetupUtilities::AddVariablesToEntity(p_condition, solution_step_variables, degrees_of_freedom);

    const auto dummy_process_info = ProcessInfo{};
    p_condition->Initialize(dummy_process_info);

    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_condition);
    auto p_loaded_condition = Condition::Pointer{nullptr};
    serializer.load("test_tag"s, p_loaded_condition);

    // Assert
    ASSERT_NE(p_loaded_condition, nullptr);

    auto dofs = Condition::DofsVectorType{};
    p_loaded_condition->GetDofList(dofs, dummy_process_info);
    const auto expected_dof_variables = Geo::ConstVariableRefs{
        std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y), std::cref(DISPLACEMENT_X),
        std::cref(DISPLACEMENT_Y), std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y),
        std::cref(WATER_PRESSURE), std::cref(WATER_PRESSURE)};
    ASSERT_EQ(dofs.size(), expected_dof_variables.size());
    for (auto i = std::size_t{0}; i < dofs.size(); ++i) {
        KRATOS_EXPECT_EQ(dofs[i]->GetVariable(), expected_dof_variables[i].get());
    }
}

} // namespace Kratos::Testing