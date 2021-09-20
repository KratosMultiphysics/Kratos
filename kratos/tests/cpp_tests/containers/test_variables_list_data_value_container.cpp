//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerHas, KratosCoreFastSuite) {
    VariablesList::Pointer p_variables_list=Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(NODAL_AREA);
    VariablesListDataValueContainer container(p_variables_list);

    KRATOS_CHECK(container.Has(NODAL_AREA));
    KRATOS_CHECK_IS_FALSE(container.Has(DISPLACEMENT));
}

KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerHasComponent, KratosCoreFastSuite) {
    VariablesList::Pointer p_variables_list=Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(DISPLACEMENT);
    VariablesListDataValueContainer container(p_variables_list);
 
    KRATOS_CHECK_IS_FALSE(container.Has(VELOCITY_X));
    KRATOS_CHECK_IS_FALSE(container.Has(VELOCITY_Y));
    KRATOS_CHECK_IS_FALSE(container.Has(VELOCITY_Z));

    KRATOS_CHECK(container.Has(DISPLACEMENT_X));
    KRATOS_CHECK(container.Has(DISPLACEMENT_Y));
    KRATOS_CHECK(container.Has(DISPLACEMENT_Z));
}

KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerGetValue, KratosCoreFastSuite) {
    VariablesList::Pointer p_variables_list=Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(ELEMENTAL_DISTANCES);
    VariablesListDataValueContainer container(p_variables_list);
    Vector original_distances(4);
    original_distances[0] = 0.00;
    original_distances[1] = 0.10;
    original_distances[2] = 0.20;
    original_distances[3] = 0.30;
    container.SetValue(ELEMENTAL_DISTANCES, original_distances);
    auto& distances = container.GetValue(ELEMENTAL_DISTANCES);

    for (std::size_t i = 0; i < distances.size(); i++)
        KRATOS_CHECK_EQUAL(distances[i], original_distances[i]);
}

KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerAddComponent, KratosCoreFastSuite) {
    VariablesList::Pointer p_variables_list=Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(VELOCITY_Y);
    VariablesListDataValueContainer container(p_variables_list);
    array_1d<double, 3> original_velocity;
    original_velocity[0] = 0.00;
    original_velocity[1] = 0.10;
    original_velocity[2] = 0.20;

    container.SetValue(VELOCITY, original_velocity);

    KRATOS_CHECK_EQUAL(container.GetValue(VELOCITY_X), original_velocity[0]);
    KRATOS_CHECK_EQUAL(container.GetValue(VELOCITY_Y), original_velocity[1]);
    KRATOS_CHECK_EQUAL(container.GetValue(VELOCITY_Z), original_velocity[2]);
}

KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerGetComponentValue, KratosCoreFastSuite) {
    VariablesList::Pointer p_variables_list=Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(VELOCITY);
    VariablesListDataValueContainer container(p_variables_list);
    array_1d<double, 3> original_velocity;
    original_velocity[0] = 0.00;
    original_velocity[1] = 0.10;
    original_velocity[2] = 0.20;

    container.SetValue(VELOCITY, original_velocity);

    KRATOS_CHECK_EQUAL(container.GetValue(VELOCITY_X), original_velocity[0]);
    KRATOS_CHECK_EQUAL(container.GetValue(VELOCITY_Y), original_velocity[1]);
    KRATOS_CHECK_EQUAL(container.GetValue(VELOCITY_Z), original_velocity[2]);
}

}
} // namespace Kratos.
