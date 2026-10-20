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
#include "includes/constitutive_law.h"
#include "includes/stream_serializer.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerHas, KratosCoreFastSuite) {
    VariablesList::Pointer p_variables_list=Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(NODAL_AREA);
    VariablesListDataValueContainer container(p_variables_list);

    KRATOS_EXPECT_TRUE(container.Has(NODAL_AREA));
    KRATOS_EXPECT_FALSE(container.Has(DISPLACEMENT));
}

KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerHasComponent, KratosCoreFastSuite) {
    VariablesList::Pointer p_variables_list=Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(DISPLACEMENT);
    VariablesListDataValueContainer container(p_variables_list);
 
    KRATOS_EXPECT_FALSE(container.Has(VELOCITY_X));
    KRATOS_EXPECT_FALSE(container.Has(VELOCITY_Y));
    KRATOS_EXPECT_FALSE(container.Has(VELOCITY_Z));

    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_X));
    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_Y));
    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_Z));
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
        KRATOS_EXPECT_EQ(distances[i], original_distances[i]);
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

    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_X), original_velocity[0]);
    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_Y), original_velocity[1]);
    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_Z), original_velocity[2]);
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

    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_X), original_velocity[0]);
    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_Y), original_velocity[1]);
    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_Z), original_velocity[2]);
}

// A load into a populated container must destruct and free its old block; the law's use count shows
// whether a leaked block still holds copies of its pointer (one per buffer step).
KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerLoadIntoPopulatedContainer, KratosCoreFastSuite) {
    VariablesList::Pointer p_variables_list = Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(CONSTITUTIVE_LAW);
    p_variables_list->Add(ELEMENTAL_DISTANCES);
    p_variables_list->Add(DISPLACEMENT);
    VariablesListDataValueContainer container(p_variables_list, 2);

    auto p_law = Kratos::make_shared<ConstitutiveLaw>();
    const std::vector<array_1d<double, 3>> displacements{
        array_1d<double, 3>{0.1, 0.2, 0.3}, array_1d<double, 3>{0.4, 0.5, 0.6}};
    for (std::size_t step = 0; step < 2; ++step) {
        container.SetValue(CONSTITUTIVE_LAW, p_law, step);
        container.SetValue(ELEMENTAL_DISTANCES, Vector(4, step + 1.0), step);
        container.SetValue(DISPLACEMENT, displacements[step], step);
    }
    KRATOS_EXPECT_EQ(p_law.use_count(), 3);

    for (std::size_t i = 0; i < 2; ++i) {
        StreamSerializer serializer;
        serializer.save("Container", container);
        container.SetValue(DISPLACEMENT_X, -1.0);
        serializer.load("Container", container);

        KRATOS_EXPECT_EQ(container.QueueSize(), 2);
        for (std::size_t step = 0; step < 2; ++step) {
            KRATOS_EXPECT_NE(container.GetValue(CONSTITUTIVE_LAW, step), nullptr);
            KRATOS_EXPECT_NE(container.GetValue(CONSTITUTIVE_LAW, step), p_law);
            KRATOS_EXPECT_VECTOR_EQ(container.GetValue(ELEMENTAL_DISTANCES, step), Vector(4, step + 1.0));
            KRATOS_EXPECT_VECTOR_EQ(container.GetValue(DISPLACEMENT, step), displacements[step]);
        }
    }

    KRATOS_EXPECT_EQ(p_law.use_count(), 1);
}

}
} // namespace Kratos.
