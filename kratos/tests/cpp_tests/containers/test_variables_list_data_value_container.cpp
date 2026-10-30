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
#include <vector>

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

// Nodes share one VariablesList. An in-place load of an archive with other historical variables must throw without
// growing the shared list, so every container still destructs its block with the layout it was allocated with. The
// rejected container keeps its block and values.
KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerLoadDifferentListIntoSharedList, KratosCoreFastSuite) {
    auto p_saved_list = Kratos::make_intrusive<VariablesList>();
    p_saved_list->Add(DISTANCE);
    p_saved_list->Add(ELEMENTAL_DISTANCES);
    p_saved_list->Add(LOCAL_AXES_MATRIX);
    VariablesListDataValueContainer saved_container(p_saved_list, 2);
    saved_container.SetValue(LOCAL_AXES_MATRIX, Matrix(2, 2, 1.0), 0);

    auto p_variables_list = Kratos::make_intrusive<VariablesList>();
    p_variables_list->Add(DISTANCE);
    p_variables_list->Add(ELEMENTAL_DISTANCES);
    std::vector<VariablesListDataValueContainer> containers;
    containers.reserve(2);
    for (std::size_t j = 0; j < 2; ++j) {
        containers.emplace_back(p_variables_list, 2);
        for (std::size_t step = 0; step < 2; ++step) {
            containers[j].SetValue(DISTANCE, -1.0, step);
            containers[j].SetValue(ELEMENTAL_DISTANCES, Vector(5, -1.0), step);
        }
    }

    StreamSerializer serializer;
    serializer.save("Container", saved_container);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(serializer.load("Container", containers[0]), "The archive does not match this VariablesList");
    KRATOS_EXPECT_EQ(p_variables_list->size(), 2);
    for (std::size_t j = 0; j < 2; ++j) {
        for (std::size_t step = 0; step < 2; ++step) {
            KRATOS_EXPECT_EQ(containers[j].GetValue(DISTANCE, step), -1.0);
            KRATOS_EXPECT_VECTOR_EQ(containers[j].GetValue(ELEMENTAL_DISTANCES, step), Vector(5, -1.0));
        }
    }
    // the containers destruct their blocks (Vector values) with the unchanged list when they go out of scope
}

// Nodes of a model part without historical variables share an empty list and hold empty blocks. An in-place load of an
// archive with historical variables must throw without adding them to the list, since no block can hold them.
KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerLoadIntoSharedEmptyList, KratosCoreFastSuite) {
    auto p_saved_list = Kratos::make_intrusive<VariablesList>();
    p_saved_list->Add(DISTANCE);
    p_saved_list->Add(LOCAL_AXES_MATRIX);
    VariablesListDataValueContainer saved_container(p_saved_list, 2);
    saved_container.SetValue(LOCAL_AXES_MATRIX, Matrix(3, 3, 1.0), 0);
    StreamSerializer serializer;
    serializer.save("Container", saved_container);

    auto p_variables_list = Kratos::make_intrusive<VariablesList>();
    std::vector<VariablesListDataValueContainer> containers;
    containers.reserve(2);
    containers.emplace_back(p_variables_list, 2);
    containers.emplace_back(p_variables_list, 2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(serializer.load("Container", containers[0]), "The archive does not match this VariablesList");
    KRATOS_EXPECT_TRUE(p_variables_list->IsEmpty());
    KRATOS_EXPECT_EQ(p_variables_list->DataSize(), 0);
    KRATOS_EXPECT_FALSE(containers[1].Has(LOCAL_AXES_MATRIX));

    // the same holds for a container that is the only owner of its list
    VariablesListDataValueContainer container(Kratos::make_intrusive<VariablesList>(), 2);
    StreamSerializer other_serializer;
    other_serializer.save("Container", saved_container);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(other_serializer.load("Container", container), "The archive does not match this VariablesList");
    KRATOS_EXPECT_TRUE(container.GetVariablesList().IsEmpty());
    // destructing the containers with a grown list would destruct a Matrix in each of their empty blocks
}

// A saved list is loaded once per archive, so a container that shares it in the archive gets the list it was loaded into.
// An in-place load accepts that list only if it has the container's layout.
KRATOS_TEST_CASE_IN_SUITE(VariablesListDataValueContainerLoadListLoadedIntoAnotherList, KratosCoreFastSuite) {
    auto p_saved_list = Kratos::make_intrusive<VariablesList>();
    p_saved_list->Add(DISTANCE);
    VariablesListDataValueContainer saved_container(p_saved_list);
    saved_container.SetValue(DISTANCE, 2.0);
    const auto load_twice = [&saved_container](VariablesListDataValueContainer& rFirst, VariablesListDataValueContainer& rSecond) {
        StreamSerializer serializer;
        serializer.save("First", saved_container);
        serializer.save("Second", saved_container);
        serializer.load("First", rFirst);
        serializer.load("Second", rSecond);
    };

    auto p_list = Kratos::make_intrusive<VariablesList>();
    p_list->Add(DISTANCE);
    VariablesListDataValueContainer first(p_list);
    auto p_same_layout_list = Kratos::make_intrusive<VariablesList>();
    p_same_layout_list->Add(DISTANCE);
    VariablesListDataValueContainer second(p_same_layout_list);
    load_twice(first, second);
    KRATOS_EXPECT_EQ(second.GetValue(DISTANCE), 2.0);

    auto p_other_list = Kratos::make_intrusive<VariablesList>();
    p_other_list->Add(DISTANCE);
    p_other_list->Add(TEMPERATURE);
    VariablesListDataValueContainer other(p_other_list);
    other.SetValue(TEMPERATURE, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(load_twice(first, other), "An in-place load cannot change the historical variables");
    KRATOS_EXPECT_EQ(other.pGetVariablesList(), p_other_list);
    KRATOS_EXPECT_EQ(other.GetValue(TEMPERATURE), -1.0);
}

}
} // namespace Kratos.
