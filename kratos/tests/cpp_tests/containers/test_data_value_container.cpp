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
#include "containers/data_value_container.h"
#include "includes/variables.h"
#include "testing/testing.h"
#include "includes/expect.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerHas, KratosCoreFastSuite) {
    DataValueContainer container;
    double area = 0.0;
    container.SetValue(NODAL_AREA, area);
    const bool check_true = container.Has(NODAL_AREA);
    const bool check_false = container.Has(DISPLACEMENT);

    KRATOS_EXPECT_EQ(check_true, true);
    KRATOS_EXPECT_EQ(check_false, false);
}

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerHasComponent, KratosCoreFastSuite) {
    DataValueContainer container;
    array_1d<double, 3> displacement;
    displacement[0] = 0.00;
    displacement[1] = 0.10;
    displacement[2] = 0.20;

    KRATOS_EXPECT_FALSE(container.Has(DISPLACEMENT_X));
    KRATOS_EXPECT_FALSE(container.Has(DISPLACEMENT_Y));
    KRATOS_EXPECT_FALSE(container.Has(DISPLACEMENT_Z));

    container.SetValue(DISPLACEMENT, displacement);

    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_X));
    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_Y));
    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_Z));
}

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerSetComponentCheck, KratosCoreFastSuite) {
    DataValueContainer container;

    KRATOS_EXPECT_FALSE(container.Has(DISPLACEMENT_X));
    KRATOS_EXPECT_FALSE(container.Has(DISPLACEMENT_Y));
    KRATOS_EXPECT_FALSE(container.Has(DISPLACEMENT_Z));

    container.SetValue(DISPLACEMENT_X, 10.0);

    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_X));
    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_Y));
    KRATOS_EXPECT_TRUE(container.Has(DISPLACEMENT_Z));

    KRATOS_EXPECT_EQ(container.GetValue(DISPLACEMENT_X), 10.0);
    KRATOS_EXPECT_EQ(container.GetValue(DISPLACEMENT_Y), 0.0);
    KRATOS_EXPECT_EQ(container.GetValue(DISPLACEMENT_Z), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(DataValueContainer, KratosCoreFastSuite) {
    DataValueContainer container;
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

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerComponent, KratosCoreFastSuite) {
    DataValueContainer container;
    array_1d<double, 3> original_velocity;
    original_velocity[0] = 0.00;
    original_velocity[1] = 0.10;
    original_velocity[2] = 0.20;

    container.SetValue(VELOCITY, original_velocity);

    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_X), original_velocity[0]);
    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_Y), original_velocity[1]);
    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_Z), original_velocity[2]);

    container.SetValue(DISPLACEMENT_Y, 1.23);

    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_X), original_velocity[0]);
    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_Y), original_velocity[1]);
    KRATOS_EXPECT_EQ(container.GetValue(VELOCITY_Z), original_velocity[2]);

    KRATOS_EXPECT_EQ(container.GetValue(DISPLACEMENT_X), 0.00);
    KRATOS_EXPECT_EQ(container.GetValue(DISPLACEMENT_Y), 1.23);
    KRATOS_EXPECT_EQ(container.GetValue(DISPLACEMENT_Z), 0.00);

}

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerMerge, KratosCoreFastSuite) {
    DataValueContainer container_origin;
    DataValueContainer container_target;

    const double density = 1000.0;
    const double viscosity_1 = 1e-3;
    const double viscosity_2 = 2e-3;

    container_origin.SetValue(DENSITY, density);
    container_origin.SetValue(VISCOSITY, viscosity_1);
    container_target.SetValue(VISCOSITY, viscosity_2);
    Flags options;
    options.Set(DataValueContainer::OVERWRITE_OLD_VALUES, false);
    container_target.Merge(container_origin, options);

    KRATOS_EXPECT_EQ(container_target.GetValue(DENSITY), density);
    KRATOS_EXPECT_EQ(container_target.GetValue(VISCOSITY), viscosity_2);
}

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerMergeOverride, KratosCoreFastSuite) {
    DataValueContainer container_origin;
    DataValueContainer container_target;

    const double density = 1000.0;
    const double viscosity_1 = 1e-3;
    const double viscosity_2 = 2e-3;

    container_origin.SetValue(DENSITY, density);
    container_origin.SetValue(VISCOSITY, viscosity_1);
    container_target.SetValue(VISCOSITY, viscosity_2);
    Flags options;
    options.Set(DataValueContainer::OVERWRITE_OLD_VALUES, true);
    container_target.Merge(container_origin, options);
    KRATOS_EXPECT_EQ(container_target.GetValue(DENSITY), density);
    KRATOS_EXPECT_EQ(container_target.GetValue(VISCOSITY), viscosity_1);

}

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerGetOrCreateValue, KratosCoreFastSuite) {
    DataValueContainer container;

    const double density = 1000.0;
    const double viscosity_1 = 1e-3;
    const double viscosity_2 = 2e-3;
    const array_1d<double, 3> velocity{0.1, 0.2, 0.3};

    const Vector vec_1 = Vector(5, 0.1);

    container.SetValue(DENSITY, density);
    container.SetValue(VISCOSITY, viscosity_1);
    container.SetValue(VISCOSITY, viscosity_2);
    container.SetValue(VELOCITY, velocity);

    KRATOS_EXPECT_EQ(container.Emplace(DENSITY), density);
    KRATOS_EXPECT_EQ(container.Emplace(VISCOSITY), viscosity_2);
    KRATOS_EXPECT_EQ(container.Emplace(PRESSURE), 0.0);
    KRATOS_EXPECT_VECTOR_EQ(container.Emplace(VELOCITY), velocity);
    KRATOS_EXPECT_EQ(container.Emplace(VELOCITY_Y), velocity[1]);


    KRATOS_EXPECT_EQ(container.Emplace(DISTANCE), 0.0);
    const array_1d<double, 3> momentum{0.9, 1.1, 1.2}, zero{};
    KRATOS_EXPECT_VECTOR_EQ(container.Emplace(MOMENTUM, momentum), momentum);

    KRATOS_EXPECT_VECTOR_EQ(container.Emplace(ACCELERATION), zero);
    KRATOS_EXPECT_EQ(container.Emplace(ACCELERATION_X), 0.0);

    KRATOS_EXPECT_EQ(container.Emplace(DISPLACEMENT_Z, 0.7), 0.7);
    KRATOS_EXPECT_EQ(container.Emplace(DISPLACEMENT_X, 0.7), 0.0);
    KRATOS_EXPECT_EQ(container.Emplace(DISPLACEMENT_Y, 0.7), 0.0);

    KRATOS_EXPECT_FALSE(container.Has(CONSTRAINT_LABELS))
    KRATOS_EXPECT_VECTOR_EQ(container.Emplace(CONSTRAINT_LABELS), vec_1);
    KRATOS_EXPECT_TRUE(container.Has(CONSTRAINT_LABELS))
}
}
} // namespace Kratos.
