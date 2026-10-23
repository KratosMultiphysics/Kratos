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
#include <sstream>
#include <unordered_set>

// External includes

// Project includes
#include "containers/data_value_container.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"
#include "includes/stream_serializer.h"
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

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerEmplace, KratosCoreFastSuite) {
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
    const array_1d<double, 3> momentum{0.9, 1.1, 1.2};
    const array_1d<double, 3> zero{0.0, 0.0, 0.0};
    KRATOS_EXPECT_VECTOR_EQ(container.Emplace(MOMENTUM, momentum), momentum);

    KRATOS_EXPECT_FALSE(container.Has(ACCELERATION))
    KRATOS_EXPECT_VECTOR_EQ(container.Emplace(ACCELERATION), zero);
    KRATOS_EXPECT_TRUE(container.Has(ACCELERATION))
    KRATOS_EXPECT_EQ(container.Emplace(ACCELERATION_X), 0.0);

    KRATOS_EXPECT_EQ(container.Emplace(DISPLACEMENT_Z, 0.7), 0.7);
    KRATOS_EXPECT_EQ(container.Emplace(DISPLACEMENT_X, 0.7), 0.0);
    KRATOS_EXPECT_EQ(container.Emplace(DISPLACEMENT_Y, 0.7), 0.0);

    KRATOS_EXPECT_FALSE(container.Has(CONSTRAINT_LABELS))
    KRATOS_EXPECT_VECTOR_EQ(container.Emplace(CONSTRAINT_LABELS, vec_1), vec_1);
    KRATOS_EXPECT_TRUE(container.Has(CONSTRAINT_LABELS))

    container.Emplace(CONSTRAINT_LABELS)[2] = 0.2;
    Vector vec_2 = Vector(5, 0.1);
    vec_2[2] = 0.2;
    KRATOS_EXPECT_VECTOR_EQ(container.GetValue(CONSTRAINT_LABELS), vec_2);
}

// A non-DataOnly load into a populated container must delete the old values; the law's use count
// shows whether the container still holds a leaked copy of its pointer.
KRATOS_TEST_CASE_IN_SUITE(DataValueContainerLoadIntoPopulatedContainer, KratosCoreFastSuite) {
    auto p_law = Kratos::make_shared<ConstitutiveLaw>();
    const Vector distances(3, 1.5);

    DataValueContainer container;
    container.SetValue(CONSTITUTIVE_LAW, p_law);
    container.SetValue(DENSITY, 2.0);
    container.SetValue(ELEMENTAL_DISTANCES, distances);
    KRATOS_EXPECT_EQ(p_law.use_count(), 2);

    for (std::size_t i = 0; i < 2; ++i) {
        StreamSerializer serializer;
        serializer.save("Container", container);
        container.SetValue(DENSITY, -1.0);
        serializer.load("Container", container);

        KRATOS_EXPECT_EQ(container.Size(), 3);
        KRATOS_EXPECT_EQ(container.GetValue(DENSITY), 2.0);
        KRATOS_EXPECT_VECTOR_EQ(container.GetValue(ELEMENTAL_DISTANCES), distances);
        KRATOS_EXPECT_NE(container.GetValue(CONSTITUTIVE_LAW), nullptr);
        KRATOS_EXPECT_NE(container.GetValue(CONSTITUTIVE_LAW), p_law);
    }

    KRATOS_EXPECT_EQ(p_law.use_count(), 1);
}

namespace {

/// In-memory serializer in DataOnly mode (saves and loads only values, into existing objects).
class DataOnlyStreamSerializer : public Serializer
{
public:
    DataOnlyStreamSerializer()
        : Serializer(new std::stringstream(std::ios::binary | std::ios::in | std::ios::out), SERIALIZER_NO_TRACE, true)
    {}
};

} // namespace

// roigcarlo's case 1: DataOnly into a container with the same variables in the same order.
KRATOS_TEST_CASE_IN_SUITE(DataValueContainerDataOnlyLoadSameLayout, KratosCoreFastSuite) {
    DataValueContainer saved;
    saved.SetValue(DENSITY, 2.0);
    saved.SetValue(ELEMENTAL_DISTANCES, Vector(3, 1.5));

    DataValueContainer container;
    container.SetValue(DENSITY, -1.0);
    container.SetValue(ELEMENTAL_DISTANCES, Vector(3, -1.0));

    DataOnlyStreamSerializer serializer;
    serializer.save("Container", saved);
    serializer.load("Container", container);

    KRATOS_EXPECT_EQ(container.Size(), 2);
    KRATOS_EXPECT_EQ(container.GetValue(DENSITY), 2.0);
    KRATOS_EXPECT_VECTOR_EQ(container.GetValue(ELEMENTAL_DISTANCES), Vector(3, 1.5));
}

// roigcarlo's case 2: DataOnly into a container holding a value added after the save; restoring removes and
// frees it (the law's use count shows a leaked copy).
KRATOS_TEST_CASE_IN_SUITE(DataValueContainerDataOnlyLoadSmallerArchive, KratosCoreFastSuite) {
    auto p_law = Kratos::make_shared<ConstitutiveLaw>();

    DataValueContainer container;
    container.SetValue(DENSITY, 2.0);
    container.SetValue(ELEMENTAL_DISTANCES, Vector(3, 1.5));

    DataOnlyStreamSerializer serializer;
    serializer.save("Container", container);
    container.SetValue(DENSITY, -1.0);
    container.SetValue(CONSTITUTIVE_LAW, p_law);
    KRATOS_EXPECT_EQ(p_law.use_count(), 2);
    serializer.load("Container", container);

    KRATOS_EXPECT_EQ(container.Size(), 2);
    KRATOS_EXPECT_FALSE(container.Has(CONSTITUTIVE_LAW));
    KRATOS_EXPECT_EQ(container.GetValue(DENSITY), 2.0);
    KRATOS_EXPECT_VECTOR_EQ(container.GetValue(ELEMENTAL_DISTANCES), Vector(3, 1.5));
    KRATOS_EXPECT_EQ(p_law.use_count(), 1);
}

// DataOnly with the same variables in another order (and so other types per position): values are restored by name.
KRATOS_TEST_CASE_IN_SUITE(DataValueContainerDataOnlyLoadReordered, KratosCoreFastSuite) {
    DataValueContainer saved;
    saved.SetValue(DENSITY, 2.0);
    saved.SetValue(ELEMENTAL_DISTANCES, Vector(3, 1.5));

    DataValueContainer container;
    container.SetValue(ELEMENTAL_DISTANCES, Vector(3, -1.0));
    container.SetValue(DENSITY, -1.0);

    DataOnlyStreamSerializer serializer;
    serializer.save("Container", saved);
    serializer.load("Container", container);

    KRATOS_EXPECT_EQ(container.Size(), 2);
    KRATOS_EXPECT_EQ(container.GetValue(DENSITY), 2.0);
    KRATOS_EXPECT_VECTOR_EQ(container.GetValue(ELEMENTAL_DISTANCES), Vector(3, 1.5));
}

// DataOnly cannot create values, so a saved variable missing from the container is an error.
KRATOS_TEST_CASE_IN_SUITE(DataValueContainerDataOnlyLoadMissingVariable, KratosCoreFastSuite) {
    DataValueContainer saved;
    saved.SetValue(DENSITY, 2.0);
    saved.SetValue(ELEMENTAL_DISTANCES, Vector(3, 1.5));

    DataValueContainer container;
    container.SetValue(ELEMENTAL_DISTANCES, Vector(3, -1.0));

    DataOnlyStreamSerializer serializer;
    serializer.save("Container", saved);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(serializer.load("Container", container), "DataOnly loading restores values into an existing container, but the variable DENSITY is not in it.");
}

}
} // namespace Kratos.
