//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "future/containers/entity_data_container.h"
#include "includes/expect.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(FutureEntityDataContainerRegistration, KratosCoreFutureSuite)
{
    Future::EntityDataContainer entity_data(1, 4);

    KRATOS_EXPECT_EQ(entity_data.NumberOfEntities(), 0);
    KRATOS_EXPECT_EQ(entity_data.Capacity(), 4);
    KRATOS_EXPECT_FALSE(entity_data.HasEntity(11));

    // Registration assigns slots in registration order and is idempotent
    KRATOS_EXPECT_EQ(entity_data.RegisterEntityId(11), 0);
    KRATOS_EXPECT_EQ(entity_data.RegisterEntityId(7), 1);
    KRATOS_EXPECT_EQ(entity_data.RegisterEntityId(11), 0); // idempotent
    KRATOS_EXPECT_EQ(entity_data.NumberOfEntities(), 2);
    KRATOS_EXPECT_TRUE(entity_data.HasEntity(11));
    KRATOS_EXPECT_EQ(entity_data.Index(7), 1);

    // Unknown ids raise clear errors
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        entity_data.Index(99),
        "Entity id 99 is not registered in this EntityDataContainer.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        entity_data.UnregisterEntityId(99),
        "Entity id 99 is not registered in this EntityDataContainer.");
}

KRATOS_TEST_CASE_IN_SUITE(FutureEntityDataContainerGrowthPreservesValues, KratosCoreFutureSuite)
{
    Future::EntityDataContainer entity_data(1, 2); // deliberately tiny chunk size
    auto accessor = entity_data.AddVariable(PRESSURE);

    entity_data.RegisterEntityId(1);
    entity_data.RegisterEntityId(2);
    entity_data.SetValue(std::size_t(1), PRESSURE, 10.0);
    entity_data.SetValue(std::size_t(2), PRESSURE, 20.0);

    // Registering past the chunk size grows the dense chunks preserving the values
    entity_data.RegisterEntityId(3);
    entity_data.RegisterEntityId(4);
    entity_data.RegisterEntityId(5);

    KRATOS_EXPECT_GE(entity_data.Capacity(), 5);
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(1), PRESSURE), 10.0);
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(2), PRESSURE), 20.0);
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(5), PRESSURE), 0.0); // new slots zero-initialized

    // The span covers the whole capacity; entity slots index into it
    auto span = entity_data.GetDataSpan(accessor);
    KRATOS_EXPECT_EQ(span.size(), entity_data.Capacity());
    KRATOS_EXPECT_EQ(span[entity_data.Index(std::size_t(2))], 20.0);
}

KRATOS_TEST_CASE_IN_SUITE(FutureEntityDataContainerLateVariableSizing, KratosCoreFutureSuite)
{
    Future::EntityDataContainer entity_data(1, 2);

    // Grow first, add the variable afterwards: the new chunk must match the capacity
    for (std::size_t id = 1; id <= 6; ++id) {
        entity_data.RegisterEntityId(id);
    }
    auto accessor = entity_data.AddVariable(TEMPERATURE);

    KRATOS_EXPECT_EQ(entity_data.GetDataSpan(accessor).size(), entity_data.Capacity());
    entity_data.SetValue(std::size_t(6), TEMPERATURE, 6.5);
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(6), TEMPERATURE), 6.5);
}

KRATOS_TEST_CASE_IN_SUITE(FutureEntityDataContainerUnregisterLeavesHole, KratosCoreFutureSuite)
{
    Future::EntityDataContainer entity_data(1, 4);
    entity_data.AddVariable(PRESSURE);

    entity_data.RegisterEntityId(1);
    entity_data.RegisterEntityId(2);
    entity_data.SetValue(std::size_t(2), PRESSURE, 2.5);
    const auto old_capacity = entity_data.Capacity();

    entity_data.UnregisterEntityId(2);
    KRATOS_EXPECT_FALSE(entity_data.HasEntity(2));
    KRATOS_EXPECT_EQ(entity_data.NumberOfEntities(), 1);
    KRATOS_EXPECT_EQ(entity_data.Capacity(), old_capacity); // the slot is a hole, not reclaimed

    // Re-registering assigns a NEW slot (zero-initialized when past the old data)
    const auto new_slot = entity_data.RegisterEntityId(2);
    KRATOS_EXPECT_EQ(new_slot, 2); // slot 1 is the hole; next fresh slot is 2
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(2), PRESSURE), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(FutureEntityDataContainerHistoricalSteps, KratosCoreFutureSuite)
{
    Future::EntityDataContainer entity_data(3, 4); // 3 buffer steps
    KRATOS_EXPECT_EQ(entity_data.GetBufferSize(), 3);

    entity_data.RegisterEntityId(1);
    entity_data.RegisterEntityId(2);
    entity_data.AddHistoricalVariable(PRESSURE);

    entity_data.SetValue(std::size_t(1), PRESSURE, 1.0);
    entity_data.SetValue(std::size_t(2), PRESSURE, 2.0);

    // Cloning the time step copies the current values onto the new current step
    entity_data.CloneStepData(StepCategory::TimeStep);
    entity_data.SetValue(std::size_t(1), PRESSURE, 10.0);

    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(1), PRESSURE), 10.0);
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(1), PRESSURE, 1), 1.0);  // previous step
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(2), PRESSURE), 2.0);     // cloned, unmodified
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(2), PRESSURE, 1), 2.0);
}

KRATOS_TEST_CASE_IN_SUITE(FutureEntityDataContainerErrorPaths, KratosCoreFutureSuite)
{
    Future::EntityDataContainer entity_data(1, 4);
    entity_data.RegisterEntityId(1);
    entity_data.AddVariable(PRESSURE);

    // Un-added variable
    KRATOS_EXPECT_FALSE(entity_data.Has(TEMPERATURE));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        entity_data.GetValue(std::size_t(1), TEMPERATURE),
        "Variable TEMPERATURE has not been added to this EntityDataContainer.");

    // Unknown entity
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        entity_data.GetValue(std::size_t(9), PRESSURE),
        "Entity id 9 is not registered in this EntityDataContainer.");

    // Step request on a non-historical variable
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        entity_data.GetValue(std::size_t(1), PRESSURE, 1),
        "Variable PRESSURE is not historical");

    KRATOS_EXPECT_TRUE(entity_data.Has(PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(FutureEntityDataContainerAccessorCoherence, KratosCoreFutureSuite)
{
    Future::EntityDataContainer entity_data(1, 4);
    entity_data.RegisterEntityId(3);
    auto add_accessor = entity_data.AddVariable(PRESSURE);

    entity_data.SetValue(std::size_t(3), PRESSURE, 7.5);

    // The accessor returned by Add, the one rebuilt by GetAccessor and GetValue all see the same storage
    auto rebuilt_accessor = entity_data.GetAccessor(PRESSURE);
    KRATOS_EXPECT_EQ(entity_data.GetDataSpan(add_accessor)[entity_data.Index(std::size_t(3))], 7.5);
    KRATOS_EXPECT_EQ(entity_data.GetDataSpan(rebuilt_accessor)[entity_data.Index(std::size_t(3))], 7.5);
    KRATOS_EXPECT_EQ(entity_data.GetValue(std::size_t(3), PRESSURE), 7.5);
}

KRATOS_TEST_CASE_IN_SUITE(FutureEntityDataContainerCopyIsShallow, KratosCoreFutureSuite)
{
    Future::EntityDataContainer entity_data(1, 4);
    entity_data.RegisterEntityId(1);
    entity_data.AddVariable(PRESSURE);

    Future::EntityDataContainer copy(entity_data);
    entity_data.SetValue(std::size_t(1), PRESSURE, 3.5);

    // Copies share the chunk storage (DataContainer shallow-copy semantics)
    KRATOS_EXPECT_EQ(copy.GetValue(std::size_t(1), PRESSURE), 3.5);
}

} // namespace Kratos::Testing
