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
#include "containers/data_container/data_chunk.h"
#include "includes/expect.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DataChunkConstructionAndZeroInitialization, KratosCoreFastSuite)
{
    // Default zero
    DataChunk<double> chunk(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy(), 8);

    KRATOS_EXPECT_EQ(chunk.NumberOfEntitiesPerStep(), 8);
    KRATOS_EXPECT_TRUE(chunk.GetVariableData() == PRESSURE);
    for (std::size_t i = 0; i < 8; ++i) {
        KRATOS_EXPECT_EQ(chunk.GetData()[i], 0.0);
    }

    // Custom policy zero initializes every entry
    DataChunk<double> custom_zero_chunk(PRESSURE, DataValuePolicy<double>(3.5), NonHistoricalDataPolicy(), 8);
    for (std::size_t i = 0; i < 8; ++i) {
        KRATOS_EXPECT_EQ(custom_zero_chunk.GetData()[i], 3.5);
    }

    // A historical chunk allocates and zero-fills all step slots
    DataChunk<double> historical_chunk(PRESSURE, DataValuePolicy<double>(1.25), HistoricalDataPolicy(StepCategory::TimeStep, 3), 4);
    for (std::size_t i = 0; i < 4 * 3; ++i) {
        KRATOS_EXPECT_EQ(historical_chunk.GetData()[i], 1.25);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataChunkZeroEntities, KratosCoreFastSuite)
{
    // Sparse chunks legitimately start with zero entities: construction, CreateNew(0) and
    // growing from empty must all work (also in debug builds)
    DataChunk<double> empty_chunk(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy(), 0);
    KRATOS_EXPECT_EQ(empty_chunk.NumberOfEntitiesPerStep(), 0);
    KRATOS_EXPECT_EQ(empty_chunk.GetData(), nullptr);

    auto p_new_empty = empty_chunk.CreateNew(0);
    KRATOS_EXPECT_EQ(p_new_empty->NumberOfEntitiesPerStep(), 0);

    // Growing from empty zero-fills the new entries
    empty_chunk.ResizeData(5);
    KRATOS_EXPECT_EQ(empty_chunk.NumberOfEntitiesPerStep(), 5);
    for (std::size_t i = 0; i < 5; ++i) {
        KRATOS_EXPECT_EQ(empty_chunk.GetData()[i], 0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataChunkMultiStepLayout, KratosCoreFastSuite)
{
    DataChunk<double> chunk(PRESSURE, DataValuePolicy<double>(), HistoricalDataPolicy(StepCategory::TimeStep, 3), 4);

    // Write distinct values per step slot through the per-step pointers
    for (std::size_t step = 0; step < 3; ++step) {
        double* p_step_data = chunk.GetDataPointer(step);
        for (std::size_t i = 0; i < 4; ++i) {
            p_step_data[i] = 10.0 * step + i;
        }
    }

    // The storage is laid out step-by-step
    for (std::size_t step = 0; step < 3; ++step) {
        for (std::size_t i = 0; i < 4; ++i) {
            KRATOS_EXPECT_EQ(chunk.GetData()[step * 4 + i], 10.0 * step + i);
            KRATOS_EXPECT_EQ(chunk.GetDataIterator(step)[i], 10.0 * step + i);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataChunkCloneStepData, KratosCoreFastSuite)
{
    DataChunk<double> chunk(PRESSURE, DataValuePolicy<double>(), HistoricalDataPolicy(StepCategory::TimeStep, 3), 4);

    // Fill the current step (slot 0)
    for (std::size_t i = 0; i < 4; ++i) {
        chunk.GetDataPointer(0)[i] = 1.0;
    }

    // Cloning a non-matching category is a no-op
    chunk.CloneStepData(StepCategory::IterationStep);
    KRATOS_EXPECT_EQ(chunk.GetHistoryPolicy().GetLatestStepIndex(), 0);
    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_EXPECT_EQ(chunk.GetDataPointer(1)[i], 0.0);
    }

    // Cloning the matching category advances the ring and copies the current data
    chunk.CloneStepData(StepCategory::TimeStep);
    KRATOS_EXPECT_EQ(chunk.GetHistoryPolicy().GetLatestStepIndex(), 1);
    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_EXPECT_EQ(chunk.GetDataPointer(1)[i], 1.0); // new current step is a copy
        KRATOS_EXPECT_EQ(chunk.GetDataPointer(0)[i], 1.0); // previous step untouched
    }

    // Two more clones wrap around: slot 0 becomes the current step again and is
    // overwritten with the data of slot 2
    chunk.GetDataPointer(1)[0] = 7.0;
    chunk.CloneStepData(StepCategory::TimeStep); // latest: 2, copy of slot 1
    chunk.CloneStepData(StepCategory::TimeStep); // latest: 0 (wrap), copy of slot 2
    KRATOS_EXPECT_EQ(chunk.GetHistoryPolicy().GetLatestStepIndex(), 0);
    KRATOS_EXPECT_EQ(chunk.GetDataPointer(0)[0], 7.0);
}

KRATOS_TEST_CASE_IN_SUITE(DataChunkResizePreservesPerStepData, KratosCoreFastSuite)
{
    // Regression test: the resize copy has to work per step slot, not on the total size
    DataChunk<double> chunk(PRESSURE, DataValuePolicy<double>(), HistoricalDataPolicy(StepCategory::TimeStep, 3), 8);

    for (std::size_t step = 0; step < 3; ++step) {
        for (std::size_t i = 0; i < 8; ++i) {
            chunk.GetDataPointer(step)[i] = 100.0 * step + i;
        }
    }

    // Shrink: each step keeps its first 5 values
    chunk.ResizeData(5);
    KRATOS_EXPECT_EQ(chunk.NumberOfEntitiesPerStep(), 5);
    for (std::size_t step = 0; step < 3; ++step) {
        for (std::size_t i = 0; i < 5; ++i) {
            KRATOS_EXPECT_EQ(chunk.GetDataPointer(step)[i], 100.0 * step + i);
        }
    }

    // Grow: each step keeps its 5 values and zero-fills the new tail
    chunk.ResizeData(12);
    KRATOS_EXPECT_EQ(chunk.NumberOfEntitiesPerStep(), 12);
    for (std::size_t step = 0; step < 3; ++step) {
        for (std::size_t i = 0; i < 5; ++i) {
            KRATOS_EXPECT_EQ(chunk.GetDataPointer(step)[i], 100.0 * step + i);
        }
        for (std::size_t i = 5; i < 12; ++i) {
            KRATOS_EXPECT_EQ(chunk.GetDataPointer(step)[i], 0.0);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataChunkResizeToZeroAndBack, KratosCoreFastSuite)
{
    // Regression test: resizing to zero releases the storage without leaving a dangling
    // pointer (the destructor must not double-free), and the chunk can grow again
    DataChunk<double> chunk(PRESSURE, DataValuePolicy<double>(2.0), NonHistoricalDataPolicy(), 6);

    chunk.ResizeData(0);
    KRATOS_EXPECT_EQ(chunk.NumberOfEntitiesPerStep(), 0);
    KRATOS_EXPECT_EQ(chunk.GetData(), nullptr);

    chunk.ResizeData(4);
    KRATOS_EXPECT_EQ(chunk.NumberOfEntitiesPerStep(), 4);
    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_EXPECT_EQ(chunk.GetData()[i], 2.0); // zero-filled with the policy zero
    }

    // Destruction after ResizeData(0) alone must also be safe
    {
        DataChunk<double> to_destroy(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy(), 3);
        to_destroy.ResizeData(0);
    } // no double delete here
}

KRATOS_TEST_CASE_IN_SUITE(DataChunkCreateNewPreservesPolicies, KratosCoreFastSuite)
{
    DataValuePolicy<double> value_policy(4.5);
    HistoricalDataPolicy history_policy(StepCategory::TimeStep, 2);
    DataChunk<double> chunk(PRESSURE, value_policy, history_policy, 3);

    auto p_new = chunk.CreateNew(7);

    KRATOS_EXPECT_EQ(p_new->NumberOfEntitiesPerStep(), 7);
    KRATOS_EXPECT_TRUE(p_new->GetVariableData() == PRESSURE);
    KRATOS_EXPECT_TRUE(p_new->GetValuePolicy().IsSame(value_policy));
    KRATOS_EXPECT_TRUE(p_new->GetHistoryPolicy().IsSame(history_policy));

    // The new chunk is zero-initialized with the (cloned) policy zero
    auto p_typed = std::dynamic_pointer_cast<DataChunk<double>>(p_new);
    KRATOS_EXPECT_TRUE(p_typed != nullptr);
    for (std::size_t i = 0; i < 7 * 2; ++i) {
        KRATOS_EXPECT_EQ(p_typed->GetData()[i], 4.5);
    }
}

} // namespace Testing
} // namespace Kratos
