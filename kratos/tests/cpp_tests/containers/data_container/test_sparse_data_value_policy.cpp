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
#include "containers/data_container/data_container.h"
#include "containers/data_container/sparse_data_value_policy.h"
#include "includes/expect.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {
namespace {

// Number of entities used in the sparse scenarios (mirrors the 5x5x5 = 125 nodes of one
// voxel block in the original prototype tests)
constexpr std::size_t SparseTestNumberOfEntities = 125;

// Mark every even entity as active with a running sparse index and every odd one as
// inactive (-1): 63 active entries for 125 entities
void SetEvenIndexPattern(DataContainer& rContainer, const DataAccessor<int>& rIndexAccessor)
{
    int index = 0;
    auto index_span = rContainer.GetDataSpan(rIndexAccessor);
    for (std::size_t i = 0; i < index_span.size(); ++i) {
        if ((i % 2) == 0) {
            index_span[i] = index++;
        } else {
            index_span[i] = -1;
        }
    }
}

} // anonymous namespace

KRATOS_TEST_CASE_IN_SUITE(SparseDataValuePolicyIsSparse, KratosCoreFastSuite)
{
    DataValuePolicy<double> test_policy;
    KRATOS_EXPECT_FALSE(test_policy.IsSparse());

    DataContainer container(SparseTestNumberOfEntities);
    auto index_accessor = container.Add(STEP, DataValuePolicy<int>());
    SparseDataValuePolicy<double> test_sparse_policy(index_accessor);

    KRATOS_EXPECT_TRUE(test_sparse_policy.IsSparse());

    // A chunk created with a sparse policy starts empty
    auto value_accessor = container.Add(PRESSURE, test_sparse_policy);
    KRATOS_EXPECT_EQ(container.GetDataSpan(value_accessor).size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(SparseDataValuePolicyUpdateSparseStorage, KratosCoreFastSuite)
{
    DataContainer container(SparseTestNumberOfEntities);
    auto index_accessor = container.Add(STEP, DataValuePolicy<int>());
    auto value_accessor = container.Add(PRESSURE, SparseDataValuePolicy<double>(index_accessor));

    auto span_before_size_update = container.GetDataSpan(value_accessor);
    KRATOS_EXPECT_EQ(span_before_size_update.size(), 0);

    SetEvenIndexPattern(container, index_accessor);

    container.UpdateSparseStorage(index_accessor);

    auto span_after_size_update = container.GetDataSpan(value_accessor);
    KRATOS_EXPECT_EQ(span_after_size_update.size(), 125/2 + 1);
}

KRATOS_TEST_CASE_IN_SUITE(SparseDataValuePolicyUpdateMultiple, KratosCoreFastSuite)
{
    DataContainer container(SparseTestNumberOfEntities);

    // Two independent sparse index variables: value1 keyed on the first one,
    // value2 and value3 keyed on the second one
    auto index1_accessor = container.Add(STEP, DataValuePolicy<int>());
    auto index2_accessor = container.Add(DOMAIN_SIZE, DataValuePolicy<int>());
    auto value1_accessor = container.Add(PRESSURE, SparseDataValuePolicy<double>(index1_accessor));
    auto value2_accessor = container.Add(TEMPERATURE, SparseDataValuePolicy<double>(index2_accessor));
    auto value3_accessor = container.Add(TEMPERATURE_OLD_IT, SparseDataValuePolicy<double>(index2_accessor));

    KRATOS_EXPECT_EQ(container.GetDataSpan(value1_accessor).size(), 0);
    KRATOS_EXPECT_EQ(container.GetDataSpan(value2_accessor).size(), 0);
    KRATOS_EXPECT_EQ(container.GetDataSpan(value3_accessor).size(), 0);

    // Update the first index: only value1 grows
    SetEvenIndexPattern(container, index1_accessor);
    container.UpdateSparseStorage(index1_accessor);

    KRATOS_EXPECT_EQ(container.GetDataSpan(value1_accessor).size(), 125/2 + 1);
    KRATOS_EXPECT_EQ(container.GetDataSpan(value2_accessor).size(), 0);
    KRATOS_EXPECT_EQ(container.GetDataSpan(value3_accessor).size(), 0);

    // Second index: the first 32 entities are active, the rest inactive
    {
        auto index2_span = container.GetDataSpan(index2_accessor);
        for (std::size_t i = 0; i < 32; ++i) {
            index2_span[i] = static_cast<int>(i);
        }
        for (std::size_t i = 32; i < index2_span.size(); ++i) {
            index2_span[i] = -1;
        }
    }

    container.UpdateSparseStorage(index2_accessor);

    // value1 remains untouched, value2 and value3 grow to 32
    KRATOS_EXPECT_EQ(container.GetDataSpan(value1_accessor).size(), 125/2 + 1);
    KRATOS_EXPECT_EQ(container.GetDataSpan(value2_accessor).size(), 32);
    KRATOS_EXPECT_EQ(container.GetDataSpan(value3_accessor).size(), 32);
}

KRATOS_TEST_CASE_IN_SUITE(SparseDataValuePolicyResize, KratosCoreFastSuite)
{
    DataContainer container(SparseTestNumberOfEntities);
    auto index1_accessor = container.Add(STEP, DataValuePolicy<int>());
    auto value1_accessor = container.Add(PRESSURE, SparseDataValuePolicy<double>(index1_accessor));
    auto value2_accessor = container.Add(TEMPERATURE, SparseDataValuePolicy<double>(index1_accessor));

    SetEvenIndexPattern(container, index1_accessor);
    container.UpdateSparseStorage(index1_accessor);

    auto value1_span_after_update = container.GetDataSpan(value1_accessor);
    KRATOS_EXPECT_EQ(value1_span_after_update.size(), 125/2 + 1);

    auto value2_span_after_update = container.GetDataSpan(value2_accessor);
    KRATOS_EXPECT_EQ(value2_span_after_update.size(), 125/2 + 1);

    std::fill(value1_span_after_update.begin(), value1_span_after_update.end(), 1.0);

    std::vector<int> additional_entities;
    for (int i = 0; i < 125; i += 5) {
        additional_entities.push_back(i);
    }

    container.AddToSparseStorage(index1_accessor, additional_entities);

    // we added 12 new terms (indices 5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115)
    // to both containers indexed by index1
    auto value1_span_after_addition = container.GetDataSpan(value1_accessor);
    KRATOS_EXPECT_EQ(value1_span_after_addition.size(), 125/2 + 13);

    auto value2_span_after_addition = container.GetDataSpan(value2_accessor);
    KRATOS_EXPECT_EQ(value2_span_after_addition.size(), 125/2 + 13);

    auto index_span = container.GetDataSpan(index1_accessor);
    for (std::size_t i = 0; i < index_span.size(); ++i) {
        if (i % 2 == 0) {
            // original sparse entries kept their value
            KRATOS_EXPECT_EQ(value1_span_after_addition[index_span[i]], 1.0);
        } else if (i % 5 == 0) {
            // added entries are set to the default
            KRATOS_EXPECT_EQ(value1_span_after_addition[index_span[i]], 0.0);
        } else {
            // not in sparse container
            KRATOS_EXPECT_EQ(index_span[i], -1);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SparseDataValuePolicyResizeFromZero, KratosCoreFastSuite)
{
    DataContainer container(SparseTestNumberOfEntities);

    // The index variable uses -1 as its zero, so the whole index span starts inactive
    auto index1_accessor = container.Add(STEP, DataValuePolicy<int>(-1));
    auto value1_accessor = container.Add(PRESSURE, SparseDataValuePolicy<double>(index1_accessor));
    auto value2_accessor = container.Add(TEMPERATURE, SparseDataValuePolicy<double>(index1_accessor));

    KRATOS_EXPECT_EQ(container.GetDataSpan(value1_accessor).size(), 0);

    std::vector<int> additional_entities;
    for (int i = 0; i < 125; i += 5) {
        additional_entities.push_back(i);
    }

    container.AddToSparseStorage(index1_accessor, additional_entities);

    KRATOS_EXPECT_EQ(container.GetDataSpan(value1_accessor).size(), 25);
    KRATOS_EXPECT_EQ(container.GetDataSpan(value2_accessor).size(), 25);
}

} // namespace Testing
} // namespace Kratos
