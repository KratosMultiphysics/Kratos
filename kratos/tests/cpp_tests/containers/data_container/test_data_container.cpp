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
#include "containers/array_1d.h"
#include "containers/data_container/data_container.h"
#include "containers/data_container/sparse_data_value_policy.h"
#include "includes/expect.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {
namespace {

// Custom user type stored in the container: needs default construction, operator== and
// operator<<; save/load are required by the Serializer support of the real Variable<T>
class CustomClass
{
public:
    int value = 0;
    bool operator==(const CustomClass& rOther) const { return value == rOther.value; }
    bool operator!=(const CustomClass& rOther) const { return !(*this == rOther); }
    CustomClass& operator=(const CustomClass& rOther) { value = rOther.value; return *this; }
    void save(Serializer& rSerializer) const { rSerializer.save("value", value); }
    void load(Serializer& rSerializer) { rSerializer.load("value", value); }
};

std::ostream& operator<<(std::ostream& rOStream, const CustomClass& rObject)
{
    rOStream << "CustomClass(value=" << rObject.value << ")";
    return rOStream;
}

// Generic round trip: adding a variable creates a zero-initialized span of the default
// chunk size; writing through the span is visible through a fresh accessor-based span
template<typename T, typename TVariableValueType, typename TPolicyType>
void DataContainerGetDataTest(
    const Variable<TVariableValueType>& rVariable,
    TPolicyType const& rValuePolicy,
    T const& rZero,
    T const& rOne)
{
    DataContainer container;

    auto accessor = container.Add(rVariable, rValuePolicy, NonHistoricalDataPolicy());

    // data container initializes values to zero
    auto data_span = container.GetDataSpan(rVariable, rValuePolicy);
    KRATOS_EXPECT_EQ(data_span.size(), 256); // default chunk size
    for (auto& r_value : data_span) {
        KRATOS_EXPECT_EQ(r_value, rZero);
    }

    for (auto& r_value : data_span) {
        r_value = rOne; // modify the data
    }

    auto data_span2 = container.GetDataSpan(accessor);
    KRATOS_EXPECT_EQ(data_span2.size(), 256); // default chunk size
    for (auto& r_value : data_span2) {
        KRATOS_EXPECT_EQ(r_value, rOne); // check that the data was modified
    }
}

// Helper for DataContainerVariableLifeTime: defines a variable on the spot, registers it
// and uses it to store data; the variable object dies when the function returns
void AddDataForLocalVariable(DataContainer& rContainer)
{
    Variable<double> local{"DATA_CONTAINER_TEST_LOCAL"};
    local.Register();
    auto accessor = rContainer.Add(local, DataValuePolicy<double>(), NonHistoricalDataPolicy());
    for (auto& r_value : rContainer.GetDataSpan(accessor)) {
        r_value = 5.7;
    }
}

} // anonymous namespace

KRATOS_TEST_CASE_IN_SUITE(DataContainerAdd, KratosCoreFastSuite)
{
    DataContainer container;

    KRATOS_EXPECT_EQ(container.Has(PRESSURE), false);
    KRATOS_EXPECT_EQ(container.Has(PRESSURE, DataValuePolicy<double>()), false);

    container.Add(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy());

    KRATOS_EXPECT_EQ(container.Has(PRESSURE), true);
    KRATOS_EXPECT_EQ(container.Has(PRESSURE, DataValuePolicy<double>()), true);
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerGetData, KratosCoreFastSuite)
{
    // Test with double data type
    DataContainerGetDataTest<double>(PRESSURE, DataValuePolicy<double>(), 0.0, 1.0);

    // Test with array_1d<double, 3> data type
    // Note: array_1d's default constructor does not zero-initialize its storage (unlike
    // std::array), so the policy zero must be passed explicitly here.
    const array_1d<double, 3> array_zero{{0.0, 0.0, 0.0}};
    DataContainerGetDataTest(VELOCITY, DataValuePolicy<array_1d<double, 3>>(array_zero),
                             array_zero, array_1d<double, 3>{{1.0, 1.0, 1.0}});

    // Test with a layered policy over a double variable
    const std::size_t number_of_layers = 10;
    DataContainerGetDataTest(TEMPERATURE, LayeredDataValuePolicy<double>(number_of_layers),
                             std::vector<double>(number_of_layers, 0.0), std::vector<double>(number_of_layers, 1.0));

    // Test with int data type
    DataContainerGetDataTest<int>(STEP, DataValuePolicy<int>(), 0, 1);

    // Test with std::string data type
    DataContainerGetDataTest<std::string>(IDENTIFIER, DataValuePolicy<std::string>(), "", "Test");

    // Test with a custom user type (registered on the spot)
    Variable<CustomClass> custom_variable("DATA_CONTAINER_TEST_CUSTOM");
    custom_variable.Register();
    DataContainerGetDataTest(custom_variable, DataValuePolicy<CustomClass>(), CustomClass(), CustomClass{1});
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerTimeStep, KratosCoreFastSuite)
{
    DataContainer container;
    DataValuePolicy<double> value_policy;
    HistoricalDataPolicy historical_policy(StepCategory::TimeStep, 3); // Storing 3 time steps including the current one

    auto double_accessor = container.Add(PRESSURE, value_policy, historical_policy);

    // data container initializes values to zero
    auto current_step_data = container.GetDataSpan(PRESSURE, value_policy, StepCategory::TimeStep);
    KRATOS_EXPECT_EQ(current_step_data.size(), 256); // default chunk size
    for (auto& r_value : current_step_data) {
        KRATOS_EXPECT_EQ(r_value, 0.0);
    }

    for (auto& r_value : current_step_data) {
        r_value = 1.0; // modify the data
    }

    auto previous_step_data = container.GetDataSpan(PRESSURE, value_policy, StepCategory::TimeStep, 1);
    KRATOS_EXPECT_EQ(previous_step_data.size(), 256); // default chunk size
    for (auto& r_value : previous_step_data) {
        KRATOS_EXPECT_EQ(r_value, 0.0); // previous step should still be zero
    }

    for (auto& r_value : previous_step_data) {
        r_value = 2.0; // modify the data again
    }

    auto current_accessed = container.GetDataSpan(double_accessor);
    KRATOS_EXPECT_EQ(current_accessed.size(), 256); // default chunk size
    for (auto& r_value : current_accessed) {
        KRATOS_EXPECT_EQ(r_value, 1.0); // check that the data was modified
    }

    auto previous_accessed = container.GetDataSpan(double_accessor.GetStepAccessor(StepCategory::TimeStep, 1));
    KRATOS_EXPECT_EQ(previous_accessed.size(), 256); // default chunk size
    for (auto& r_value : previous_accessed) {
        KRATOS_EXPECT_EQ(r_value, 2.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerCloneStep, KratosCoreFastSuite)
{
    DataContainer container;
    DataValuePolicy<double> value_policy;
    HistoricalDataPolicy historical_policy(StepCategory::TimeStep, 3); // Storing 3 time steps including the current one

    auto current_data_accessor = container.Add(PRESSURE, value_policy, historical_policy);
    auto previous_data_accessor = current_data_accessor.GetStepAccessor(StepCategory::TimeStep, 1);

    // initialize current step data
    auto current_step_data = container.GetDataSpan(current_data_accessor);
    for (auto& r_value : current_step_data) {
        r_value = 1.0;
    }

    // initialize previous step data
    auto previous_step_data = container.GetDataSpan(previous_data_accessor);
    for (auto& r_value : previous_step_data) {
        r_value = 2.0;
    }

    container.CloneStepData(StepCategory::TimeStep);

    // now both the current and the previous step hold the pre-clone current data
    auto new_current_step_data = container.GetDataSpan(current_data_accessor);
    for (auto& r_value : new_current_step_data) {
        KRATOS_EXPECT_EQ(r_value, 1.0); // new current step is a copy of the pre-clone current step
    }

    auto new_previous_step_data = container.GetDataSpan(previous_data_accessor);
    for (auto& r_value : new_previous_step_data) {
        KRATOS_EXPECT_EQ(r_value, 1.0); // the pre-clone current step is now the previous one
    }

    auto two_steps_ago_data = container.GetDataSpan(current_data_accessor.GetStepAccessor(StepCategory::TimeStep, 2));
    for (auto& r_value : two_steps_ago_data) {
        KRATOS_EXPECT_EQ(r_value, 2.0); // the pre-clone previous step is now two steps ago
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerAccessorRetrieval, KratosCoreFastSuite)
{
    DataContainer container;

    auto original_double_accessor = container.Add(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy());
    container.Add(STEP, LayeredDataValuePolicy<int>(10), HistoricalDataPolicy(StepCategory::TimeStep, 3));

    for (auto& r_value : container.GetDataSpan(original_double_accessor)) {
        r_value = 3.14;
    }

    // Correct data
    KRATOS_EXPECT_TRUE(container.Has(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy()));

    // For a valid combination, we can get a new accessor to existing data
    auto new_accessor = container.GetAccessor(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy());
    for (const auto& r_value : container.GetDataSpan(new_accessor)) {
        KRATOS_EXPECT_EQ(r_value, 3.14);
    }

    // Variable mismatch
    KRATOS_EXPECT_FALSE(container.Has(TEMPERATURE, DataValuePolicy<double>(), NonHistoricalDataPolicy()));

    // Value policy mismatch
    KRATOS_EXPECT_FALSE(container.Has(PRESSURE, LayeredDataValuePolicy<double>(3), NonHistoricalDataPolicy()));

    // History policy mismatch
    KRATOS_EXPECT_FALSE(container.Has(PRESSURE, DataValuePolicy<double>(), HistoricalDataPolicy(StepCategory::TimeStep, 3)));

    // We do not need to match the policy details (layer number, num steps), knowing the types is enough
    KRATOS_EXPECT_TRUE(container.Has(STEP, LayeredDataValuePolicy<int>(1), HistoricalDataPolicy(StepCategory::TimeStep, 2)));

    // GetAccessor for data not in the container raises an error
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        container.GetAccessor(TEMPERATURE, DataValuePolicy<double>(), NonHistoricalDataPolicy()),
        "No data found in container for TEMPERATURE");
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerVariableLifeTime, KratosCoreFastSuite)
{
    // Ensure that variables defined on the spot can be used to define stored data
    // and that the data is still retrievable even if the original variable goes out of scope
    DataContainer container;

    // This function defines and registers a local variable and uses it to store data in the container
    AddDataForLocalVariable(container);

    // Even though the original variable has gone out of scope, we can still access the associated data
    Variable<double> local_var{"DATA_CONTAINER_TEST_LOCAL"};
    KRATOS_EXPECT_TRUE(container.Has(local_var));

    auto accessor = container.GetAccessor(local_var, DataValuePolicy<double>());
    for (const auto& r_value : container.GetDataSpan(accessor)) {
        KRATOS_EXPECT_EQ(r_value, 5.7);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerAddErrors, KratosCoreFastSuite)
{
    DataContainer container;

    // Value policy incompatible with the variable type
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        container.Add(STEP, DataValuePolicy<double>(), NonHistoricalDataPolicy()),
        "STEP is not compatible with the given value policy.");

    // Unregistered variables are rejected (the container canonicalizes through the Registry)
    Variable<double> unregistered("DATA_CONTAINER_TEST_UNREGISTERED");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        container.Add(unregistered, DataValuePolicy<double>(), NonHistoricalDataPolicy()),
        "DATA_CONTAINER_TEST_UNREGISTERED is not registered in the Registry.");

    // Re-adding with the same policy returns an accessor to the existing chunk
    auto first_accessor = container.Add(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy());
    auto second_accessor = container.Add(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy());
    KRATOS_EXPECT_TRUE(first_accessor == second_accessor);

    // Re-adding with a different policy configuration is an error
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        container.Add(PRESSURE, DataValuePolicy<double>(5.0), NonHistoricalDataPolicy()),
        "Variable PRESSURE exists with a different policy.");
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerForeignAccessorErrors, KratosCoreFastSuite)
{
    DataContainer container;
    container.Add(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy());

    // Accessor pointing beyond the chunks of this container
    DataAccessor<double> out_of_bounds_accessor(PRESSURE, 5, StepCategory::AnyStep);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        container.GetDataSpan(out_of_bounds_accessor),
        "Accessor index out of bounds");

    // Accessor of another variable at an existing chunk index
    DataAccessor<double> mismatched_accessor(TEMPERATURE, 0, StepCategory::AnyStep);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        container.GetDataSpan(mismatched_accessor),
        "Variable mismatch in DataContainer.");
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerCopyIsShallow, KratosCoreFastSuite)
{
    DataContainer container;
    auto accessor = container.Add(PRESSURE, DataValuePolicy<double>(), NonHistoricalDataPolicy());

    // Copies share the chunks with the original
    DataContainer copy(container);
    KRATOS_EXPECT_TRUE(copy.Has(PRESSURE));

    for (auto& r_value : container.GetDataSpan(accessor)) {
        r_value = 4.2;
    }
    for (const auto& r_value : copy.GetDataSpan(accessor)) {
        KRATOS_EXPECT_EQ(r_value, 4.2);
    }

    // Copy assignment behaves the same
    DataContainer assigned;
    assigned = container;
    for (const auto& r_value : assigned.GetDataSpan(accessor)) {
        KRATOS_EXPECT_EQ(r_value, 4.2);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataContainerInitialize, KratosCoreFastSuite)
{
    DataContainer container;
    auto index_accessor = container.Add(STEP, DataValuePolicy<int>(-1), NonHistoricalDataPolicy());
    container.Add(PRESSURE, SparseDataValuePolicy<double>(index_accessor), NonHistoricalDataPolicy());
    container.Add(TEMPERATURE, DataValuePolicy<double>(), NonHistoricalDataPolicy());

    // Initialize mirrors the chunk structure with fresh, zero-initialized chunks:
    // dense chunks at the new chunk size, sparse chunks empty
    DataContainer initialized;
    initialized.Initialize(container, 64);

    KRATOS_EXPECT_TRUE(initialized.Has(STEP));
    KRATOS_EXPECT_TRUE(initialized.Has(PRESSURE));
    KRATOS_EXPECT_TRUE(initialized.Has(TEMPERATURE));

    auto index_span = initialized.GetDataSpan(index_accessor);
    KRATOS_EXPECT_EQ(index_span.size(), 64);
    for (const auto& r_value : index_span) {
        KRATOS_EXPECT_EQ(r_value, -1); // zero-initialized with the policy zero
    }

    auto sparse_accessor = initialized.GetAccessor(PRESSURE, SparseDataValuePolicy<double>(index_accessor));
    KRATOS_EXPECT_EQ(initialized.GetDataSpan(sparse_accessor).size(), 0);

    // The chunk size argument is stored: chunks added afterwards use it
    auto late_accessor = initialized.Add(TEMPERATURE_OLD_IT, DataValuePolicy<double>(), NonHistoricalDataPolicy());
    KRATOS_EXPECT_EQ(initialized.GetDataSpan(late_accessor).size(), 64);

    // A zero chunk size is rejected
    DataContainer invalid;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        invalid.Initialize(container, 0),
        "Chunk size must be greater than zero.");
}

} // namespace Testing
} // namespace Kratos
