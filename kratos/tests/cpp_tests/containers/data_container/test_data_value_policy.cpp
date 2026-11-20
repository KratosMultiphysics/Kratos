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
#include "containers/data_container/data_value_policy.h"
#include "containers/data_container/sparse_data_value_policy.h"
#include "includes/expect.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DataValuePolicyComparison, KratosCoreFastSuite)
{
    DataValuePolicy<double> double_value_policy;

    KRATOS_EXPECT_TRUE(DataValuePolicy<double>().IsSameType(double_value_policy));
    KRATOS_EXPECT_FALSE(DataValuePolicy<int>().IsSameType(double_value_policy));

    LayeredDataValuePolicy<double> double_layer_policy(3);

    KRATOS_EXPECT_FALSE(double_layer_policy.IsSameType(double_value_policy));

    LayeredDataValuePolicy<double> other_layer_policy(10);
    KRATOS_EXPECT_TRUE(double_layer_policy.IsSameType(other_layer_policy)); // same type
    KRATOS_EXPECT_FALSE(double_layer_policy.IsSame(other_layer_policy)); // different number of layers
}

KRATOS_TEST_CASE_IN_SUITE(DataValuePolicyIsCompatible, KratosCoreFastSuite)
{
    // Real core variables of the three value types
    const auto& r_int_variable = STEP;
    const auto& r_double_variable = PRESSURE;
    const auto& r_array_variable = VELOCITY;

    DataValuePolicy<int> int_value_policy;

    KRATOS_EXPECT_TRUE(int_value_policy.IsCompatible(r_int_variable));
    KRATOS_EXPECT_FALSE(int_value_policy.IsCompatible(r_double_variable));
    KRATOS_EXPECT_FALSE(int_value_policy.IsCompatible(r_array_variable));

    DataValuePolicy<double> double_value_policy;

    KRATOS_EXPECT_FALSE(double_value_policy.IsCompatible(r_int_variable));
    KRATOS_EXPECT_TRUE(double_value_policy.IsCompatible(r_double_variable));
    KRATOS_EXPECT_FALSE(double_value_policy.IsCompatible(r_array_variable));

    DataValuePolicy<array_1d<double, 3>> array_value_policy;

    KRATOS_EXPECT_FALSE(array_value_policy.IsCompatible(r_int_variable));
    KRATOS_EXPECT_FALSE(array_value_policy.IsCompatible(r_double_variable));
    KRATOS_EXPECT_TRUE(array_value_policy.IsCompatible(r_array_variable));

    LayeredDataValuePolicy<int> int_layer_policy(3);

    KRATOS_EXPECT_TRUE(int_layer_policy.IsCompatible(r_int_variable));
    KRATOS_EXPECT_FALSE(int_layer_policy.IsCompatible(r_double_variable));
    KRATOS_EXPECT_FALSE(int_layer_policy.IsCompatible(r_array_variable));

    LayeredDataValuePolicy<double> double_layer_policy(3);

    KRATOS_EXPECT_FALSE(double_layer_policy.IsCompatible(r_int_variable));
    KRATOS_EXPECT_TRUE(double_layer_policy.IsCompatible(r_double_variable));
    KRATOS_EXPECT_FALSE(double_layer_policy.IsCompatible(r_array_variable));

    LayeredDataValuePolicy<array_1d<double, 3>> array_layer_policy(3);

    KRATOS_EXPECT_FALSE(array_layer_policy.IsCompatible(r_int_variable));
    KRATOS_EXPECT_FALSE(array_layer_policy.IsCompatible(r_double_variable));
    KRATOS_EXPECT_TRUE(array_layer_policy.IsCompatible(r_array_variable));
}

KRATOS_TEST_CASE_IN_SUITE(DataValuePolicyZero, KratosCoreFastSuite)
{
    // Default zero values
    KRATOS_EXPECT_EQ(DataValuePolicy<double>().Zero(), 0.0);
    KRATOS_EXPECT_EQ(DataValuePolicy<int>().Zero(), 0);
    KRATOS_EXPECT_EQ(DataValuePolicy<std::string>().Zero(), std::string(""));

    // Custom zero values (as used for sparse index variables)
    DataValuePolicy<int> minus_one_policy(-1);
    KRATOS_EXPECT_EQ(minus_one_policy.Zero(), -1);

    // Policies with different zeros are the same type but not the same policy
    KRATOS_EXPECT_TRUE(minus_one_policy.IsSameType(DataValuePolicy<int>()));
    KRATOS_EXPECT_FALSE(minus_one_policy.IsSame(DataValuePolicy<int>()));
    KRATOS_EXPECT_TRUE(minus_one_policy.IsSame(DataValuePolicy<int>(-1)));

    // The layered policy zero is a vector of per-layer zeros
    LayeredDataValuePolicy<double> layer_policy(4, 1.5);
    KRATOS_EXPECT_EQ(layer_policy.Zero(), std::vector<double>(4, 1.5));
    KRATOS_EXPECT_EQ(layer_policy.NumberOfLayers(), 4);
}

KRATOS_TEST_CASE_IN_SUITE(DataValuePolicyRawOperations, KratosCoreFastSuite)
{
    // Exercise the type-erased operations through the base class interface
    DataValuePolicy<double> typed_policy(2.5);
    const DataValuePolicyBase& r_policy = typed_policy;

    KRATOS_EXPECT_EQ(r_policy.SizeOfData(), sizeof(double));
    KRATOS_EXPECT_FALSE(r_policy.IsSparse());

    // Allocate a default-constructed value
    void* p_data = nullptr;
    r_policy.Allocate(&p_data);
    KRATOS_EXPECT_TRUE(p_data != nullptr);

    // AssignZero constructs the policy zero in place
    r_policy.AssignZero(p_data);
    KRATOS_EXPECT_EQ(*static_cast<double*>(p_data), 2.5);

    // Assign copies an existing value
    const double source = 7.25;
    r_policy.Assign(&source, p_data);
    KRATOS_EXPECT_EQ(*static_cast<double*>(p_data), 7.25);

    // CloneData allocates a copy
    void* p_clone = r_policy.CloneData(p_data);
    KRATOS_EXPECT_EQ(*static_cast<double*>(p_clone), 7.25);

    // Copy placement-constructs into already allocated memory
    double destination = 0.0;
    r_policy.Copy(&source, &destination);
    KRATOS_EXPECT_EQ(destination, 7.25);

    // Print writes the plain value, PrintData adds the label
    std::stringstream print_stream;
    r_policy.Print(p_data, print_stream);
    KRATOS_EXPECT_EQ(print_stream.str(), "7.25");

    std::stringstream print_data_stream;
    r_policy.PrintData(p_data, print_data_stream);
    KRATOS_EXPECT_EQ(print_data_stream.str(), "Data: 7.25");

    r_policy.Delete(p_data);
    r_policy.Delete(p_clone);
}

KRATOS_TEST_CASE_IN_SUITE(SparseDataValuePolicyBasics, KratosCoreFastSuite)
{
    DataValuePolicy<double> dense_policy;
    KRATOS_EXPECT_FALSE(dense_policy.IsSparse());

    // The base class refuses to provide a sparse index accessor
    const DataValuePolicyBase& r_dense_base = dense_policy;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        r_dense_base.GetSparseIndexAccessor(),
        "Calling base class DataValuePolicyBase::GetSparseIndexAccessor");

    // A sparse policy reports IsSparse and returns its index accessor
    DataAccessor<int> index_accessor(STEP, 0, StepCategory::AnyStep);
    SparseDataValuePolicy<double> sparse_policy(index_accessor);

    KRATOS_EXPECT_TRUE(sparse_policy.IsSparse());
    KRATOS_EXPECT_TRUE(sparse_policy.GetSparseIndexAccessor() == index_accessor);

    // SparseDataValuePolicy does not override IsSameType: sparse and dense policies of the
    // same value type compare as the same type (sparseness is queried via IsSparse instead)
    KRATOS_EXPECT_TRUE(sparse_policy.IsSameType(dense_policy));
    KRATOS_EXPECT_TRUE(dense_policy.IsSameType(sparse_policy));
    KRATOS_EXPECT_TRUE(sparse_policy.IsSameType(SparseDataValuePolicy<double>(index_accessor)));
    KRATOS_EXPECT_FALSE(sparse_policy.IsSameType(DataValuePolicy<int>()));

    // Clone preserves sparseness and the index accessor
    auto p_clone = sparse_policy.Clone();
    KRATOS_EXPECT_TRUE(p_clone->IsSparse());
    KRATOS_EXPECT_TRUE(p_clone->GetSparseIndexAccessor() == index_accessor);
}

KRATOS_TEST_CASE_IN_SUITE(DataAccessorEquality, KratosCoreFastSuite)
{
    // Default-constructed accessors compare equal to each other
    KRATOS_EXPECT_TRUE(DataAccessor<double>() == DataAccessor<double>());
    KRATOS_EXPECT_EQ(DataAccessor<double>().pGetVariableData(), nullptr);

    DataAccessor<double> accessor(PRESSURE, 3, StepCategory::TimeStep, 0);
    KRATOS_EXPECT_EQ(accessor.GetIndex(), 3);
    KRATOS_EXPECT_TRUE(accessor.GetStepCategory() == StepCategory::TimeStep);
    KRATOS_EXPECT_EQ(accessor.GetStepBeforeCurrent(), 0);
    KRATOS_EXPECT_TRUE(accessor.GetVariableData() == PRESSURE);

    // Same variable, index and step: equal
    KRATOS_EXPECT_TRUE(accessor == DataAccessor<double>(PRESSURE, 3, StepCategory::TimeStep, 0));

    // Different index, variable, category or step: not equal
    KRATOS_EXPECT_FALSE(accessor == DataAccessor<double>(PRESSURE, 4, StepCategory::TimeStep, 0));
    KRATOS_EXPECT_FALSE(accessor == DataAccessor<double>(TEMPERATURE, 3, StepCategory::TimeStep, 0));
    KRATOS_EXPECT_FALSE(accessor == DataAccessor<double>(PRESSURE, 3, StepCategory::IterationStep, 0));
    KRATOS_EXPECT_FALSE(accessor == DataAccessor<double>(PRESSURE, 3, StepCategory::TimeStep, 1));

    // GetStepAccessor re-parameterizes step information keeping variable and index
    auto previous_step_accessor = accessor.GetStepAccessor(StepCategory::TimeStep, 1);
    KRATOS_EXPECT_EQ(previous_step_accessor.GetIndex(), 3);
    KRATOS_EXPECT_EQ(previous_step_accessor.GetStepBeforeCurrent(), 1);
    KRATOS_EXPECT_TRUE(previous_step_accessor.GetVariableData() == PRESSURE);
    KRATOS_EXPECT_FALSE(previous_step_accessor == accessor);
    KRATOS_EXPECT_TRUE(previous_step_accessor.GetStepAccessor(StepCategory::TimeStep, 0) == accessor);
}

} // namespace Testing
} // namespace Kratos
