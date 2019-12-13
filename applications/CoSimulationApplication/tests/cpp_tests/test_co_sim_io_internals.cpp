// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
#include "co_simulation_io/impl/co_sim_io_internals.h"

namespace Kratos {
namespace Testing {

typedef double ValueType;
typedef CoSimIO::Internals::DataContainer<ValueType> DataContainerBase;
typedef Kratos::unique_ptr<DataContainerBase> DataContainerBasePointer;
typedef CoSimIO::Internals::DataContainerStdVector<ValueType> DataContainerStdVectorType;
typedef CoSimIO::Internals::DataContainerRawMemory<ValueType> DataContainerRawMemoryType;

namespace {

void TestDataContainerBasics(const std::vector<ValueType>& rRefValues, DataContainerBase& rDataContainer)
{
    // checking size
    KRATOS_CHECK_EQUAL(rRefValues.size(), rDataContainer.size());

    // checking values
    KRATOS_CHECK_VECTOR_NEAR(rRefValues, rDataContainer, 1e-12)
}

void TestDataContainerDifferentValues(const std::vector<std::vector<ValueType>>& rRefValues, DataContainerBase& rDataContainer)
{
    for (const auto& r_current_ref_vals : rRefValues) {
        const std::size_t current_size(r_current_ref_vals.size());
        if (rDataContainer.size() != current_size) {
            rDataContainer.resize(current_size);
        }

        for (std::size_t i=0; i<current_size; ++i) {
            rDataContainer[i] = r_current_ref_vals[i];
        }

        // checking size
        KRATOS_CHECK_EQUAL(current_size, rDataContainer.size());

        // checking values
        KRATOS_CHECK_VECTOR_NEAR(r_current_ref_vals, rDataContainer, 1e-12)
    }
}

} // helpers namespace

KRATOS_TEST_CASE_IN_SUITE(DataContainers_RawMemory_empty, KratosCoSimulationFastSuite)
{
    ValueType** values_raw = (ValueType**)malloc(sizeof(ValueType*)*1);
    values_raw[0] = NULL;

    DataContainerBasePointer p_container(Kratos::make_unique<DataContainerRawMemoryType>(values_raw, 0));

    KRATOS_CHECK_EQUAL(0, p_container->size());

    // deallocating memory
    free(*values_raw);
    free(values_raw);
}

KRATOS_TEST_CASE_IN_SUITE(DataContainers_StdVector_basics, KratosCoSimulationFastSuite)
{
    const std::vector<ValueType> ref_values {
        1.0, -2.333, 15.88, 14.7, -99.6
    };

    std::vector<ValueType> values_vec(ref_values);

    DataContainerBasePointer p_container(Kratos::make_unique<DataContainerStdVectorType>(values_vec));

    TestDataContainerBasics(ref_values, *p_container);
}

KRATOS_TEST_CASE_IN_SUITE(DataContainers_RawMemory_basics, KratosCoSimulationFastSuite)
{
    const std::vector<ValueType> ref_values {
        1.0, -2.333, 15.88, 14.7, -99.6
    };

    const std::size_t size(ref_values.size());

    ValueType** values_raw = (ValueType**)malloc(sizeof(ValueType*)*1);
    values_raw[0]= (ValueType*)malloc(sizeof(ValueType)*size);

    for (std::size_t i=0; i<size; ++i) {
        (*values_raw)[i] = ref_values[i];
    }

    DataContainerBasePointer p_container(Kratos::make_unique<DataContainerRawMemoryType>(values_raw, size));

    TestDataContainerBasics(ref_values, *p_container);
}

KRATOS_TEST_CASE_IN_SUITE(DataContainers_StdVector_multiple_resizes, KratosCoSimulationFastSuite)
{
    const std::vector<std::vector<ValueType>> ref_values {
        {1.0, -2.333, 15.88, 14.7, -99.6},
        {2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.5},
        {-88.66, 77.9}
    };

    std::vector<ValueType> values_vec;

    DataContainerBasePointer p_container(Kratos::make_unique<DataContainerStdVectorType>(values_vec));

    TestDataContainerDifferentValues(ref_values, *p_container);
}

KRATOS_TEST_CASE_IN_SUITE(DataContainers_RawMemory_resize_larger, KratosCoSimulationFastSuite)
{
    const std::vector<std::vector<ValueType>> ref_values {
        {1.0, -2.333, 15.88, 14.7, -99.6},
        {2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.5}
    };

    ValueType** values_raw = (ValueType**)malloc(sizeof(ValueType*)*1);
    values_raw[0] = NULL;
    DataContainerBasePointer p_container(Kratos::make_unique<DataContainerRawMemoryType>(values_raw, 0));

    TestDataContainerDifferentValues(ref_values, *p_container);

    // deallocating memory
    free(*values_raw);
    free(values_raw);
}

KRATOS_TEST_CASE_IN_SUITE(DataContainers_RawMemory_resize_smaller, KratosCoSimulationFastSuite)
{
    const std::vector<std::vector<ValueType>> ref_values {
        {1.0, -2.333, 15.88, 14.7, -99.6},
        {-88.66, 77.9}
    };

    ValueType** values_raw = (ValueType**)malloc(sizeof(ValueType*)*1);
    values_raw[0] = NULL;
    DataContainerBasePointer p_container(Kratos::make_unique<DataContainerRawMemoryType>(values_raw, 0));

    TestDataContainerDifferentValues(ref_values, *p_container);

    // deallocating memory
    free(*values_raw);
    free(values_raw);
}

KRATOS_TEST_CASE_IN_SUITE(DataContainers_RawMemory_multiple_resizes, KratosCoSimulationFastSuite)
{
    const std::vector<std::vector<ValueType>> ref_values {
        {1.0, -2.333, 15.88, 14.7, -99.6},
        {-88.66, 77.9},
        {1.0, -2.333, 15.88, 14.7, -99.6},
        {2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.5},
        {1.0, -2.333, 15.88, 14.7, -99.6},
        {-88.66, 77.9}
    };

    ValueType** values_raw = (ValueType**)malloc(sizeof(ValueType*)*1);
    values_raw[0] = NULL;
    DataContainerBasePointer p_container(Kratos::make_unique<DataContainerRawMemoryType>(values_raw, 0));

    TestDataContainerDifferentValues(ref_values, *p_container);

    // deallocating memory
    free(*values_raw);
    free(values_raw);
}

} // namespace Testing
}  // namespace Kratos.
