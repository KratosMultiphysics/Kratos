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
// #include <cstdlib>
#include <cmath>

// External includes

// Project includes
#include "testing/testing.h"
#include "co_simulation_io/internals/co_sim_io_internals.h"

namespace Kratos {
namespace Testing {

typedef double ValueType;
typedef CoSimIO::Internals::DataContainer<ValueType> DataContainerBase;
typedef Kratos::unique_ptr<DataContainerBase> DataContainerBasePointer;
typedef CoSimIO::Internals::DataContainerStdVector<ValueType> DataContainerStdVectorType;
typedef CoSimIO::Internals::DataContainerRawMemory<ValueType> DataContainerRawMemoryType;

KRATOS_TEST_CASE_IN_SUITE(DataContainers, KratosCoSimulationFastSuite)
{
    const std::size_t init_size(15);
    const double init_val(1.23);

    std::vector<ValueType> ref_values(init_size);
    std::vector<ValueType> values_vec(init_size);
    // double* values_raw = new double[init_size];
    double** values_raw = new ValueType*[1];
    *values_raw = new ValueType[init_size];

    for (std::size_t i=0; i<init_size; ++i) {
        ref_values[i] = i*init_val;
        values_vec[i] = i*init_val;
        (*values_raw)[i] = i*init_val;
    }

    DataContainerBasePointer ptr_std_vec(Kratos::make_unique<DataContainerStdVectorType>(values_vec));
    DataContainerBasePointer ptr_raw_mem(Kratos::make_unique<DataContainerRawMemoryType>(values_raw, init_size));

    // checking size
    KRATOS_CHECK_EQUAL(init_size, ptr_std_vec->size());
    KRATOS_CHECK_EQUAL(init_size, ptr_raw_mem->size());

    // checking values
    KRATOS_CHECK_VECTOR_NEAR(ref_values, (*ptr_std_vec), 1e-12)
    KRATOS_CHECK_VECTOR_NEAR(ref_values, (*ptr_raw_mem), 1e-12)

    // checking resize (increasing size)
    std::size_t size_after_resize(init_size+10);
    ptr_std_vec->resize(size_after_resize);
    ptr_raw_mem->resize(size_after_resize);

    KRATOS_CHECK_EQUAL(size_after_resize, ptr_std_vec->size());
    KRATOS_CHECK_EQUAL(size_after_resize, ptr_raw_mem->size());

    ref_values.resize(size_after_resize);
    for (std::size_t i=0; i<size_after_resize; ++i) {
        ref_values[i] = i*init_val+2.558;
        (*ptr_std_vec)[i] = i*init_val+2.558;
        (*ptr_raw_mem)[i] = i*init_val+2.558;
    }

    // deallocating memory
    delete [] values_raw[0];
    delete [] values_raw;
}

} // namespace Testing
}  // namespace Kratos.
