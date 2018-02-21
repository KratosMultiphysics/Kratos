#include "tests/test_utils.h"

#include "containers/array_1d.h"

namespace Kratos
{
namespace Testing
{

template<>
HDF5::File::Vector<array_1d<double,3>> TestVector(std::size_t n)
{
    HDF5::File::Vector<array_1d<double,3>> vec(n);
    for (std::size_t i = 0; i < n; ++i)
        vec(i) = array_1d<double,3>(3, i);
    return vec;
}

} // namespace Testing
} // namespace Kratos.
