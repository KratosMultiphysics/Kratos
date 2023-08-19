//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Suneth Warnakulasuriya
//

// System includes
#include "hdf5.h"
#include <type_traits>

// Project includes
#include "utilities/data_type_traits.h"

// Application includes


namespace Kratos
{
namespace HDF5
{
namespace Internals
{

// H5 data types
template <class TDataType> hid_t GetPrimitiveH5Type()
{
    using primitive_type = typename DataTypeTraits<TDataType>::PrimitiveType;

    if constexpr(std::is_same_v<primitive_type, char>) {
        return H5T_NATIVE_CHAR;
    } else if constexpr(std::is_same_v<primitive_type, int>) {
        return H5T_NATIVE_INT;
    } else if constexpr(std::is_same_v<primitive_type, double>) {
        return H5T_NATIVE_DOUBLE;
    } else if constexpr(std::is_same_v<primitive_type, hsize_t>) {
        return H5T_NATIVE_HSIZE;
    } else {
        static_assert(!std::is_same_v<primitive_type, primitive_type>, "Unsupported data type.");
    }
}

} // namespace Internals
} // namespace HDF5
} // namespace Kratos