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

// Project includes

// Application includes


namespace Kratos
{
namespace HDF5
{
namespace Internals
{

// H5 data types
template <class TDataType> hid_t inline GetH5DataType() { static_assert(true, "Unsupported data type."); return H5T_NATIVE_INT; }
template <> hid_t inline GetH5DataType<hsize_t>() { return H5T_NATIVE_HSIZE; }
template <> hid_t inline GetH5DataType<int>() { return H5T_NATIVE_INT; }
template <> hid_t inline GetH5DataType<double>() { return H5T_NATIVE_DOUBLE; }

} // namespace Internals
} // namespace HDF5
} // namespace Kratos