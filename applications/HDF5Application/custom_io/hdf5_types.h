//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

#pragma once

// --- HDF5 Includes ---
#include "hdf5_application_define.h"


namespace Kratos::HDF5 {


template <class TValue>
struct HDF5ScalarTypeID {};


template <>
struct HDF5ScalarTypeID<int>
{static inline hid_t value = H5T_NATIVE_INT;};


template <>
struct HDF5ScalarTypeID<double>
{static inline hid_t value = H5T_NATIVE_DOUBLE;};


template <>
struct HDF5ScalarTypeID<Vector<int>>
{static inline hid_t value = H5T_NATIVE_INT;};


template <>
struct HDF5ScalarTypeID<Vector<double>>
{static inline hid_t value = H5T_NATIVE_DOUBLE;};


template <>
struct HDF5ScalarTypeID<Vector<array_1d<double, 3>>>
{static inline hid_t value = H5T_NATIVE_DOUBLE;};


template <>
struct HDF5ScalarTypeID<Matrix<int>>
{static inline hid_t value = H5T_NATIVE_INT;};


template <>
struct HDF5ScalarTypeID<Matrix<double>>
{static inline hid_t value = H5T_NATIVE_DOUBLE;};


template <class TValue>
const inline hid_t HDF5ScalarTypeID_v = HDF5ScalarTypeID<TValue>::value;


} // namespace Kratos::HDF5
