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

#pragma once

// System includes
#include "hdf5.h"
#include <algorithm>
#include <type_traits>

// Project includes
#include "containers/flags.h"
#include "containers/variable.h"
#include "expression/container_data_io.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"

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

template<class TIndexType>
void LegacyGetShapeFromAttributes(
    std::vector<TIndexType>& rShape,
    const Parameters Attributes)
{
    if (Attributes.Has("Size1") && Attributes.Has("Size2") && Attributes["Size1"].IsInt() && Attributes["Size2"].IsInt()) {
        rShape.resize(2);
        rShape[0] = Attributes["Size1"].GetInt();
        rShape[1] = Attributes["Size2"].GetInt();
    } else {
        KRATOS_ERROR << "Shape information not found in attributes.";
    }
}

template<class TIndexType>
void GetShapeFromAttributes(
    std::vector<TIndexType>& rShape,
    const Parameters Attributes)
{
    if (Attributes.Has("__data_shape") && Attributes["__data_shape"].IsVector()) {
        const auto& shape_param = Attributes["__data_shape"];
        rShape.resize(shape_param.size());
        for (unsigned int i = 0; i < shape_param.size(); ++i) {
            KRATOS_ERROR_IF_NOT(shape_param.GetArrayItem(i).IsInt())
                << "Invalid type found in \"__data_shape\" attribute.";
            rShape[i] = shape_param.GetArrayItem(i).GetInt();
        }
    } else {
        LegacyGetShapeFromAttributes(rShape, Attributes);
    }
}



template<class TDataType>
struct ComponentTraits {};

template<> struct ComponentTraits<Flags> { using ValueType = char; };
template<class TDataType> struct ComponentTraits<Variable<TDataType>> { using ValueType = TDataType; };

} // namespace Internals
} // namespace HDF5
} // namespace Kratos