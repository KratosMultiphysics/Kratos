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
#include <unordered_map>
#include <sstream>

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

    if constexpr(std::is_same_v<primitive_type, bool>) {
        return H5T_NATIVE_HBOOL;
    } else if constexpr(std::is_same_v<primitive_type, char>) {
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

template<class TContainerType>
[[nodiscard]] std::vector<std::string> GetListOfAvailableVariables(
    TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    using map_type = std::unordered_map<VariableData::KeyType, std::string>;

    // Map to Map reduction class
    class MapReduction
    {
    public:
        using return_type = map_type;

        return_type mValue;

        /// access to reduced value
        return_type GetValue() const
        {
            return mValue;
        }

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(map_type rValue){
            mValue.merge(rValue);
        }

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(MapReduction& rOther)
        {
            KRATOS_CRITICAL_SECTION
            mValue.merge(rOther.mValue);
        }
    };

    // collect all the unique key, variable names map from all the entities
    const auto& key_variable_names_map = block_for_each<MapReduction>(rContainer, map_type(), [](auto& rEntity, auto& rTLS) {
        const auto& r_data = rEntity.GetData();
        for (auto it = r_data.begin(); it != r_data.end(); ++it) {
            auto itr = rTLS.find(it->first->Key());
            if (itr == rTLS.end()) {
                rTLS.insert(std::make_pair(it->first->Key(), it->first->Name()));
            }
        }

        return rTLS;
    });

    // serialize the collected map with only the variable names
    std::vector<char> serialized_char_array;
    std::for_each(key_variable_names_map.begin(), key_variable_names_map.end(),
                  [&serialized_char_array](const auto& rPair) {
                      std::for_each(rPair.second.begin(), rPair.second.end(),
                                    [&serialized_char_array](const auto CurrentChar) {
                                        serialized_char_array.push_back(CurrentChar);
                                    });
                      serialized_char_array.push_back(';');
                  });

    // do the mpi communication
    const auto& global_serialized_char_array = rDataCommunicator.AllGatherv(serialized_char_array);

    std::vector<std::string> result;
    for (const auto& r_serialized_char_array : global_serialized_char_array) {
        std::stringstream temp;
        for (const auto current_char : r_serialized_char_array) {
            if (current_char != ';') {
                temp << current_char;
            } else {
                result.push_back(temp.str());
                temp.str(std::string()); // clear the string stream
            }
        }
    }

    std::sort(result.begin(), result.end());
    auto last = std::unique(result.begin(), result.end());
    result.erase(last, result.end());
    return result;

    KRATOS_CATCH("");
}

template<class TDataType>
struct ComponentTraits {};

template<> struct ComponentTraits<Flags> { using ValueType = char; };
template<class TDataType> struct ComponentTraits<Variable<TDataType>> { using ValueType = TDataType; };

} // namespace Internals
} // namespace HDF5
} // namespace Kratos