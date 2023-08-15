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
template <class TDataType> hid_t inline GetPrimitiveH5Type()
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
    if (Attributes.Has("__shape") && Attributes["__shape"].IsVector()) {
        const auto& shape_param = Attributes["__shape"];
        rShape.resize(shape_param.size());
        for (unsigned int i = 0; i < shape_param.size(); ++i) {
            KRATOS_ERROR_IF_NOT(shape_param.GetArrayItem(i).IsInt())
                << "Invalid type found in \"__shape\" attribute.";
            rShape[i] = shape_param.GetArrayItem(i).GetInt();
        }
    } else {
        LegacyGetShapeFromAttributes(rShape, Attributes);
    }
}

template<class TContainerType>
std::string GetContainerType()
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return "NODES";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return "CONDITIONS";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return "ELEMENTS";
    } else {
        static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
    }
}

struct FlagIO
{
    static constexpr std::string_view mInfo = "Flag";

    template<class TEntityType>
    static char GetValue(
        const TEntityType& rEntity,
        const Flags& rFlag)
    {
        return rEntity.Is(rFlag);
    }

    template<class TEntityType>
    static void SetValue(
        TEntityType& rEntity,
        const Flags& rVariable,
        const char rValue)
    {
        rEntity.Set(rVariable, static_cast<bool>(rValue));
    }
};

template<class TContainerIOType>
std::string GetContainerIOType()
{
    if constexpr(std::is_same_v<TContainerIOType, ContainerDataIO<ContainerDataIOTags::Historical>>) {
        return "HISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
        return "NONHISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, FlagIO>) {
        return "FLAGS";
    } else {
        static_assert(!std::is_same_v<TContainerIOType, TContainerIOType>, "Unsupported container io type.");
    }
}

template<class TDataType>
struct ComponentTraits {};

template<> struct ComponentTraits<Flags> { using ValueType = char; };
template<class TDataType> struct ComponentTraits<Variable<TDataType>> { using ValueType = TDataType; };

template<class TContainerDataIO, class TContainerType, class TComponentType, class TDataType>
void CopyToContiguousDataArray(
    const TContainerType& rContainer,
    const TComponentType& rComponent,
    TDataType* pBegin)
{
    KRATOS_TRY

    using value_type_traits = DataTypeTraits<typename ComponentTraits<TComponentType>::ValueType>;

    static_assert(value_type_traits::IsContiguous, "Only contiguous data types are supported.");

    if (rContainer.empty()) {
        // do nothing if the container is empty.
        return;
    }

    // get the first item for sizing.
    const auto& initial_value = TContainerDataIO::GetValue(rContainer.front(), rComponent);

    // get the stride from the first element to support dynamic types.
    const auto stride = value_type_traits::Size(initial_value);

    IndexPartition<unsigned int>(rContainer.size()).for_each([&rContainer, &rComponent, pBegin, stride](const auto Index) {
        const auto& value = TContainerDataIO::GetValue(*(rContainer.begin() + Index), rComponent);
        TDataType const* p_value_begin = value_type_traits::GetContiguousData(value);
        auto p_array_start = pBegin + Index * stride;
        std::copy(p_value_begin, p_value_begin + stride, p_array_start);
    });

    KRATOS_CATCH("");
}

template<class TContainerDataIO, class TContainerType, class TComponentType, class TDataType>
void CopyFromContiguousDataArray(
    TContainerType& rContainer,
    const TComponentType& rComponent,
    TDataType const* pBegin,
    const std::vector<unsigned int>& rShape)
{
    KRATOS_TRY

    using value_type_traits = DataTypeTraits<typename ComponentTraits<TComponentType>::ValueType>;

    static_assert(value_type_traits::IsContiguous, "Only contiguous data types are supported.");

    using value_type = typename value_type_traits::ContainerType;

    if (rContainer.empty()) {
        // do nothing if the container is empty.
        return;
    }

    value_type tls_prototype;
    value_type_traits::Reshape(tls_prototype, rShape);

    const auto stride = value_type_traits::Size(tls_prototype);

    IndexPartition<unsigned int>(rContainer.size()).for_each(tls_prototype, [&rContainer, &rComponent, pBegin, stride](const auto Index, auto& rTLS) {
        TDataType * p_value_begin = value_type_traits::GetContiguousData(rTLS);
        TDataType const * p_array_start = pBegin + Index * stride;
        std::copy(p_array_start, p_array_start + stride, p_value_begin);
        TContainerDataIO::SetValue(*(rContainer.begin() + Index), rComponent, rTLS);
    });

    KRATOS_CATCH("");
}

} // namespace Internals
} // namespace HDF5
} // namespace Kratos