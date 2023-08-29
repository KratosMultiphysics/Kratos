//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya, https://github.com/sunethwarna
//

// System includes
#include <numeric>
#include <type_traits>

// Project includes
#include "expression/container_data_io.h"

// Application includes
#include "custom_utilities/container_io_utils.h"
#include "custom_utilities/data_type_utilities.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "hdf5_application_variables.h"

// Include base h
#include "custom_io/hdf5_container_component_io.h"

namespace Kratos
{
namespace HDF5
{

namespace NewContainerComponentIOUtilities
{
template<class TContainerType>
struct SynchronizeComponent {
    template<class TContainerDataIO, class TComponent>
    static void Execute(
        Communicator& rCommunicator,
        const TComponent& rComponent)
    {}
};

template<>
struct SynchronizeComponent<ModelPart::NodesContainerType>
{
    template<class TContainerDataIO, class TComponent>
    static void Execute(
        Communicator& rCommunicator,
        const TComponent& rComponent)
    {
        if constexpr(std::is_same_v<TContainerDataIO, Internals::FlagIO>) {
            rCommunicator.SynchronizeOrNodalFlags(rComponent);
        } else if constexpr(std::is_same_v<TContainerDataIO, Internals::HistoricalIO> || std::is_same_v<TContainerDataIO, Internals::BossakIO>) {
            rCommunicator.SynchronizeVariable(rComponent);
        } else if constexpr(std::is_same_v<TContainerDataIO, Internals::NonHistoricalIO>) {
            rCommunicator.SynchronizeNonHistoricalVariable(rComponent);
        } else {
            static_assert(!std::is_same_v<TContainerDataIO, TContainerDataIO>, "Unsupported TContainerDataIO.");
        }
    }
};
} // namespace NewContainerComponentIOUtilities

template <class TContainerType, class TContainerDataIO, class... TComponents>
ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::ContainerComponentIO(
    Parameters Settings,
    File::Pointer pFile,
    const std::string& rLegacySuffix)
    : mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix"           : "",
            "list_of_variables": []
        })");

    Settings.AddMissingParameters(default_params);

    mComponentPrefix = Settings["prefix"].GetString() + rLegacySuffix;
    mComponentNames = Settings["list_of_variables"].GetStringArray();

    KRATOS_ERROR_IF(mComponentPrefix == "" || mComponentPrefix == "/")
        << "The prefix should not be blank or \"/\" [ prefix = " << mComponentPrefix << " ].\n";

    // Sort component names to make sure they're in the same order on each rank.
    // The basic assumption is that the set of components is identical on every
    // rank, but they may not be in the same order (which would lead to ranks
    // trying to write different variables at the same time, resulting in
    // a deadlock), hence the sorting.
    std::sort(mComponentNames.begin(), mComponentNames.end());

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
void ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::Write(
    const ModelPart& rModelPart,
    const TContainerDataIO& rContainerDataIO,
    const Parameters Attributes)
{
    KRATOS_TRY

    auto current_attributes = Attributes.Clone();

    if constexpr(
        std::is_same_v<TContainerType, NodesContainerType> ||
        std::is_same_v<TContainerType, ConditionsContainerType> ||
        std::is_same_v<TContainerType, ElementsContainerType>) {
        if (rModelPart.Has(HDF5_MESH_LOCATION_INFO)) {
            KRATOS_ERROR_IF(current_attributes.Has("__mesh_location"))
                << "The reserved keyword \"__mesh_location\" is found. Please remove it from attributes.";
            current_attributes.AddString("__mesh_location", rModelPart.GetValue(HDF5_MESH_LOCATION_INFO));
        }
        Write(Internals::GetLocalContainer<TContainerType>(rModelPart), rContainerDataIO, current_attributes);
    } else {
        KRATOS_ERROR << "Unsupported container type.";
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
std::map<std::string, Parameters> ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::Read(
    ModelPart& rModelPart,
    const TContainerDataIO& rContainerDataIO)
{
    KRATOS_TRY

    if constexpr(
        std::is_same_v<TContainerType, NodesContainerType> ||
        std::is_same_v<TContainerType, ConditionsContainerType> ||
        std::is_same_v<TContainerType, ElementsContainerType>) {
        return Read(Internals::GetLocalContainer<TContainerType>(rModelPart), rContainerDataIO, rModelPart.GetCommunicator());
    } else {
        KRATOS_ERROR << "Unsupported container type.";
        return std::map<std::string, Parameters>{};
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
void ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::Write(
    const TContainerType& rLocalContainer,
    const TContainerDataIO& rContainerDataIO,
    const Parameters Attributes)
{
    KRATOS_TRY;

    if (mComponentNames.size() == 0) {
        return;
    }

    WriteInfo info;

    // Write each variable.
    for (const auto& r_component_name : mComponentNames) {
        const bool is_component_written = (... || (WriteComponentData<TComponents>(
                        r_component_name, rContainerDataIO, rLocalContainer, Attributes.Clone(), info)));
        KRATOS_ERROR_IF_NOT(is_component_written)
            << "Component \"" << r_component_name << "\" is not found in registered components.";
    }

    // Write block partition.
    WritePartitionTable(*mpFile, mComponentPrefix, info);

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
std::map<std::string, Parameters> ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::Read(
    TContainerType& rLocalContainer,
    const TContainerDataIO& rContainerDataIO,
    Communicator& rCommunicator)
{
    KRATOS_TRY;

    std::map<std::string, Parameters> attributes;

    if (mComponentNames.size() == 0) {
        return attributes;
    }

    IndexType start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, mComponentPrefix);

    // Write each variable.
    for (const auto& r_component_name : mComponentNames) {
        const bool is_component_written =
            (... || (ReadComponentData<TComponents>(
                        r_component_name, rContainerDataIO, rLocalContainer, rCommunicator, attributes, start_index, block_size)));
        KRATOS_ERROR_IF_NOT(is_component_written)
            << "Component \"" << r_component_name << "\" is not found in registered components.";
    }

    return attributes;

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
template <class TComponentType>
bool ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::WriteComponentData(
    const std::string& rComponentName,
    const TContainerDataIO& rContainerDataIO,
    const TContainerType& rLocalContainer,
    Parameters Attributes,
    WriteInfo& rInfo)
{
    KRATOS_TRY

    using component_data_type = typename Internals::template ComponentTraits<TComponentType>::ValueType;

    using value_type = typename TContainerDataIO::template ComponentDataType<component_data_type>;

    using value_type_traits = DataTypeTraits<value_type>;

    using value_primitive_type = typename value_type_traits::PrimitiveType;

    if (KratosComponents<TComponentType>::Has(rComponentName)) {
        const auto& r_data_set_path = mComponentPrefix + rComponentName;
        const auto& r_component = KratosComponents<TComponentType>::Get(rComponentName);

        std::vector<int> shape(value_type_traits::Dimension);
        if constexpr(value_type_traits::IsDynamic) {
            // retrieving dynamic shape
            if (!rLocalContainer.empty()) {
                typename TContainerDataIO::template TLSType<component_data_type> tls;
                const auto& value_prototype = rContainerDataIO.GetValue(rLocalContainer.front(), r_component, tls);
                value_type_traits::Shape(value_prototype, shape.data(), shape.data() + value_type_traits::Dimension);
            }
            shape = mpFile->GetDataCommunicator().MaxAll(shape);

            Matrix<value_primitive_type> values;
            values.resize(rLocalContainer.size(), std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<int>{}));
            Internals::CopyToContiguousDataArray(rLocalContainer, r_component, rContainerDataIO, DataTypeTraits<Matrix<value_primitive_type>>::GetContiguousData(values));
            mpFile->WriteDataSet(r_data_set_path, values, rInfo);
        } else {
            value_type_traits::Shape(value_type{}, shape.data(), shape.data() + value_type_traits::Dimension);
            Vector<value_type> values;
            values.resize(rLocalContainer.size());
            Internals::CopyToContiguousDataArray(rLocalContainer, r_component, rContainerDataIO, DataTypeTraits<Vector<value_type>>::GetContiguousData(values));
            mpFile->WriteDataSet(r_data_set_path, values, rInfo);
        }

        // add the shape to attributes.
        KRATOS_ERROR_IF(Attributes.Has("__data_shape"))
            << "The reserved keyword \"__data_shape\" is found. Please remove it from attributes.";
        if (shape.size() > 0) {
            Attributes.AddEmptyArray("__data_shape");
            for (const auto v : shape) {
                Attributes["__data_shape"].Append(v);
            }
        }

        KRATOS_ERROR_IF(Attributes.Has("__data_dimension"))
            << "The reserved keyword \"__data_dimension\" is found. Please remove it from attributes.";
        Attributes.AddInt("__data_dimension", value_type_traits::Dimension);

        KRATOS_ERROR_IF(Attributes.Has("__container_type"))
            << "The reserved keyword \"__container_type\" is found. Please remove it from attributes.";
        Attributes.AddString("__container_type", Internals::GetContainerType<TContainerType>());

        KRATOS_ERROR_IF(Attributes.Has("__data_location"))
            << "The reserved keyword \"__data_location\" is found. Please remove it from attributes.";
        Attributes.AddString("__data_location", Internals::GetContainerIOType<TContainerDataIO>());

        KRATOS_ERROR_IF(Attributes.Has("__data_name"))
            << "The reserved keyword \"__data_name\" is found. Please remove it from attributes.";
        Attributes.AddString("__data_name", rComponentName);

        // now write the attributes
        mpFile->WriteAttribute(r_data_set_path, Attributes);
        return true;
    } else {
        // if the component is not a match.
        return false;
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
template <class TComponentType>
bool ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::ReadComponentData(
    const std::string& rComponentName,
    const TContainerDataIO& rContainerDataIO,
    TContainerType& rLocalContainer,
    Communicator& rCommunicator,
    std::map<std::string, Parameters>& rAttributesMap,
    const IndexType StartIndex,
    const IndexType BlockSize)
{
    KRATOS_TRY

    using value_type = typename TContainerDataIO::template ComponentDataType<typename Internals::template ComponentTraits<TComponentType>::ValueType>;

    using value_type_traits = DataTypeTraits<value_type>;

    using value_primitive_type = typename value_type_traits::PrimitiveType;

    if (KratosComponents<TComponentType>::Has(rComponentName)) {
        const auto& r_data_set_path = mComponentPrefix + rComponentName;
        const auto& r_component = KratosComponents<TComponentType>::Get(rComponentName);

        KRATOS_ERROR_IF_NOT(mpFile->IsDataSet(r_data_set_path))
            << "Dataset \"" << r_data_set_path << "\" not found to read \""
            << rComponentName << "\".";

        KRATOS_ERROR_IF_NOT(BlockSize == rLocalContainer.size())
            << "Container size and the file block size mismatch [ rLocalContainer.size() = "
            << rLocalContainer.size() << ", file block size = " << BlockSize << " ].\n";

        KRATOS_ERROR_IF(rAttributesMap.find(rComponentName) != rAttributesMap.end())
            << "The attributes for component \"" << rComponentName << "\" already exists.";

        auto attributes = mpFile->ReadAttribute(r_data_set_path);

        if constexpr(value_type_traits::IsDynamic) {
            std::vector<unsigned int> shape;
            Internals::GetShapeFromAttributes(shape, attributes);

            Matrix<value_primitive_type> values;
            values.resize(BlockSize, std::accumulate(shape.begin(), shape.end(), 1U, std::multiplies<unsigned int>{}));
            mpFile->ReadDataSet(r_data_set_path, values, StartIndex, BlockSize);
            Internals::CopyFromContiguousDataArray(rLocalContainer, r_component, rContainerDataIO, DataTypeTraits<Matrix<value_primitive_type>>::GetContiguousData(values), shape);
        } else {
            Vector<value_type> values;
            values.resize(BlockSize);
            mpFile->ReadDataSet(r_data_set_path, values, StartIndex, BlockSize);
            Internals::CopyFromContiguousDataArray(rLocalContainer, r_component, rContainerDataIO, DataTypeTraits<Vector<value_type>>::GetContiguousData(values), value_type_traits::Shape(value_type{}));
        }

        NewContainerComponentIOUtilities::SynchronizeComponent<TContainerType>::template Execute<TContainerDataIO>(rCommunicator, r_component);

        if (attributes.Has("__data_dimension")) attributes.RemoveValue("__data_dimension");
        if (attributes.Has("__data_shape")) attributes.RemoveValue("__data_shape");
        if (attributes.Has("__container_type")) attributes.RemoveValue("__container_type");
        if (attributes.Has("__data_name")) attributes.RemoveValue("__data_name");
        if (attributes.Has("__data_location")) attributes.RemoveValue("__data_location");
        if (attributes.Has("__mesh_location")) attributes.RemoveValue("__mesh_location");

        rAttributesMap[rComponentName] = attributes;

        return true;
    } else {
        // if the component is not a match.
        return false;
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
std::map<std::string, Parameters> ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::ReadAttributes()
{
    KRATOS_TRY

    std::map<std::string, Parameters> attributes_map;

    for (const auto& r_component_name : mComponentNames) {
        KRATOS_ERROR_IF(attributes_map.find(r_component_name) != attributes_map.end())
            << "The attributes for component \"" << r_component_name << "\" already exists.";

        const auto& r_data_set_path = mComponentPrefix + r_component_name;
        auto attributes = mpFile->ReadAttribute(r_data_set_path);
        if (attributes.Has("__data_dimension")) attributes.RemoveValue("__data_dimension");
        if (attributes.Has("__data_shape")) attributes.RemoveValue("__data_shape");
        if (attributes.Has("__container_type")) attributes.RemoveValue("__container_type");
        if (attributes.Has("__data_name")) attributes.RemoveValue("__data_name");
        if (attributes.Has("__data_location")) attributes.RemoveValue("__data_location");
        if (attributes.Has("__mesh_location")) attributes.RemoveValue("__mesh_location");

        attributes_map[r_component_name] = attributes;
    }

    return attributes_map;

    KRATOS_CATCH("")
}


// template instantiations
#ifndef KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO
#define KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO(CONTAINER_TYPE, CONTAINER_DATA_IO, ...)                      \
template class KRATOS_API(HDF5_APPLICATION) ContainerComponentIO<CONTAINER_TYPE, CONTAINER_DATA_IO, __VA_ARGS__>;
#endif

#ifndef KRATOS_HDF5_INSTANTIATE_FLAGS_CONTAINER_COMPONENT_IO
#define KRATOS_HDF5_INSTANTIATE_FLAGS_CONTAINER_COMPONENT_IO(CONTAINER_TYPE)                \
KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO(CONTAINER_TYPE, Internals::FlagIO, Flags);
#endif

#ifndef KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO
#define KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(CONTAINER_TYPE, CONTAINER_DATA_IO)  \
KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO(CONTAINER_TYPE, CONTAINER_DATA_IO,                   \
                                               Variable<int>,                                       \
                                               Variable<double>,                                    \
                                               Variable<array_1d<double, 3>>,                       \
                                               Variable<array_1d<double, 4>>,                       \
                                               Variable<array_1d<double, 6>>,                       \
                                               Variable<array_1d<double, 9>>,                       \
                                               Variable<Kratos::Vector>,                            \
                                               Variable<Kratos::Matrix>);
#endif

#ifndef KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO
#define KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(CONTAINER_TYPE)                                              \
KRATOS_HDF5_INSTANTIATE_FLAGS_CONTAINER_COMPONENT_IO(CONTAINER_TYPE)                                                        \
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(CONTAINER_TYPE, Internals::NonHistoricalIO)
#endif

KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(ModelPart::NodesContainerType, Internals::HistoricalIO);
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(ModelPart::NodesContainerType, Internals::BossakIO);
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(Detail::VertexContainerType, Internals::VertexValueIO);
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(ModelPart::ConditionsContainerType, Internals::GaussPointValueIO);
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(ModelPart::ElementsContainerType, Internals::GaussPointValueIO);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::NodesContainerType);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::ConditionsContainerType);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::ElementsContainerType);


#undef KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_FLAGS_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO

} // namespace HDF5.
} // namespace Kratos.
