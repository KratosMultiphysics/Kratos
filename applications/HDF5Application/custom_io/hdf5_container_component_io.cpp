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
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/container_io_utils.h"
#include "custom_utilities/data_type_utilities.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"

// Include base h
#include "custom_io/hdf5_container_component_io.h"

namespace Kratos
{
namespace HDF5
{

// TODO: rename this namespace once HDF5 refactoring is complete
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
void ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::CheckReservedAttributes(const Parameters Attributes)
{
    KRATOS_TRY

    for (const auto& r_attribute_key : ReservedAttributeKeys) {
        KRATOS_ERROR_IF(Attributes.Has(r_attribute_key))
            << "The reserved keyword \"" << r_attribute_key << "\" is found. Please remove it from attributes."
            << "Followings are the given attributes:\n" << Attributes << std::endl;
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
void ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::RemoveReservedAttributes(Parameters Attributes)
{
    KRATOS_TRY

    for (const auto& r_attribute_key : ReservedAttributeKeys) {
        if (Attributes.Has(r_attribute_key)) Attributes.RemoveValue(r_attribute_key);
    }

    KRATOS_CATCH("");
}

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
        << "The prefix must not be empty or \"/\" [ prefix = " << mComponentPrefix << " ].\n";

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

    if constexpr(
        std::is_same_v<TContainerType, NodesContainerType> ||
        std::is_same_v<TContainerType, ConditionsContainerType> ||
        std::is_same_v<TContainerType, ElementsContainerType>) {
        Write(Internals::GetLocalContainer<TContainerType>(rModelPart), rContainerDataIO, Attributes);
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
        return {};
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
        const bool is_component_written = (... || (WriteComponents<TComponents>(
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
            (... || (ReadComponents<TComponents>(
                        r_component_name, rContainerDataIO, rLocalContainer, rCommunicator, attributes, start_index, block_size)));
        KRATOS_ERROR_IF_NOT(is_component_written)
            << "Component \"" << r_component_name << "\" is not found in registered components.";
    }

    return attributes;

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
template <class TComponentType>
bool ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::WriteComponents(
    const std::string& rComponentName,
    const TContainerDataIO& rContainerDataIO,
    const TContainerType& rLocalContainer,
    Parameters Attributes,
    WriteInfo& rInfo)
{
    KRATOS_TRY

    using component_type = typename Internals::template ComponentTraits<TComponentType>::ValueType;

    using value_type = typename TContainerDataIO::template ComponentType<component_type>;

    using value_type_traits = DataTypeTraits<value_type>;

    using value_primitive_type = typename value_type_traits::PrimitiveType;

    if (KratosComponents<TComponentType>::Has(rComponentName)) {
        const auto& r_data_set_path = mComponentPrefix + rComponentName;
        const auto& r_component = KratosComponents<TComponentType>::Get(rComponentName);

        std::vector<int> shape(value_type_traits::Dimension);
        if constexpr(value_type_traits::Dimension == 0) {
            Vector<value_type> values(rLocalContainer.size());
            Internals::CopyToContiguousArray(rLocalContainer, r_component, rContainerDataIO, DataTypeTraits<Vector<value_type>>::GetContiguousData(values), values.size());
            mpFile->WriteDataSet(r_data_set_path, values, rInfo);
        } else {
            // get the correct size
            value_type value_prototype{};
            if (!rLocalContainer.empty()) {
                // if the value type is not static, then we need to get the shape
                // of the first element assuming all the entities will have the same
                // shape.
                typename TContainerDataIO::template TLSType<component_type> tls;
                value_prototype = rContainerDataIO.GetValue(rLocalContainer.front(), r_component, tls);
            }

            // now retrieve the shape
            value_type_traits::Shape(value_prototype, shape.data(), shape.data() + value_type_traits::Dimension);

            // communicate between ranks(to provide shape for the ranks which are empty.)
            shape = mpFile->GetDataCommunicator().MaxAll(shape);

            Matrix<value_primitive_type> values;
            values.resize(rLocalContainer.size(), std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<int>{}));
            Internals::CopyToContiguousArray(rLocalContainer, r_component, rContainerDataIO, DataTypeTraits<Matrix<value_primitive_type>>::GetContiguousData(values), values.size1() * values.size2());
            mpFile->WriteDataSet(r_data_set_path, values, rInfo);
        }

        // now writing the availability of each value
        // it is only required to calculate if the components are optional
        if constexpr(TContainerDataIO::RequireAvailabilityCheck) {
            Vector<bool> availability(rLocalContainer.size());

            const auto availability_local_counts_pair = IndexPartition<IndexType>(rLocalContainer.size()).for_each<CombinedReduction<SumReduction<IndexType>, SumReduction<IndexType>>>([&rContainerDataIO, &rLocalContainer, &r_component, &availability](const auto Index) {
                availability[Index] = rContainerDataIO.HasValue(*(rLocalContainer.begin() + Index), r_component);
                return std::make_tuple<IndexType, IndexType>(availability[Index] == true, availability[Index] == false);
            });

            const auto& availability_global_counts_pair = mpFile->GetDataCommunicator().SumAll(std::vector<IndexType>{std::get<0>(availability_local_counts_pair), std::get<1>(availability_local_counts_pair), rLocalContainer.size()});

            if (availability_global_counts_pair[0] != availability_global_counts_pair[2] && availability_global_counts_pair[1] != availability_global_counts_pair[2]) {
                // number of entities having the component defined is not equal to total number of entities or, number of entities
                // not having the component defined is not equal to total number of entities.
                // therefore, it is essential to write the availability to record which entity has the component.
                mpFile->WriteDataSet(r_data_set_path + "_availability", availability, rInfo);
            }
        }

        // add the shape to attributes.
        CheckReservedAttributes(Attributes);

        if (shape.size() > 0) {
            Attributes.AddEmptyArray("__data_shape");
            for (const auto v : shape) {
                Attributes["__data_shape"].Append(v);
            }
        }
        Attributes.AddInt("__data_dimension", value_type_traits::Dimension);
        Attributes.AddString("__container_type", Internals::GetContainerName<TContainerType>());
        Attributes.AddString("__data_location", Internals::GetContainerIOName<TContainerDataIO>());
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
bool ContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::ReadComponents(
    const std::string& rComponentName,
    const TContainerDataIO& rContainerDataIO,
    TContainerType& rLocalContainer,
    Communicator& rCommunicator,
    std::map<std::string, Parameters>& rAttributesMap,
    const IndexType StartIndex,
    const IndexType BlockSize)
{
    KRATOS_TRY

    using value_type = typename TContainerDataIO::template ComponentType<typename Internals::template ComponentTraits<TComponentType>::ValueType>;

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

        RemoveReservedAttributes(attributes);
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
        RemoveReservedAttributes(attributes);
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
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(Detail::VertexContainerType, Internals::VertexHistoricalValueIO);
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(Detail::VertexContainerType, Internals::VertexNonHistoricalValueIO);
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(ModelPart::ConditionsContainerType, Internals::GaussPointIO);
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(ModelPart::ElementsContainerType, Internals::GaussPointIO);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::NodesContainerType);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::ConditionsContainerType);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::ElementsContainerType);


#undef KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_FLAGS_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO

} // namespace HDF5.
} // namespace Kratos.
