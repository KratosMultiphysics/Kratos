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
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "custom_utilities/data_type_utilities.h"

// Include base h
#include "custom_io/hdf5_new_container_component_io.h"

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
        } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>) {
            rCommunicator.SynchronizeVariable(rComponent);
        } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
            rCommunicator.SynchronizeNonHistoricalVariable(rComponent);
        } else {
            static_assert(!std::is_same_v<TContainerDataIO, TContainerDataIO>, "Unsupported TContainerDataIO.");
        }
    }
};
} // namespace NewContainerComponentIOUtilities

template <class TContainerType, class TContainerDataIO, class... TComponents>
NewContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::NewContainerComponentIO(
    Parameters Settings,
    File::Pointer pFile)
    : mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix"           : "",
            "list_of_variables": []
        })");

    Settings.ValidateAndAssignDefaults(default_params);

    mComponentPath = Settings["prefix"].GetString();
    mComponentNames = Settings["list_of_variables"].GetStringArray();

    // Sort component names to make sure they're in the same order on each rank.
    // The basic assumption is that the set of components is identical on every
    // rank, but they may not be in the same order (which would lead to ranks
    // trying to write different variables at the same time, resulting in
    // a deadlock), hence the sorting.
    std::sort(mComponentNames.begin(), mComponentNames.end());

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
void NewContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::Write(
    TContainerType const& rLocalContainer,
    Parameters Attributes)
{
    KRATOS_TRY;

    if (mComponentNames.size() == 0) {
        return;
    }

    WriteInfo info;

    // Write each variable.
    for (const auto& r_component_name : mComponentNames) {
        const bool is_component_written = (... || (WriteComponentData<TComponents>(
                        r_component_name, rLocalContainer, Attributes.Clone(), info)));
        KRATOS_ERROR_IF_NOT(is_component_written)
            << "Component \"" << r_component_name << "\" is not found in registered components.";
    }

    // Write block partition.
    WritePartitionTable(*mpFile, mComponentPath, info);

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
std::map<std::string, Parameters> NewContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::Read(
    TContainerType& rLocalContainer,
    Communicator& rCommunicator)
{
    KRATOS_TRY;

    std::map<std::string, Parameters> attributes;

    if (mComponentNames.size() == 0) {
        return attributes;
    }

    IndexType start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, mComponentPath);

    // Write each variable.
    for (const auto& r_component_name : mComponentNames) {
        const bool is_component_written =
            (... || (ReadComponentData<TComponents>(
                        r_component_name, rLocalContainer, rCommunicator, attributes, start_index, block_size)));
        KRATOS_ERROR_IF_NOT(is_component_written)
            << "Component \"" << r_component_name << "\" is not found in registered components.";
    }

    return attributes;

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
template <class TComponentType>
bool NewContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::WriteComponentData(
    const std::string& rComponentName,
    const TContainerType& rLocalContainer,
    Parameters Attributes,
    WriteInfo& rInfo)
{
    using value_type = typename Internals::ComponentTraits<TComponentType>::ValueType;

    using value_type_traits = DataTypeTraits<value_type>;

    using value_primitive_type = typename value_type_traits::PrimitiveType;

    if (KratosComponents<TComponentType>::Has(rComponentName)) {
        const auto& r_data_set_path = mComponentPath + "/" + rComponentName;
        const auto& r_component = KratosComponents<TComponentType>::Get(rComponentName);

        if constexpr(value_type_traits::IsDynamic) {
            // dynamic shape writing
            value_type value_prototype;
            if (!rLocalContainer.empty()) {
                value_prototype = TContainerDataIO::GetValue(rLocalContainer.front(), r_component);
            }
            mpFile->GetDataCommunicator().SynchronizeShape(value_prototype);

            const auto& shape = value_type_traits::template Shape<int>(value_prototype);

            // add the shape to attributes.
            KRATOS_ERROR_IF(Attributes.Has("__shape"))
                << "The reserved keyword \"__shape\" is found. Please remove it from attributes.";
            Attributes.AddEmptyArray("__shape");
            for (const auto v : shape) {
                Attributes["__shape"].Append(v);
            }

            Matrix<value_primitive_type> values;
            values.resize(rLocalContainer.size(), value_type_traits::Size(value_prototype));
            Internals::CopyToContiguousDataArray<TContainerDataIO>(rLocalContainer, r_component, DataTypeTraits<Matrix<value_primitive_type>>::GetContiguousData(values));
            mpFile->WriteDataSet(r_data_set_path, values, rInfo);
        } else {
            Vector<value_type> values;
            values.resize(rLocalContainer.size());
            Internals::CopyToContiguousDataArray<TContainerDataIO>(rLocalContainer, r_component, DataTypeTraits<Vector<value_type>>::GetContiguousData(values));
            mpFile->WriteDataSet(r_data_set_path, values, rInfo);
        }

        KRATOS_ERROR_IF(Attributes.Has("__container_type"))
            << "The reserved keyword \"__container_type\" is found. Please remove it from attributes.";
        Attributes.AddString("__container_type", Internals::GetContainerType<TContainerType>());

        KRATOS_ERROR_IF(Attributes.Has("__component_name"))
            << "The reserved keyword \"__component_name\" is found. Please remove it from attributes.";
        Attributes.AddString("__component_name", rComponentName);

        KRATOS_ERROR_IF(Attributes.Has("__data_location"))
            << "The reserved keyword \"__data_location\" is found. Please remove it from attributes.";
        Attributes.AddString("__data_location", Internals::GetContainerIOType<TContainerDataIO>());

        // now write the attributes
        mpFile->WriteAttribute(r_data_set_path, Attributes);
        return true;
    } else {
        // if the component is not a match.
        return false;
    }
}

template <class TContainerType, class TContainerDataIO, class... TComponents>
template <class TComponentType>
bool NewContainerComponentIO<TContainerType, TContainerDataIO, TComponents...>::ReadComponentData(
    const std::string& rComponentName,
    TContainerType& rLocalContainer,
    Communicator& rCommunicator,
    std::map<std::string, Parameters>& rAttributesMap,
    const IndexType StartIndex,
    const IndexType BlockSize)
{
    using value_type = typename Internals::ComponentTraits<TComponentType>::ValueType;

    using value_type_traits = DataTypeTraits<value_type>;

    using value_primitive_type = typename value_type_traits::PrimitiveType;

    if (KratosComponents<TComponentType>::Has(rComponentName)) {
        const auto& r_data_set_path = mComponentPath + "/" + rComponentName;
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
            Internals::CopyFromContiguousDataArray<TContainerDataIO>(rLocalContainer, r_component, DataTypeTraits<Matrix<value_primitive_type>>::GetContiguousData(values), shape);
        } else {
            Vector<value_type> values;
            values.resize(BlockSize);
            mpFile->ReadDataSet(r_data_set_path, values, StartIndex, BlockSize);
            Internals::CopyFromContiguousDataArray<TContainerDataIO>(rLocalContainer, r_component, DataTypeTraits<Vector<value_type>>::GetContiguousData(values), value_type_traits::Shape(value_type{}));
        }

        NewContainerComponentIOUtilities::SynchronizeComponent<TContainerType>::template Execute<TContainerDataIO>(rCommunicator, r_component);

        if (attributes.Has("__shape")) {
            attributes.RemoveValue("__shape");
        }

        rAttributesMap[rComponentName] = attributes;

        return true;
    } else {
        // if the component is not a match.
        return false;
    }
}


// template instantiations
#ifndef KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO
#define KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO(CONTAINER_TYPE, CONTAINER_DATA_IO, ...)                      \
template class KRATOS_API(HDF5_APPLICATION) NewContainerComponentIO<CONTAINER_TYPE, CONTAINER_DATA_IO, __VA_ARGS__>;
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
KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(CONTAINER_TYPE, ContainerDataIO<ContainerDataIOTags::NonHistorical>)
#endif

KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO(ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::NodesContainerType);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::ConditionsContainerType);
KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO(ModelPart::ElementsContainerType);


#undef KRATOS_HDF5_INSTANTIATE_GENERIC_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_VARIABLE_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_FLAGS_CONTAINER_COMPONENT_IO
#undef KRATOS_HDF5_INSTANTIATE_CONTAINER_COMPONENT_IO

} // namespace HDF5.
} // namespace Kratos.
