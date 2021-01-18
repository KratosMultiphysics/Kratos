#include "custom_io/hdf5_container_component_io.h"

#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "custom_utilities/local_ghost_splitting_utility.h"
#include "custom_utilities/registered_component_lookup.h"
#include "includes/communicator.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
namespace HDF5
{
namespace
{
template <typename TContainerItemType, typename TDataType>
void SetDataBuffer(Variable<TDataType> const&,
                   std::vector<TContainerItemType*> const&,
                   Vector<TDataType>&);

template <typename TContainerItemType>
void SetDataBuffer(Flags const&, std::vector<TContainerItemType*> const&, Vector<int>&);

template <typename TContainerItemType>
void SetDataBuffer(Variable<Vector<double>> const&,
                   std::vector<TContainerItemType*> const&,
                   Matrix<double>&);

template <typename TContainerItemType>
void SetDataBuffer(Variable<Matrix<double>> const&,
                   std::vector<TContainerItemType*> const&,
                   Matrix<double>&);

template <typename TContainerItemType, typename TDataType>
void SetItemDataValues(Variable<TDataType> const&,
                       Vector<TDataType> const&,
                       std::vector<TContainerItemType*>&);

template <typename TContainerItemType>
void SetItemDataValues(Flags const&, Vector<int> const&, std::vector<TContainerItemType*>&);

template <typename TContainerItemType>
void SetItemDataValues(Variable<Vector<double>> const&,
                       Matrix<double> const&,
                       std::vector<TContainerItemType*>&);

template <typename TContainerItemType>
void SetItemDataValues(Variable<Vector<double>> const&,
                       Matrix<double> const&,
                       std::vector<TContainerItemType*>&);

template <typename TContainerItemType>
void SetItemDataValues(Variable<Matrix<double>> const&,
                       Matrix<double> const&,
                       std::vector<TContainerItemType*>&,
                       int,
                       int);

template <typename TContainerItemType>
class ContainerItemComponentIO
{
public:
    template <typename TComponent, typename dummy = void>
    class WriteComponentFunctor
    {
    public:
        void operator()(TComponent const& rComponent,
                        std::vector<TContainerItemType*>& rContainerItems,
                        File& rFile,
                        std::string const& rPath,
                        WriteInfo& rInfo,
                        std::string const&)
        {
            Vector<typename TComponent::Type> data;
            SetDataBuffer(rComponent, rContainerItems, data);
            rFile.WriteDataSet(rPath + "/" + rComponent.Name(), data, rInfo);
        }
    };

    template <typename dummy>
    class WriteComponentFunctor<Flags, dummy>
    {
    public:
        void operator()(Flags const& rComponent,
                        std::vector<TContainerItemType*>& rContainerItems,
                        File& rFile,
                        std::string const& rPath,
                        WriteInfo& rInfo,
                        std::string const& rComponentName)
        {
            Vector<int> data;
            SetDataBuffer(rComponent, rContainerItems, data);
            rFile.WriteDataSet(rPath + "/" + rComponentName, data, rInfo);
        }
    };

    template <typename dummy>
    class WriteComponentFunctor<Variable<Vector<double>>, dummy>
    {
    public:
        void operator()(Variable<Vector<double>> const& rComponent,
                        std::vector<TContainerItemType*>& rContainerItems,
                        File& rFile,
                        std::string const& rPath,
                        WriteInfo& rInfo,
                        std::string const&)
        {
            Matrix<double> data;
            SetDataBuffer(rComponent, rContainerItems, data);
            rFile.WriteDataSet(rPath + "/" + rComponent.Name(), data, rInfo);
        }
    };

    template <typename dummy>
    class WriteComponentFunctor<Variable<Matrix<double>>, dummy>
    {
    public:
        void operator()(Variable<Matrix<double>> const& rComponent,
                        std::vector<TContainerItemType*>& rContainerItems,
                        File& rFile,
                        std::string const& rPath,
                        WriteInfo& rInfo,
                        std::string const&)
        {
            Matrix<double> data;
            SetDataBuffer(rComponent, rContainerItems, data);
            rFile.WriteDataSet(rPath + "/" + rComponent.Name(), data, rInfo);
            const int size1 = rContainerItems.front()->GetValue(rComponent).size1();
            const int size2 = rContainerItems.front()->GetValue(rComponent).size2();
            rFile.WriteAttribute(rPath + "/" + rComponent.Name(), "Size1", size1);
            rFile.WriteAttribute(rPath + "/" + rComponent.Name(), "Size2", size2);
        }
    };

    template <typename TComponent, template <typename T> class TSynchronizingFunctor>
    class ReadComponentFunctor
    {
    public:
        void operator()(TComponent const& rComponent,
                        std::vector<TContainerItemType*>& rLocalItems,
                        Communicator& rCommunicator,
                        File& rFile,
                        std::string const& rPath,
                        std::size_t StartIndex,
                        std::size_t BlockSize,
                        std::string const&)
        {
            Vector<typename TComponent::Type> data;
            rFile.ReadDataSet(rPath + "/" + rComponent.Name(), data, StartIndex, BlockSize);
            SetItemDataValues(rComponent, data, rLocalItems);
            TSynchronizingFunctor<TComponent>()(rComponent, rCommunicator);
        }
    };

    template <template <typename T> class TSynchronizingFunctor>
    class ReadComponentFunctor<Flags, TSynchronizingFunctor>
    {
    public:
        void operator()(Flags const& rComponent,
                        std::vector<TContainerItemType*>& rLocalItems,
                        Communicator& rCommunicator,
                        File& rFile,
                        std::string const& rPath,
                        std::size_t StartIndex,
                        std::size_t BlockSize,
                        std::string const& rComponentName)
        {
            Vector<int> data;
            rFile.ReadDataSet(rPath + "/" + rComponentName, data, StartIndex, BlockSize);
            SetItemDataValues(rComponent, data, rLocalItems);
        }
    };

    template <template <typename T> class TSynchronizingFunctor>
    class ReadComponentFunctor<Variable<Vector<double>>, TSynchronizingFunctor>
    {
    public:
        void operator()(Variable<Vector<double>> const& rComponent,
                        std::vector<TContainerItemType*>& rLocalItems,
                        Communicator& rCommunicator,
                        File& rFile,
                        std::string const& rPath,
                        std::size_t StartIndex,
                        std::size_t BlockSize,
                        std::string const&)
        {
            Matrix<double> data;
            rFile.ReadDataSet(rPath + "/" + rComponent.Name(), data, StartIndex, BlockSize);
            SetItemDataValues(rComponent, data, rLocalItems);
            TSynchronizingFunctor<Variable<Vector<double>>>()(rComponent, rCommunicator);
        }
    };

    template <template <typename T> class TSynchronizingFunctor>
    class ReadComponentFunctor<Variable<Matrix<double>>, TSynchronizingFunctor>
    {
    public:
        void operator()(Variable<Matrix<double>> const& rComponent,
                        std::vector<TContainerItemType*>& rLocalItems,
                        Communicator& rCommunicator,
                        File& rFile,
                        std::string const& rPath,
                        std::size_t StartIndex,
                        std::size_t BlockSize,
                        std::string const&)
        {
            Matrix<double> data;
            rFile.ReadDataSet(rPath + "/" + rComponent.Name(), data, StartIndex, BlockSize);
            int size1, size2;
            rFile.ReadAttribute(rPath + "/" + rComponent.Name(), "Size1", size1);
            rFile.ReadAttribute(rPath + "/" + rComponent.Name(), "Size2", size2);
            SetItemDataValues(rComponent, data, rLocalItems, size1, size2);
            TSynchronizingFunctor<Variable<Matrix<double>>>()(rComponent, rCommunicator);
        }
    };
};

template <typename TComponent>
class ComponentElementSynchronizingFunctor
{
public:
    void operator()(TComponent const& rComponent, Communicator& rCommunicator)
    {
    }
};

template <typename TComponent>
class ComponentConditionSynchronizingFunctor
{
public:
    void operator()(TComponent const& rComponent, Communicator& rCommunicator)
    {
    }
};

template <typename TComponent>
class ComponentNodeSynchronizingFunctor
{
public:
    void operator()(TComponent const& rComponent, Communicator& rCommunicator)
    {
        rCommunicator.SynchronizeNonHistoricalVariable(rComponent);
    }
};

template <>
class ComponentNodeSynchronizingFunctor<Flags>
{
public:
    void operator()(Flags const& rComponent, std::vector<NodeType*>& rGhostNodes, Communicator& rCommunicator)
    {
    }
};

template <typename TComponent>
class WriteElementComponent
    : public ContainerItemComponentIO<ElementType>::WriteComponentFunctor<TComponent, void>
{
};

template <typename TComponent>
class ReadElementComponent
    : public ContainerItemComponentIO<ElementType>::ReadComponentFunctor<TComponent, ComponentElementSynchronizingFunctor>
{
};

template <typename TComponent>
class WriteConditionComponent
    : public ContainerItemComponentIO<ConditionType>::WriteComponentFunctor<TComponent, void>
{
};

template <typename TComponent>
class ReadConditionComponent
    : public ContainerItemComponentIO<ConditionType>::ReadComponentFunctor<TComponent, ComponentConditionSynchronizingFunctor>
{
};

template <typename TComponent>
class WriteNodeComponent
    : public ContainerItemComponentIO<NodeType>::WriteComponentFunctor<TComponent, void>
{
};

template <typename TComponent>
class ReadNodeComponent
    : public ContainerItemComponentIO<NodeType>::ReadComponentFunctor<TComponent, ComponentNodeSynchronizingFunctor>
{
};

void GetContainerItemReferences(std::vector<ElementType*>& rLocalElements,
                                ElementsContainerType const& rElements)
{
    rLocalElements.resize(rElements.size());

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
    {
        auto it = rElements.begin() + i;
        rLocalElements[i] = (&(*it));
    }
}

void GetContainerItemReferences(std::vector<ConditionType*>& rLocalConditions,
                                ConditionsContainerType const& rConditions)
{
    rLocalConditions.resize(rConditions.size());

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rConditions.size()); ++i)
    {
        auto it = rConditions.begin() + i;
        rLocalConditions[i] = (&(*it));
    }
}

void GetContainerItemReferences(std::vector<NodeType*>& rLocalNodes,
                                NodesContainerType const& rNodes)
{
    GetLocalNodes(rNodes, rLocalNodes);
}

} // unnamed namespace

template <typename TContainerType, typename TContainerItemType, typename... TComponents>
ContainerComponentIO<TContainerType, TContainerItemType, TComponents...>::ContainerComponentIO(
    Parameters Settings, File::Pointer pFile, const std::string& rComponentPath)
    : mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix": "",
            "list_of_variables": []
        })");

    Settings.ValidateAndAssignDefaults(default_params);

    mComponentPath = Settings["prefix"].GetString() + rComponentPath;

    const std::size_t num_components = Settings["list_of_variables"].size();

    if (mComponentNames.size() != num_components)
        mComponentNames.resize(num_components);

    for (std::size_t i = 0; i < num_components; ++i)
        mComponentNames[i] = Settings["list_of_variables"].GetArrayItem(i).GetString();

    KRATOS_CATCH("");
}

// Adding element flags WriteRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<ElementsContainerType, ElementType, Flags>::WriteRegisteredComponent(
    const std::string& rComponentName, Targs&... args)
{
    RegisteredComponentLookup<Flags>(rComponentName).Execute<WriteElementComponent>(args...);
}

// Adding element other variable WriteRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<ElementsContainerType,
                          ElementType,
                          Variable<array_1d<double, 3>>,
                          Variable<double>,
                          Variable<int>,
                          Variable<Vector<double>>,
                          Variable<Matrix<double>>>::WriteRegisteredComponent(const std::string& rComponentName,
                                                                              Targs&... args)
{
    RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                              Variable<Vector<double>>, Variable<Matrix<double>>>(rComponentName)
        .Execute<WriteElementComponent>(args...);
}

// Adding element flags ReadRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<ElementsContainerType, ElementType, Flags>::ReadRegisteredComponent(
    const std::string& rComponentName, Targs&... args)
{
    RegisteredComponentLookup<Flags>(rComponentName).Execute<ReadElementComponent>(args...);
}

// Adding element other variable ReadRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<ElementsContainerType,
                          ElementType,
                          Variable<array_1d<double, 3>>,
                          Variable<double>,
                          Variable<int>,
                          Variable<Vector<double>>,
                          Variable<Matrix<double>>>::ReadRegisteredComponent(const std::string& rComponentName,
                                                                             Targs&... args)
{
    RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                              Variable<Vector<double>>, Variable<Matrix<double>>>(rComponentName)
        .Execute<ReadElementComponent>(args...);
}

// Adding condition flags WriteRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<ConditionsContainerType, ConditionType, Flags>::WriteRegisteredComponent(
    const std::string& rComponentName, Targs&... args)
{
    RegisteredComponentLookup<Flags>(rComponentName).Execute<WriteConditionComponent>(args...);
}

// Adding condition other variable WriteRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<ConditionsContainerType,
                          ConditionType,
                          Variable<array_1d<double, 3>>,
                          Variable<double>,
                          Variable<int>,
                          Variable<Vector<double>>,
                          Variable<Matrix<double>>>::WriteRegisteredComponent(const std::string& rComponentName,
                                                                              Targs&... args)
{
    RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                              Variable<Vector<double>>, Variable<Matrix<double>>>(rComponentName)
        .Execute<WriteConditionComponent>(args...);
}

// Adding condition flags ReadRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<ConditionsContainerType, ConditionType, Flags>::ReadRegisteredComponent(
    const std::string& rComponentName, Targs&... args)
{
    RegisteredComponentLookup<Flags>(rComponentName).Execute<ReadConditionComponent>(args...);
}

// Adding condition other variable ReadRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<ConditionsContainerType,
                          ConditionType,
                          Variable<array_1d<double, 3>>,
                          Variable<double>,
                          Variable<int>,
                          Variable<Vector<double>>,
                          Variable<Matrix<double>>>::ReadRegisteredComponent(const std::string& rComponentName,
                                                                             Targs&... args)
{
    RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                              Variable<Vector<double>>, Variable<Matrix<double>>>(rComponentName)
        .Execute<ReadConditionComponent>(args...);
}

// Adding node flags WriteRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<NodesContainerType, NodeType, Flags>::WriteRegisteredComponent(
    const std::string& rComponentName, Targs&... args)
{
    RegisteredComponentLookup<Flags>(rComponentName).Execute<WriteNodeComponent>(args...);
}

// Adding node other variable WriteRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<NodesContainerType,
                          NodeType,
                          Variable<array_1d<double, 3>>,
                          Variable<double>,
                          Variable<int>,
                          Variable<Vector<double>>,
                          Variable<Matrix<double>>>::WriteRegisteredComponent(const std::string& rComponentName,
                                                                              Targs&... args)
{
    RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                              Variable<Vector<double>>, Variable<Matrix<double>>>(rComponentName)
        .Execute<WriteNodeComponent>(args...);
}

// Adding node flags ReadRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<NodesContainerType, NodeType, Flags>::ReadRegisteredComponent(
    const std::string& rComponentName, Targs&... args)
{
    RegisteredComponentLookup<Flags>(rComponentName).Execute<ReadNodeComponent>(args...);
}

// Adding node other variable ReadRegisteredComponent template definition
template <>
template <typename... Targs>
void ContainerComponentIO<NodesContainerType,
                          NodeType,
                          Variable<array_1d<double, 3>>,
                          Variable<double>,
                          Variable<int>,
                          Variable<Vector<double>>,
                          Variable<Matrix<double>>>::ReadRegisteredComponent(const std::string& rComponentName,
                                                                             Targs&... args)
{
    RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                              Variable<Vector<double>>, Variable<Matrix<double>>>(rComponentName)
        .Execute<ReadNodeComponent>(args...);
}

template <typename TContainerType, typename TContainerItemType, typename... TComponents>
void ContainerComponentIO<TContainerType, TContainerItemType, TComponents...>::WriteContainerComponents(
    TContainerType const& rContainer)
{
    KRATOS_TRY;

    if (mComponentNames.size() == 0)
        return;

    WriteInfo info;
    std::vector<TContainerItemType*> local_items;
    GetContainerItemReferences(local_items, rContainer);

    // Write each variable.
    for (const std::string& r_component_name : mComponentNames)
        WriteRegisteredComponent(r_component_name, local_items, *mpFile,
                                 mComponentPath, info, r_component_name);

    // Write block partition.
    WritePartitionTable(*mpFile, mComponentPath, info);

    KRATOS_CATCH("");
}

template <typename TContainerType, typename TContainerItemType, typename... TComponents>
void ContainerComponentIO<TContainerType, TContainerItemType, TComponents...>::ReadContainerComponents(
    TContainerType& rContainer, Communicator& rCommunicator)
{
    KRATOS_TRY;

    if (mComponentNames.size() == 0)
        return;

    std::vector<TContainerItemType*> local_items;
    GetContainerItemReferences(local_items, rContainer);
    std::size_t start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, mComponentPath);

    // Read local data for each variable.
    for (const std::string& r_component_name : mComponentNames)
        ReadRegisteredComponent(r_component_name, local_items, rCommunicator, *mpFile,
                                mComponentPath, start_index, block_size, r_component_name);

    KRATOS_CATCH("");
}

namespace
{
template <typename TContainerItemType, typename TDataType>
void SetDataBuffer(Variable<TDataType> const& rVariable,
                   std::vector<TContainerItemType*> const& rContainerItems,
                   Vector<TDataType>& rData)
{
    KRATOS_TRY;

    rData.resize(rContainerItems.size(), false);
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rContainerItems.size()); ++i)
    {
        const TContainerItemType& r_item = *rContainerItems[i];
        rData[i] = r_item.GetValue(rVariable);
    }

    KRATOS_CATCH("");
}

template <typename TContainerItemType>
void SetDataBuffer(Flags const& rVariable,
                   std::vector<TContainerItemType*> const& rContainerItems,
                   Vector<int>& rData)
{
    KRATOS_TRY;

    rData.resize(rContainerItems.size(), false);
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rContainerItems.size()); ++i)
    {
        const TContainerItemType& r_item = *rContainerItems[i];
        rData[i] = static_cast<int>(r_item.Is(rVariable));
    }

    KRATOS_CATCH("");
}

template <typename TContainerItemType>
void SetDataBuffer(Variable<Vector<double>> const& rVariable,
                   std::vector<TContainerItemType*> const& rContainerItems,
                   Matrix<double>& rData)
{
    KRATOS_TRY;

    rData.resize(rContainerItems.size(),
                 rContainerItems.front()->GetValue(rVariable).size(), false);
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rContainerItems.size()); ++i)
    {
        const TContainerItemType& r_item = *rContainerItems[i];
        Vector<double> const& r_vec = r_item.GetValue(rVariable);
        KRATOS_ERROR_IF(r_vec.size() != rData.size2())
            << "Invalid dimension.\n";
        for (std::size_t j = 0; j < rData.size2(); ++j)
            rData(i, j) = r_vec(j);
    }

    KRATOS_CATCH("");
}

template <typename TContainerItemType>
void SetDataBuffer(Variable<Matrix<double>> const& rVariable,
                   std::vector<TContainerItemType*> const& rContainerItems,
                   Matrix<double>& rData)
{
    KRATOS_TRY;

    rData.resize(rContainerItems.size(),
                 rContainerItems.front()->GetValue(rVariable).data().size(), false);
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rContainerItems.size()); ++i)
    {
        const TContainerItemType& r_item = *rContainerItems[i];
        const auto& r_data = r_item.GetValue(rVariable).data();
        KRATOS_ERROR_IF(r_data.size() != rData.size2())
            << "Invalid dimension.\n";
        for (std::size_t j = 0; j < rData.size2(); ++j)
            rData(i, j) = r_data[j];
    }

    KRATOS_CATCH("");
}

template <typename TContainerItemType, typename TDataType>
void SetItemDataValues(Variable<TDataType> const& rVariable,
                       Vector<TDataType> const& rData,
                       std::vector<TContainerItemType*>& rContainerItems)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rContainerItems.size() != rData.size())
        << "Number of items does not equal data set dimension\n";
    for (std::size_t i = 0; i < rContainerItems.size(); ++i)
        rContainerItems[i]->GetValue(rVariable) = rData[i];

    KRATOS_CATCH("");
}

template <typename TContainerItemType>
void SetItemDataValues(Flags const& rVariable,
                       Vector<int> const& rData,
                       std::vector<TContainerItemType*>& rContainerItems)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rContainerItems.size() != rData.size())
        << "Number of items does not equal data set dimension\n";
    for (std::size_t i = 0; i < rContainerItems.size(); ++i)
        rContainerItems[i]->Set(rVariable, static_cast<bool>(rData[i]));

    KRATOS_CATCH("");
}

template <typename TContainerItemType>
void SetItemDataValues(Variable<Vector<double>> const& rVariable,
                       Matrix<double> const& rData,
                       std::vector<TContainerItemType*>& rContainerItems)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rContainerItems.size() != rData.size1())
        << "Number of items does not equal data set dimension\n";
    for (std::size_t i = 0; i < rContainerItems.size(); ++i)
    {
        Vector<double>& r_vec = rContainerItems[i]->GetValue(rVariable);
        r_vec.resize(rData.size2(), false);
        for (std::size_t j = 0; j < rData.size2(); ++j)
            r_vec(j) = rData(i, j);
    }

    KRATOS_CATCH("");
}

template <typename TContainerItemType>
void SetItemDataValues(Variable<Matrix<double>> const& rVariable,
                       Matrix<double> const& rData,
                       std::vector<TContainerItemType*>& rContainerItems,
                       int size1,
                       int size2)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rContainerItems.size() != rData.size1())
        << "Number of items does not equal data set dimension\n";
    for (std::size_t k = 0; k < rContainerItems.size(); ++k)
    {
        Matrix<double>& r_mat = rContainerItems[k]->GetValue(rVariable);
        r_mat.resize(size1, size2, false);
        for (int i = 0; i < size1; ++i)
            for (int j = 0; j < size2; ++j)
                r_mat(i, j) = rData(k, i * size2 + j);
    }

    KRATOS_CATCH("");
}
} // unnamed namespace.

// template instantiations

template class ContainerComponentIO<ElementsContainerType, ElementType, Flags>;
template class ContainerComponentIO<ElementsContainerType,
                                    ElementType,
                                    Variable<array_1d<double, 3>>,
                                    Variable<double>,
                                    Variable<int>,
                                    Variable<Vector<double>>,
                                    Variable<Matrix<double>>>;

template class ContainerComponentIO<ConditionsContainerType, ConditionType, Flags>;
template class ContainerComponentIO<ConditionsContainerType,
                                    ConditionType,
                                    Variable<array_1d<double, 3>>,
                                    Variable<double>,
                                    Variable<int>,
                                    Variable<Vector<double>>,
                                    Variable<Matrix<double>>>;

template class ContainerComponentIO<NodesContainerType, NodeType, Flags>;
template class ContainerComponentIO<NodesContainerType,
                                    NodeType,
                                    Variable<array_1d<double, 3>>,
                                    Variable<double>,
                                    Variable<int>,
                                    Variable<Vector<double>>,
                                    Variable<Matrix<double>>>;

} // namespace HDF5.
} // namespace Kratos.