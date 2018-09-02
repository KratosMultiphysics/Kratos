#include "custom_io/hdf5_element_data_value_io.h"

#include "utilities/openmp_utils.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/registered_variable_lookup.h"

namespace Kratos
{
namespace HDF5
{
namespace
{
template <class TVariableType, class TFileDataType>
void SetDataBuffer(TVariableType const& rVariable,
                   std::vector<ElementType*> const& rElements,
                   Vector<TFileDataType>& rData);

template <class TVariableType, class TFileDataType>
void SetElementDataValues(TVariableType const& rVariable,
                          Vector<TFileDataType> const& rData,
                          std::vector<ElementType*>& rElements);

template <typename TVariable>
class WriteVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    std::vector<ElementType*>& rElements,
                    File& rFile,
                    std::string const& rPrefix,
                    WriteInfo& rInfo)
    {
        Vector<typename TVariable::Type> data;
        SetDataBuffer(rVariable, rElements, data);
        rFile.WriteDataSet(rPrefix + "/ElementDataValues/" + rVariable.Name(), data, rInfo);
    }
};

template <typename TVariable>
class ReadVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    std::vector<ElementType*>& rElements,
                    File& rFile,
                    std::string const& rPrefix,
                    std::size_t StartIndex,
                    std::size_t BlockSize)
    {
        Vector<typename TVariable::Type> data;
        rFile.ReadDataSet(rPrefix + "/ElementDataValues/" + rVariable.Name(), data,
                          StartIndex, BlockSize);
        SetElementDataValues(rVariable, data, rElements);
    }
};

std::vector<ElementType*> GetElementReferences(ElementsContainerType const& rElements)
{
    std::vector<ElementType*> element_pointers(rElements.size());

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
    {
        auto it = rElements.begin() + i;
        element_pointers[i] = (&(*it));
    }

    return element_pointers;
}
} // unnamed namespace

ElementDataValueIO::ElementDataValueIO(Parameters Settings, File::Pointer pFile)
    : mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix": "",
            "list_of_variables": []
        })");

    Settings.ValidateAndAssignDefaults(default_params);

    mPrefix = Settings["prefix"].GetString();

    const std::size_t num_variables = Settings["list_of_variables"].size();

    if (mVariableNames.size() != num_variables)
        mVariableNames.resize(num_variables);

    for (std::size_t i = 0; i < num_variables; ++i)
        mVariableNames[i] = Settings["list_of_variables"].GetArrayItem(i).GetString();

    KRATOS_CATCH("");
}

void ElementDataValueIO::WriteElementResults(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    WriteInfo info;
    std::vector<ElementType*> local_elements = GetElementReferences(rElements);

    // Write each variable.
    for (const std::string& r_variable_name : mVariableNames)
        RegisteredVariableLookup<Variable<array_1d<double, 3>>,
                                 VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>,
                                 Variable<double>, Variable<int>>(r_variable_name)
            .Execute<WriteVariableFunctor>(local_elements, *mpFile, mPrefix, info);

    // Write block partition.
    WritePartitionTable(*mpFile, mPrefix + "/ElementDataValues", info);

    KRATOS_CATCH("");
}

void ElementDataValueIO::ReadElementResults(ElementsContainerType& rElements)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    std::vector<ElementType*> local_elements = GetElementReferences(rElements);
    std::size_t start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, mPrefix + "/ElementDataValues");

    // Read local data for each variable.
    for (const std::string& r_variable_name : mVariableNames)
        RegisteredVariableLookup<Variable<array_1d<double, 3>>,
                                 VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>,
                                 Variable<double>, Variable<int>>(r_variable_name)
            .Execute<ReadVariableFunctor>(local_elements, *mpFile, mPrefix,
                                          start_index, block_size);

    KRATOS_CATCH("");
}

namespace
{
template <class TVariableType, class TFileDataType>
void SetDataBuffer(TVariableType const& rVariable,
                   std::vector<ElementType*> const& rElements,
                   Vector<TFileDataType>& rData
                   )
{
    KRATOS_TRY;

    if (rData.size() != rElements.size())
        rData.resize(rElements.size(), false);

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
    {
        const ElementType& r_elem = *rElements[i];
        rData[i] = r_elem.GetValue(rVariable);
    }

    KRATOS_CATCH("");
}

template <class TVariableType, class TFileDataType>
void SetElementDataValues(TVariableType const& rVariable,
                          Vector<TFileDataType> const& rData,
                          std::vector<ElementType*>& rElements)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rData.size() != rElements.size())
        << "File data block size (" << rData.size()
        << ") is not equal to number of nodes (" << rElements.size() << ")." << std::endl;

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
        rElements[i]->SetValue(rVariable, rData[i]);

    KRATOS_CATCH("");
}
} // unnamed namespace.
} // namespace HDF5.
} // namespace Kratos.