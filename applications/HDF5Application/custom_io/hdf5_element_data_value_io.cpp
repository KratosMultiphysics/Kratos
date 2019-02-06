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
template <typename TDataType>
void SetDataBuffer(Variable<TDataType> const&,
                   std::vector<ElementType*> const&,
                   Vector<TDataType>&);

void SetDataBuffer(Variable<Vector<double>> const&,
                   std::vector<ElementType*> const&,
                   Matrix<double>&);

void SetDataBuffer(Variable<Matrix<double>> const&,
                   std::vector<ElementType*> const&,
                   Matrix<double>&);


template <typename TDataType>
void SetElementDataValues(Variable<TDataType> const&,
                          Vector<TDataType> const&,
                          std::vector<ElementType*>&);

void SetElementDataValues(Variable<Vector<double>> const&,
                          Matrix<double> const&,
                          std::vector<ElementType*>&);

void SetElementDataValues(Variable<Matrix<double>> const&,
                          Matrix<double> const&,
                          std::vector<ElementType*>&,
                          int,
                          int);

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

template <>
class WriteVariableFunctor<Variable<Vector<double>>>
{
public:
    void operator()(Variable<Vector<double>> const& rVariable,
                    std::vector<ElementType*>& rElements,
                    File& rFile,
                    std::string const& rPrefix,
                    WriteInfo& rInfo)
    {
        Matrix<double> data;
        SetDataBuffer(rVariable, rElements, data);
        rFile.WriteDataSet(rPrefix + "/ElementDataValues/" + rVariable.Name(), data, rInfo);
    }
};

template <>
class WriteVariableFunctor<Variable<Matrix<double>>>
{
public:
    void operator()(Variable<Matrix<double>> const& rVariable,
                    std::vector<ElementType*>& rElements,
                    File& rFile,
                    std::string const& rPrefix,
                    WriteInfo& rInfo)
    {
        Matrix<double> data;
        SetDataBuffer(rVariable, rElements, data);
        rFile.WriteDataSet(rPrefix + "/ElementDataValues/" + rVariable.Name(), data, rInfo);
        const int size1 = rElements.front()->GetValue(rVariable).size1();
        const int size2 = rElements.front()->GetValue(rVariable).size2();
        rFile.WriteAttribute(
            rPrefix + "/ElementDataValues/" + rVariable.Name(), "Size1", size1);
        rFile.WriteAttribute(
            rPrefix + "/ElementDataValues/" + rVariable.Name(), "Size2", size2);
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

template <>
class ReadVariableFunctor<Variable<Vector<double>>>
{
public:
    void operator()(Variable<Vector<double>> const& rVariable,
                    std::vector<ElementType*>& rElements,
                    File& rFile,
                    std::string const& rPrefix,
                    std::size_t StartIndex,
                    std::size_t BlockSize)
    {
        Matrix<double> data;
        rFile.ReadDataSet(rPrefix + "/ElementDataValues/" + rVariable.Name(), data,
                          StartIndex, BlockSize);
        SetElementDataValues(rVariable, data, rElements);
    }
};

template <>
class ReadVariableFunctor<Variable<Matrix<double>>>
{
public:
    void operator()(Variable<Matrix<double>> const& rVariable,
                    std::vector<ElementType*>& rElements,
                    File& rFile,
                    std::string const& rPrefix,
                    std::size_t StartIndex,
                    std::size_t BlockSize)
    {
        Matrix<double> data;
        rFile.ReadDataSet(rPrefix + "/ElementDataValues/" + rVariable.Name(), data,
                          StartIndex, BlockSize);
        int size1, size2;
        rFile.ReadAttribute(
            rPrefix + "/ElementDataValues/" + rVariable.Name(), "Size1", size1);
        rFile.ReadAttribute(
            rPrefix + "/ElementDataValues/" + rVariable.Name(), "Size2", size2);
        SetElementDataValues(rVariable, data, rElements, size1, size2);
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
                                 Variable<double>,
                                 Variable<int>,
                                 Variable<Vector<double>>,
                                 Variable<Matrix<double>>>(r_variable_name)
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
                                 Variable<double>,
                                 Variable<int>,
                                 Variable<Vector<double>>,
                                 Variable<Matrix<double>>>(r_variable_name)
            .Execute<ReadVariableFunctor>(local_elements, *mpFile, mPrefix,
                                          start_index, block_size);

    KRATOS_CATCH("");
}

namespace
{
template <typename TDataType>
void SetDataBuffer(Variable<TDataType> const& rVariable,
                   std::vector<ElementType*> const& rElements,
                   Vector<TDataType>& rData)
{
    KRATOS_TRY;

    rData.resize(rElements.size(), false);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
    {
        const ElementType& r_element = *rElements[i];
        rData[i] = r_element.GetValue(rVariable);
    }

    KRATOS_CATCH("");
}

void SetDataBuffer(Variable<Vector<double>> const& rVariable,
                   std::vector<ElementType*> const& rElements,
                   Matrix<double>& rData)
{
    KRATOS_TRY;

    rData.resize(rElements.size(), rElements.front()->GetValue(rVariable).size(), false);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
    {
        const ElementType& r_element = *rElements[i];
        Vector<double> const& r_vec = r_element.GetValue(rVariable);
        KRATOS_ERROR_IF(r_vec.size() != rData.size2())
            << "Invalid dimension.\n";
        for (std::size_t j = 0; j < rData.size2(); ++j)
            rData(i, j) = r_vec(j);
    }

    KRATOS_CATCH("");
}

void SetDataBuffer(Variable<Matrix<double>> const& rVariable,
                   std::vector<ElementType*> const& rElements,
                   Matrix<double>& rData)
{
    KRATOS_TRY;

    rData.resize(rElements.size(), rElements.front()->GetValue(rVariable).data().size(), false);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
    {
        const ElementType& r_element = *rElements[i];
        const auto& r_data = r_element.GetValue(rVariable).data();
        KRATOS_ERROR_IF(r_data.size() != rData.size2())
            << "Invalid dimension.\n";
        for (std::size_t j = 0; j < rData.size2(); ++j)
            rData(i, j) = r_data[j];
    }

    KRATOS_CATCH("");
}

template <typename TDataType>
void SetElementDataValues(Variable<TDataType> const& rVariable,
                          Vector<TDataType> const& rData,
                          std::vector<ElementType*>& rElements)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rElements.size() != rData.size())
        << "Number of elements does not equal data set dimension\n";
    for (std::size_t i = 0; i < rElements.size(); ++i)
        rElements[i]->GetValue(rVariable) = rData[i];

    KRATOS_CATCH("");
}

void SetElementDataValues(Variable<Vector<double>> const& rVariable,
                         Matrix<double> const& rData,
                         std::vector<ElementType*>& rElements)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rElements.size() != rData.size1())
        << "Number of elements does not equal data set dimension\n";
    for (std::size_t i = 0; i < rElements.size(); ++i)
    {
        Vector<double>& r_vec = rElements[i]->GetValue(rVariable);
        r_vec.resize(rData.size2(), false);
        for (std::size_t j = 0; j < rData.size2(); ++j)
            r_vec(j) = rData(i, j);
    }

    KRATOS_CATCH("");
}

void SetElementDataValues(Variable<Matrix<double>> const& rVariable,
                          Matrix<double> const& rData,
                          std::vector<ElementType*>& rElements,
                          int size1,
                          int size2)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rElements.size() != rData.size1())
        << "Number of elements does not equal data set dimension\n";
    for (std::size_t k = 0; k < rElements.size(); ++k)
    {
        Matrix<double>& r_mat = rElements[k]->GetValue(rVariable);
        r_mat.resize(size1, size2, false);
        for (int i = 0; i < size1; ++i)
            for (int j = 0; j < size2; ++j)
                r_mat(i, j) = rData(k, i * size2 + j);
    }

    KRATOS_CATCH("");
}
} // unnamed namespace.
} // namespace HDF5.
} // namespace Kratos.