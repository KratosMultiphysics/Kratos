//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya, https://github.com/sunethwarna
//

// System includes
#include <tuple>

// External includes

// Project includes
#include "includes/communicator.h"
#include "includes/kratos_parameters.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "custom_utilities/local_ghost_splitting_utility.h"
#include "custom_utilities/registered_component_lookup.h"

// Include base h
#include "custom_io/hdf5_container_gauss_point_output.h"

namespace Kratos
{
namespace HDF5
{
namespace
{
template<typename TDataType>
std::size_t GetSize(std::vector<TDataType>& rData);

template<typename TDataType>
void FlattenData(
    Vector<double>& rFlattedData,
    const std::vector<TDataType>& rData);

template <typename TContainerType, typename TDataType>
void SetDataBuffer(Variable<TDataType> const&,
                   TContainerType&,
                   Matrix<double>&,
                   int&,
                   const DataCommunicator& rDataCommunicator,
                   const ProcessInfo&);

template <typename TContainerType>
class ContainerItemGaussPointValuesOutput
{
public:
    template <typename TVariableType>
    class WriteGaussPointValuesFunctor
    {
    public:
        void operator()(TVariableType const& rVariable,
                        TContainerType& rContainer,
                        File& rFile,
                        std::string const& rPath,
                        WriteInfo& rInfo,
                        const DataCommunicator& rDataCommunicator,
                        const ProcessInfo& rProcessInfo)
        {
            Matrix<double> data;
            int number_of_gauss_points;
            SetDataBuffer<TContainerType, typename TVariableType::Type>(rVariable, rContainer, data, number_of_gauss_points, rDataCommunicator, rProcessInfo);
            KRATOS_WARNING_IF("WriteGaussPointValuesFunctor", number_of_gauss_points == 0) << "No integration point values were found for " << rVariable.Name() << ".\n";
            rFile.WriteDataSet(rPath + "/" + rVariable.Name(), data, rInfo);
            rFile.WriteAttribute(rPath + "/" + rVariable.Name(), "NumberOfGaussPoints", number_of_gauss_points);
        }
    };
};

template <typename TVariableType>
class WriteElementGaussPointValues
    : public ContainerItemGaussPointValuesOutput<ElementsContainerType>::WriteGaussPointValuesFunctor<TVariableType>
{
};

template <typename TVariableType>
class WriteConditionGaussPointValues
    : public ContainerItemGaussPointValuesOutput<ConditionsContainerType>::WriteGaussPointValuesFunctor<TVariableType>
{
};

} // unnamed namespace

template <typename TContainerType, typename... TComponents>
ContainerGaussPointOutput<TContainerType, TComponents...>::ContainerGaussPointOutput(
    Parameters Settings, File::Pointer pFile, const std::string& rPath)
    : mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix": "",
            "list_of_variables": []
        })");

    Settings.ValidateAndAssignDefaults(default_params);

    mVariablePath = Settings["prefix"].GetString() + rPath;

    const std::size_t num_components = Settings["list_of_variables"].size();

    if (mVariableNames.size() != num_components)
        mVariableNames.resize(num_components);

    for (std::size_t i = 0; i < num_components; ++i)
        mVariableNames[i] = Settings["list_of_variables"].GetArrayItem(i).GetString();

    KRATOS_CATCH("");
}

// Adding element other variable WriteRegisteredGaussPointValues template definition
template <>
template <typename... Targs>
void ContainerGaussPointOutput<ElementsContainerType,
                          Variable<array_1d<double, 3>>,
                          Variable<double>,
                          Variable<int>,
                          Variable<Vector<double>>,
                          Variable<Matrix<double>>>::WriteRegisteredGaussPointValues(
                              const std::string& rComponentName,
                              Targs&... args)
{
    RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                              Variable<Vector<double>>, Variable<Matrix<double>>>(rComponentName)
        .Execute<WriteElementGaussPointValues>(args...);
}

// Adding condition other variable WriteRegisteredGaussPointValues template definition
template <>
template <typename... Targs>
void ContainerGaussPointOutput<ConditionsContainerType,
                          Variable<array_1d<double, 3>>,
                          Variable<double>,
                          Variable<int>,
                          Variable<Vector<double>>,
                          Variable<Matrix<double>>>::WriteRegisteredGaussPointValues(
                              const std::string& rComponentName,
                              Targs&... args)
{
    RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                              Variable<Vector<double>>, Variable<Matrix<double>>>(rComponentName)
        .Execute<WriteConditionGaussPointValues>(args...);
}


template <typename TContainerType, typename... TComponents>
void ContainerGaussPointOutput<TContainerType, TComponents...>::WriteContainerGaussPointsValues(
    TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    WriteInfo info;

    // Write each variable.
    for (const std::string& r_component_name : mVariableNames)
        WriteRegisteredGaussPointValues(
            r_component_name, rContainer, *mpFile,
            mVariablePath, info, rDataCommunicator, rProcessInfo);

    // Write block partition.
    WritePartitionTable(*mpFile, mVariablePath, info);

    KRATOS_CATCH("");
}

namespace
{

template<typename TDataType>
std::size_t GetSize(std::vector<TDataType>& rData)
{
    return rData.size();
}

template<>
std::size_t GetSize(std::vector<array_1d<double, 3>>& rData)
{
    return rData.size() * 3;
}

template<>
std::size_t GetSize(std::vector<Vector<double>>& rData)
{
    if (rData.size() > 0) {
        return rData.size() * rData[0].size();
    } else {
        return 0;
    }
}

template<>
std::size_t GetSize(std::vector<Matrix<double>>& rData)
{
    if (rData.size() > 0) {
        return rData.size() * rData[0].data().size();
    } else {
        return 0;
    }
}

template<typename TDataType>
void FlattenData(
    Vector<double>& rFlattedData,
    const std::vector<TDataType>& rData)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rFlattedData.size() != rData.size())
        << "Size mismatch. [ rFlattedData.size() = " << rFlattedData.size()
        << ", rData.size() = " << rData.size() << " ].\n";

    std::size_t local_index = 0;
    for (const auto& r_data : rData) {
        rFlattedData[local_index++] = r_data;
    }

    KRATOS_CATCH("");
}

template<>
void FlattenData(
    Vector<double>& rFlattedData,
    const std::vector<array_1d<double, 3>>& rData)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rFlattedData.size() != rData.size() * 3)
        << "Size mismatch. [ rFlattedData.size() = " << rFlattedData.size()
        << ", rData.size() * 3 = " << rData.size() * 3 << " ].\n";

    std::size_t local_index = 0;
    for (const auto& r_data : rData) {
        rFlattedData[local_index++] = r_data[0];
        rFlattedData[local_index++] = r_data[1];
        rFlattedData[local_index++] = r_data[2];
    }

    KRATOS_CATCH("");
}

template<>
void FlattenData(
    Vector<double>& rFlattedData,
    const std::vector<Vector<double>>& rData)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rFlattedData.size() != rData.size() * rData[0].size())
        << "Size mismatch. [ rFlattedData.size() = " << rFlattedData.size()
        << ", rData.size() * rData[0].size() = " << rData.size() * rData[0].size() << " ].\n";

    std::size_t local_index = 0;
    for (const auto& r_data : rData) {
        for (std::size_t i = 0; i < r_data.size(); ++i) {
            rFlattedData[local_index++] = r_data[i];
        }
    }

    KRATOS_CATCH("");
}

template<>
void FlattenData(
    Vector<double>& rFlattedData,
    const std::vector<Matrix<double>>& rData)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rFlattedData.size() != rData.size() * rData[0].data().size())
        << "Size mismatch. [ rFlattedData.size() = " << rFlattedData.size()
        << ", rData.size() * rData[0].data().size() = " << rData.size() * rData[0].data().size() << " ].\n";

    std::size_t local_index = 0;
    for (const auto& r_data : rData) {
        for (std::size_t i = 0; i < r_data.data().size(); ++i) {
            rFlattedData[local_index++] = r_data.data()[i];
        }
    }

    KRATOS_CATCH("");
}

template <typename TContainerType, typename TDataType>
void SetDataBuffer(Variable<TDataType> const& rVariable,
                   TContainerType& rContainer,
                   Matrix<double>& rData,
                   int& NumberOfGaussPoints,
                   const DataCommunicator& rDataCommunicator,
                   const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    std::size_t flat_data_size = 0;
    NumberOfGaussPoints = 0;
    if (rContainer.size() > 0) {
        std::vector<TDataType> values;
        rContainer.begin()->CalculateOnIntegrationPoints(rVariable, values, rProcessInfo);
        NumberOfGaussPoints = values.size();
        flat_data_size = GetSize(values);
    }

    NumberOfGaussPoints = rDataCommunicator.MaxAll(NumberOfGaussPoints);
    flat_data_size = rDataCommunicator.MaxAll(flat_data_size);

    rData.resize(rContainer.size(), flat_data_size);

    struct tls_type
    {
        explicit tls_type(const std::size_t FlatDataSize) : mFlattedGaussPointValues(FlatDataSize) {}
        std::vector<TDataType> mGaussPointValues;
        Vector<double> mFlattedGaussPointValues;
    };

    IndexPartition<int>(rContainer.size()).for_each(tls_type(flat_data_size), [&](const int ItemIndex, tls_type& rTLS){
        (rContainer.begin() + ItemIndex)->CalculateOnIntegrationPoints(rVariable, rTLS.mGaussPointValues, rProcessInfo);
        FlattenData(rTLS.mFlattedGaussPointValues, rTLS.mGaussPointValues);
        row(rData, ItemIndex) = rTLS.mFlattedGaussPointValues;
    });

    KRATOS_CATCH("");
}

} // unnamed namespace.

// template instantiations

template class ContainerGaussPointOutput<ElementsContainerType,
                                    Variable<array_1d<double, 3>>,
                                    Variable<double>,
                                    Variable<int>,
                                    Variable<Vector<double>>,
                                    Variable<Matrix<double>>>;

template class ContainerGaussPointOutput<ConditionsContainerType,
                                    Variable<array_1d<double, 3>>,
                                    Variable<double>,
                                    Variable<int>,
                                    Variable<Vector<double>>,
                                    Variable<Matrix<double>>>;

} // namespace HDF5.
} // namespace Kratos.
