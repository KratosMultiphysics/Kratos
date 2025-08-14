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

// External includes

// Project includes
#include "custom_utilities/hdf5_data_set_partition_utility.h"

// Application includes

// Include base h
#include "hdf5_tensor_adaptor_io.h"

namespace Kratos
{
class Parameters;

namespace HDF5
{

TensorAdaptorIO::TensorAdaptorIO(
    Parameters Settings,
    File::Pointer pFile)
    : mpFile(pFile)
{
    Parameters default_params(R"(
        {
            "prefix": ""
        })");

    Settings.AddMissingParameters(default_params);

    mPrefix = Settings["prefix"].GetString();
}

void TensorAdaptorIO::Write(
    const std::string& rTensorAdaptorName,
    TensorAdaptorPointerType pTensorAdaptor,
    const Parameters Attributes)
{
    KRATOS_TRY

    std::visit([this, &rTensorAdaptorName, Attributes](auto pContainer){
        auto appended_attribs = Attributes.Clone();
        const auto& r_tensor_shape = pContainer->Shape();
        const auto& dataset_path = mPrefix + rTensorAdaptorName;
        WriteInfo info;
        mpFile->WriteDataSet(dataset_path, pContainer->ViewData().data(), r_tensor_shape.data().begin(), r_tensor_shape.data().end(), info);
        mpFile->WriteAttribute(dataset_path, appended_attribs);
        WritePartitionTable(*mpFile, dataset_path, info);

    }, pTensorAdaptor);

    KRATOS_CATCH("");
}

std::pair<TensorAdaptorIO::TensorAdaptorPointerType, Parameters> TensorAdaptorIO::Read(const std::string& rExpressionName)
{
    KRATOS_TRY

    const auto& dataset_path = mPrefix + rExpressionName;

    KRATOS_ERROR_IF_NOT(mpFile->HasPath(dataset_path))
        << "Path \"" << dataset_path << "\" does not exist.";

    IndexType start_index{}, block_size{};
    if (HasPartitionTable(*mpFile, dataset_path)) {
        // partition table is available for the given data set.
        std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, dataset_path);
    } else {
        // check if the partition table is written as a common table for set of datasets
        // such as in the Variable data writing. This common data set should be having a partition table in the same
        // group.
        const auto group_name_pos = dataset_path.rfind("/");
        const auto& partition_path = dataset_path.substr(0, group_name_pos + 1);

        if (HasPartitionTable(*mpFile, partition_path)) {
            std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, partition_path);
        } else {
            KRATOS_ERROR << "Partition table not found for dataset at \""
                         << dataset_path << "\" either at the dataset location or at \""
                         << partition_path << "\".";
        }
    }

    // get raw data dimensions of h5.
    const auto& h5_dimensions = mpFile->GetDataDimensions(dataset_path);
    DenseVector<unsigned int> tensor_shape(h5_dimensions.size());
    tensor_shape[0] = block_size;
    std::copy(h5_dimensions.begin() + 1, h5_dimensions.end(), tensor_shape.data().begin());

    TensorAdaptorPointerType pTensorAdaptor;
    if (mpFile->HasDataType<bool>(dataset_path)) {
        pTensorAdaptor = Kratos::make_intrusive<TensorAdaptor<bool>>(tensor_shape);
    } else if (mpFile->HasDataType<int>(dataset_path)) {
        pTensorAdaptor = Kratos::make_intrusive<TensorAdaptor<int>>(tensor_shape);
    } else if (mpFile->HasDataType<double>(dataset_path)) {
        pTensorAdaptor = Kratos::make_intrusive<TensorAdaptor<double>>(tensor_shape);
    } else {
        KRATOS_ERROR
                << "Unsupported data set type found at \"" << dataset_path
                << "\". TensorAdaptors only support bool, int and double data types.\n";
    }

    std::visit([this, &dataset_path, &tensor_shape, start_index](auto p_tensor_adaptor) {
        mpFile->ReadDataSet(dataset_path, p_tensor_adaptor->ViewData().data(), tensor_shape.data().begin(),tensor_shape.data().end(), start_index);
    }, pTensorAdaptor);

    auto attributes = mpFile->ReadAttribute(dataset_path);

    return std::make_pair(pTensorAdaptor, attributes);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
