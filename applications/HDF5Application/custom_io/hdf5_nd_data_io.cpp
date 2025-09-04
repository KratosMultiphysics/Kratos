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
#include "custom_utilities/data_type_utilities.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"

// Application includes

// Include base h
#include "hdf5_nd_data_io.h"

namespace Kratos
{
class Parameters;

namespace HDF5
{

NDDataIO::NDDataIO(
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

void NDDataIO::Write(
    const std::string& rNDDataName,
    NDDataPointerType pNDData,
    const Parameters Attributes)
{
    KRATOS_TRY

    const auto& dataset_path = mPrefix + rNDDataName;

    std::visit([this, &dataset_path](auto p_nd_data){
        WriteInfo info;
        const auto& r_nd_data_shape = p_nd_data->Shape();
        if constexpr(std::is_same_v<BareType<decltype(*p_nd_data)>, NDData<bool>>) {
            // HDF5 does not allow writing bool because, HDF5 does not differentiate
            // bool and and unsigned char when written. Therefore, we only allow writing
            // unsigned char. Hence, we need to convert bool array to unsigned char array
            Vector<unsigned char> temp_bools(p_nd_data->Size());
            const auto span = p_nd_data->ViewData();
            IndexPartition<IndexType>(temp_bools.size()).for_each([&temp_bools, &span](const auto Index) {
                temp_bools[Index] = static_cast<unsigned char>(span[Index]);
            });
            this->mpFile->WriteDataSet(dataset_path, temp_bools.data().begin(), r_nd_data_shape.data().begin(), r_nd_data_shape.data().end(), info);
        } else {
            this->mpFile->WriteDataSet(dataset_path, p_nd_data->ViewData().data(), r_nd_data_shape.data().begin(), r_nd_data_shape.data().end(), info);
        }
        WritePartitionTable(*this->mpFile, dataset_path, info);
    }, pNDData);

    auto appended_attribs = Attributes.Clone();
    mpFile->WriteAttribute(dataset_path, appended_attribs);

    KRATOS_CATCH("");
}

std::pair<NDDataIO::NDDataPointerType, Parameters> NDDataIO::Read(
    const std::string& rNDDataName)
{
    KRATOS_TRY

    const auto& dataset_path = mPrefix + rNDDataName;

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

        KRATOS_ERROR_IF(group_name_pos == std::string::npos)
            << "The dataset path = \"" << dataset_path << "\" is not a valid hdf5 path. HDF5 path should always start with \"/\"\n.";

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
    std::copy(h5_dimensions.begin() + 1, h5_dimensions.end(), tensor_shape.data().begin() + 1);

    NDDataPointerType p_nd_data;
    if (mpFile->HasDataType<unsigned char>(dataset_path)) {
        p_nd_data = Kratos::make_shared<NDData<unsigned char>>(tensor_shape);
    } else if (mpFile->HasDataType<int>(dataset_path)) {
        p_nd_data = Kratos::make_shared<NDData<int>>(tensor_shape);
    } else if (mpFile->HasDataType<double>(dataset_path)) {
        p_nd_data = Kratos::make_shared<NDData<double>>(tensor_shape);
    } else {
        KRATOS_ERROR
                << "Unsupported data set type found at \"" << dataset_path
                << "\". NDDatas only support uchar, int and double data types.\n";
    }

    std::visit([this, &dataset_path, &tensor_shape, start_index](auto p_nd_data) {
        if constexpr(!std::is_same_v<BareType<decltype(*p_nd_data)>, NDData<bool>>) {
            this->mpFile->ReadDataSet(dataset_path, p_nd_data->ViewData().data(), tensor_shape.data().begin(),tensor_shape.data().end(), start_index);
        }
    }, p_nd_data);

    auto attributes = mpFile->ReadAttribute(dataset_path);
    return std::make_pair(p_nd_data, attributes);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
