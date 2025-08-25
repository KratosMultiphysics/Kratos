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

    const auto& dataset_path = mPrefix + rTensorAdaptorName;

    std::visit([this, &dataset_path](auto p_tensor_adaptor){
        WriteInfo info;
        const auto& r_tensor_shape = p_tensor_adaptor->Shape();
        if constexpr(std::is_same_v<BareType<decltype(*p_tensor_adaptor)>, TensorAdaptor<bool>>) {
            // HDF5 does not allow writing bool because, HDF5 does not differentiate
            // bool and and unsigned char when written. Therefore, we only allow writing
            // unsigned char. Hence, we need to convert bool array to unsigned char array
            Vector<unsigned char> temp_bools(p_tensor_adaptor->Size());
            const auto span = p_tensor_adaptor->ViewData();
            IndexPartition<IndexType>(temp_bools.size()).for_each([&temp_bools, &span](const auto Index) {
                temp_bools[Index] = static_cast<unsigned char>(span[Index]);
            });
            this->mpFile->WriteDataSet(dataset_path, temp_bools.data().begin(), r_tensor_shape.data().begin(), r_tensor_shape.data().end(), info);
        } else {
            this->mpFile->WriteDataSet(dataset_path, p_tensor_adaptor->ViewData().data(), r_tensor_shape.data().begin(), r_tensor_shape.data().end(), info);
        }
        WritePartitionTable(*this->mpFile, dataset_path, info);
    }, pTensorAdaptor);

    auto appended_attribs = Attributes.Clone();
    mpFile->WriteAttribute(dataset_path, appended_attribs);

    KRATOS_CATCH("");
}

Parameters TensorAdaptorIO::Read(
    const std::string& rTensorAdaptorName,
    TensorAdaptorPointerType pTensorAdaptor)
{
    KRATOS_TRY

    const auto& dataset_path = mPrefix + rTensorAdaptorName;

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
    std::copy(h5_dimensions.begin() + 1, h5_dimensions.end(), tensor_shape.data().begin() + 1);


    std::visit([this, &dataset_path, &tensor_shape, start_index](auto p_tensor_adaptor) {
        const auto& current_shape = p_tensor_adaptor->Shape();

        KRATOS_ERROR_IF_NOT(tensor_shape.size() == current_shape.size())
            << "Number of dimensions mismatch [ shape in the hdf5 = " << tensor_shape
            << ", shape of the tensor = " << current_shape << " ].\n";

        for (IndexType i_dim = 0; i_dim < tensor_shape.size(); ++i_dim) {
            KRATOS_ERROR_IF_NOT(tensor_shape[i_dim] == current_shape[i_dim])
                << "Number of components in dimensions mismatch [ shape in the hdf5 = " << tensor_shape
                << ", shape of the tensor = " << current_shape << " ].\n";
        }

        if constexpr(std::is_same_v<BareType<decltype(*p_tensor_adaptor)>, TensorAdaptor<bool>>) {
            // HDF5 does not allow reading bool because, HDF5 does not differentiate
            // bool and and unsigned char when written. Therefore, we only allow reading
            // unsigned char. Hence, we need to convert unsigned char array to bool.
            Vector<unsigned char> temp_values(p_tensor_adaptor->Size());
            this->mpFile->ReadDataSet(dataset_path, temp_values.data().begin(), tensor_shape.data().begin(),tensor_shape.data().end(), start_index);
            const auto span = p_tensor_adaptor->ViewData();
            IndexPartition<IndexType>(temp_values.size()).for_each([&temp_values, &span](const auto Index) {
                span[Index] = temp_values[Index];
            });
        } else {
            this->mpFile->ReadDataSet(dataset_path, p_tensor_adaptor->ViewData().data(), tensor_shape.data().begin(),tensor_shape.data().end(), start_index);
        }
    }, pTensorAdaptor);

    return mpFile->ReadAttribute(dataset_path);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
