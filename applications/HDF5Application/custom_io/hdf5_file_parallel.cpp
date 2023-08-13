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
//                  Suneth Warnakulasuriya
//

// System includes
#include <numeric>

// Project includes
#include "includes/kratos_parameters.h"
#include "utilities/builtin_timer.h"
#include "input_output/logger.h"
#include "utilities/data_type_traits.h"

// Application includes
#include "custom_utilities/h5_data_type_traits.h"

// Include base h
#include "hdf5_file_parallel.h"

namespace Kratos
{
namespace HDF5
{

FileParallel::FileParallel(Parameters& rSettings) : File(rSettings)
{
}

FileParallel::FileParallel(
    const DataCommunicator& rDataCommunicator,
    Parameters Settings)
    : File(rDataCommunicator, Settings)
{
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath,
                              Vector<array_1d<double, 3>>& rData,
                              unsigned StartIndex,
                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<int>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<double>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<array_1d<double, 3>>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                              Matrix<int>& rData,
                                              unsigned StartIndex,
                                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                              Matrix<double>& rData,
                                              unsigned StartIndex,
                                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

template <class TDataSetType>
void FileParallel::WriteDataSetImpl(
    const std::string& rPath,
    const TDataSetType& rData,
    DataTransferMode Mode,
    WriteInfo& rInfo)
{
    KRATOS_TRY

    using type_trait = DataTypeTraits<TDataSetType>;

    static_assert(type_trait::IsContiguous, "HDF5File can only write contiguous data sets.");

    constexpr auto local_dimension = type_trait::Dimension;

    constexpr auto global_dimension = (local_dimension == 1 ? 1 : 2);

    BuiltinTimer timer;

    // Create any missing subpaths.
    auto pos = rPath.find_last_of('/');
    if (pos != 0) {// Skip if last '/' is root.
        std::string sub_path = rPath.substr(0, pos);
        AddPath(sub_path);
    }

    // Initialize data space dimensions.
    std::vector<hsize_t> local_shape(local_dimension);
    type_trait::Shape(rData, local_shape.data(), local_shape.data() + local_dimension);

    const hsize_t number_of_local_primitive_data_values = type_trait::Size(rData);

    const auto& r_data_communicator = GetDataCommunicator();
    // get the maximized dimensions of the underlying data. Max is taken because,
    // there can be empty ranks which will give wrong sizes in the case of dynamic
    // data types.
    if constexpr(local_dimension >= 2) {
        if constexpr(type_trait::template IsDimensionDynamic<1>()) {
            // this is the matrix version. Hence it is better to get the max size
            // from all ranks.
            const auto max_size = r_data_communicator.MaxAll(local_shape[1]);

            // now check every non-empty ranks have the same sizes since this dimension
            // is a dynamic dimension.
            KRATOS_ERROR_IF(number_of_local_primitive_data_values > 0 && max_size != local_shape[1])
                << "Mismatching shapes found in different ranks. All ranks should have the same shapes in data sets.";

            local_shape[1] = max_size;
        }
    }

    // local_reduced_shape holds the max 2d flattened shape if the local_shape dimensions
    // are higher than 2.
    std::vector<hsize_t> global_shape(global_dimension, 0), local_reduced_shape(global_dimension, 0), local_shape_start(global_dimension, 0);

    // get total number of items to be written in the data set to the first dimension.
    global_shape[0] = r_data_communicator.SumAll(local_shape[0]);
    local_reduced_shape[0] = local_shape[0];

    if constexpr(global_dimension > 1) {
        // flattens higher dimensions into one since we write matrices which is the highest dimension
        // supported by paraview for visualization
        global_shape[1] = std::accumulate(local_shape.begin() + 1, local_shape.end(), hsize_t{1}, std::multiplies<hsize_t>());
        local_reduced_shape[1] = global_shape[1];
    }

    // Set the data type.
    hid_t dtype_id = Internals::GetPrimitiveH5Type<TDataSetType>();

    // Create and write the data set.
    hid_t dset_id{}, fspace_id{};
    GetDataSet<typename type_trait::PrimitiveType>(dset_id, fspace_id, global_shape, rPath);

    if (r_data_communicator.IsDistributed()) {
        if (Mode == DataTransferMode::collective) {
            local_shape_start[0] = r_data_communicator.ScanSum(local_reduced_shape[0]) - local_reduced_shape[0];
        }

        if (Mode == DataTransferMode::collective || number_of_local_primitive_data_values > 0) {
            hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
            if (Mode == DataTransferMode::collective) {
                H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
            } else {
                H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
            }

            // select the local hyperslab
            H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, local_shape_start.data(), nullptr, local_reduced_shape.data(), nullptr);
            hid_t mspace_id = H5Screate_simple(global_dimension, local_reduced_shape.data(), nullptr);
            if (number_of_local_primitive_data_values > 0) {
                KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, type_trait::GetContiguousData(rData)) < 0)
                    << "H5Dwrite failed." << std::endl;
            } else {
                KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, nullptr) < 0)
                    << "H5Dwrite failed. Please ensure global data set is non-empty."
                    << std::endl;
            }

            KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
            KRATOS_ERROR_IF(H5Sclose(mspace_id) < 0) << "H5Sclose failed." << std::endl;
        }
    } else {
        if (number_of_local_primitive_data_values > 0) {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, type_trait::GetContiguousData(rData)) < 0)
                << "H5Dwrite failed." << std::endl;
        } else {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, nullptr) < 0)
                << "H5Dwrite failed. Please ensure global data set is non-empty."
                << std::endl;
        }
    }

    KRATOS_ERROR_IF(H5Sclose(fspace_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    // Set the write info.
    rInfo.StartIndex = local_shape_start[0];
    rInfo.BlockSize = local_reduced_shape[0];
    rInfo.TotalSize = global_shape[0];

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Write time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

template <class TDataSetType>
void FileParallel::ReadDataSetImpl(
    const std::string& rPath,
    TDataSetType& rData,
    unsigned StartIndex,
    unsigned BlockSize,
    DataTransferMode Mode)
{
    KRATOS_TRY;

    using type_trait = DataTypeTraits<TDataSetType>;

    static_assert(type_trait::IsContiguous, "HDF5File can only write contiguous data sets.");

    constexpr auto local_dimension = type_trait::Dimension;

    constexpr auto global_dimension = (local_dimension == 1 ? 1 : 2);

    BuiltinTimer timer;
    // Check that full path exists.
    KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
        << "Path is not a data set: " << rPath << std::endl;

    const auto& file_space_dims = GetDataDimensions(rPath);

    // Check consistency of file's data set dimensions.
    KRATOS_ERROR_IF(file_space_dims.size() != global_dimension)
        << "Invalid data set dimension." << std::endl;
    KRATOS_ERROR_IF(StartIndex + BlockSize > file_space_dims[0])
        << "StartIndex (" << StartIndex << ") + BlockSize (" << BlockSize
        << ") > size of data set (" << file_space_dims[0] << ")." << std::endl;

    std::vector<hsize_t> memory_space_dims(local_dimension);
    type_trait::Shape(rData, memory_space_dims.data(), memory_space_dims.data() + local_dimension);

    if constexpr(local_dimension >= 2) {
        if constexpr(type_trait::template IsDimensionDynamic<1>()) {
            const auto v = std::accumulate(memory_space_dims.begin() + 2, memory_space_dims.end(), hsize_t{1}, std::multiplies<hsize_t>());
            KRATOS_ERROR_IF_NOT(file_space_dims[1] % v == 0) << "Size mismatch with memory space and file space.";
            memory_space_dims[1] = file_space_dims[1] / v;
        }
    }
    memory_space_dims[0] = BlockSize;

    // now reshape the memory space data
    type_trait::Reshape(rData, memory_space_dims.data(), memory_space_dims.data() + local_dimension);

    std::vector<hsize_t> local_reduced_space_dims(file_space_dims.begin(), file_space_dims.end()), local_space_start(global_dimension, 0);
    local_reduced_space_dims[0] = BlockSize;
    local_space_start[0] = StartIndex;

    if constexpr(std::is_same_v<typename type_trait::PrimitiveType, int>) {
        KRATOS_ERROR_IF_NOT(HasIntDataType(rPath))
            << "Data type is not int: " << rPath << std::endl;
    } else if constexpr(std::is_same_v<typename type_trait::PrimitiveType, double>) {
        KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
            << "Data type is not float: " << rPath << std::endl;
    } else {
        static_assert(!std::is_same_v<TDataSetType, TDataSetType>, "Unsupported data type.");
    }

    // Set the data type.
    hid_t dtype_id = Internals::GetPrimitiveH5Type<TDataSetType>();

    hid_t file_id = GetFileId();
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    if (Mode == DataTransferMode::collective) {
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    } else {
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
    }

    hid_t dset_id = H5Dopen(file_id, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
    hid_t file_space_id = H5Dget_space(dset_id);
    hid_t mem_space_id = H5Screate_simple(global_dimension, local_reduced_space_dims.data(), nullptr);
    KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, local_space_start.data(),
                                        nullptr, local_reduced_space_dims.data(), nullptr) < 0)
        << "H5Sselect_hyperslab failed." << std::endl;
    if (type_trait::Size(rData) > 0) {
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, dxpl_id, type_trait::GetContiguousData(rData)) < 0)
            << "H5Dread failed." << std::endl;
    } else {
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, dxpl_id, nullptr) < 0)
            << "H5Dread failed." << std::endl;
    }

    KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(file_space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(mem_space_id) < 0) << "H5Sclose failed." << std::endl;

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Read time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

// template instantiations
#ifndef KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS
#define KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(...)                                                                      \
template void FileParallel::WriteDataSetImpl(const std::string&, const __VA_ARGS__&, DataTransferMode, WriteInfo& rInfo);   \
template void FileParallel::ReadDataSetImpl(const std::string&, __VA_ARGS__&, unsigned, unsigned, DataTransferMode);
#endif

KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<int>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<double>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<array_1d<double, 3>>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<array_1d<double, 4>>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<array_1d<double, 6>>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<array_1d<double, 9>>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Matrix<int>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Matrix<double>);

#undef KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS

} // namespace HDF5.
} // namespace Kratos.