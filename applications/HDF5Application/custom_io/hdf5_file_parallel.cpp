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
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath,
                              Vector<array_1d<double, 3>>& rData,
                              unsigned StartIndex,
                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSet(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::collective);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<int>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<double>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                        Vector<array_1d<double, 3>>& rData,
                                        unsigned StartIndex,
                                        unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                              Matrix<int>& rData,
                                              unsigned StartIndex,
                                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
    KRATOS_CATCH("");
}

void FileParallel::ReadDataSetIndependent(const std::string& rPath,
                                              Matrix<double>& rData,
                                              unsigned StartIndex,
                                              unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize, DataTransferMode::independent);
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

    using dataset_type_trait = DataTypeTraits<TDataSetType>;

    using value_type_trait = DataTypeTraits<typename dataset_type_trait::ValueType>;

    constexpr auto local_dimension = dataset_type_trait::Dimension;

    static_assert(local_dimension <= 2 && local_dimension > 0, "HDF5File can only write data sets with dimension in range [1, 2].");

    constexpr auto global_dimension = (local_dimension == 1 ? 1 : 2);

    BuiltinTimer timer;

    // Expects a valid free path.
    KRATOS_ERROR_IF(HasPath(rPath)) << "Path already exists: " << rPath << std::endl;

    // Create any missing subpaths.
    auto pos = rPath.find_last_of('/');
    if (pos != 0) {// Skip if last '/' is root.
        std::string sub_path = rPath.substr(0, pos);
        AddPath(sub_path);
    }

    // Initialize data space dimensions.
    std::vector<hsize_t> local_shape(local_dimension);
    dataset_type_trait::Shape(rData, local_shape.data(), local_shape.data() + local_dimension);

    const auto& r_data_communicator = GetDataCommunicator();
    // get the maximized dimensions of the underlying data. Max is taken because,
    // there can be empty ranks which will give wrong sizes in the case of dynamic
    // data types.
    std::vector<hsize_t> max_local_shape(local_shape.begin() + 1, local_shape.end());
    if constexpr(local_dimension == 2 && value_type_trait::Dimension == 0) {
        // this is the matrix version. Hence it is better to get the max size
        // from all ranks.
        max_local_shape[0] = r_data_communicator.MaxAll(max_local_shape[0]);
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
        global_shape[1] = max_local_shape[0];
        local_reduced_shape[1] = global_shape[1];
    }

    const hsize_t number_of_local_primitive_data_values = std::accumulate(local_reduced_shape.begin(), local_reduced_shape.end(), hsize_t{1}, std::multiplies<hsize_t>());

    if (Mode == DataTransferMode::collective) {
        local_shape_start[0] = r_data_communicator.ScanSum(local_reduced_shape[0]) - local_reduced_shape[0];
    }

    // Set the data type.
    hid_t dtype_id = Internals::GetH5DataType<typename dataset_type_trait::PrimitiveType>();

    // Create and write the data set.
    hid_t file_id = GetFileId();
    hid_t fspace_id = H5Screate_simple(global_dimension, global_shape.data(), nullptr);
    // H5Dcreate() must be called collectively for both collective and
    // independent write.
    hid_t dset_id = H5Dcreate(file_id, rPath.c_str(), dtype_id, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;

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
        if (dataset_type_trait::Size(rData) > 0) {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, dataset_type_trait::GetContiguousData(rData)) < 0)
                << "H5Dwrite failed." << std::endl;
        } else {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, nullptr) < 0)
                << "H5Dwrite failed. Please ensure global data set is non-empty."
                << std::endl;
        }

        KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(mspace_id) < 0) << "H5Sclose failed." << std::endl;
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

template <class T>
void FileParallel::ReadDataSetVectorImpl(
    const std::string& rPath,
    Vector<T>& rData,
    unsigned StartIndex,
    unsigned BlockSize,
    DataTransferMode Mode)
{
    KRATOS_TRY;

    using dataset_type_trait = DataTypeTraits<Vector<T>>;

    using value_type_trait = DataTypeTraits<typename dataset_type_trait::ValueType>;

    constexpr auto local_dimension = dataset_type_trait::Dimension;

    static_assert(local_dimension <= 2 && local_dimension > 0, "HDF5File can only write data sets with dimension in range [1, 2].");

    constexpr auto global_dimension = (local_dimension == 1 ? 1 : 2);

    BuiltinTimer timer;
    // Check that full path exists.
    KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
        << "Path is not a data set: " << rPath << std::endl;

    constexpr bool is_int_type = std::is_same<int, T>::value;
    constexpr bool is_double_type = std::is_same<double, T>::value;
    constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
    constexpr unsigned ndims = (!is_array_1d_type) ? 1 : 2;

    auto file_space_dims = GetDataDimensions(rPath);

    // Check consistency of file's data set dimensions.
    KRATOS_ERROR_IF(file_space_dims.size() != global_dimension)
        << "Invalid data set dimension." << std::endl;
    KRATOS_ERROR_IF(StartIndex + BlockSize > file_space_dims[0])
        << "StartIndex (" << StartIndex << ") + BlockSize (" << BlockSize
        << ") > size of data set (" << file_space_dims[0] << ")." << std::endl;

    std::vector<hsize_t> local_space_dims(file_space_dims.begin(), file_space_dims.end()), local_space_start(global_dimension, 0);
    local_space_dims[0] = BlockSize;
    local_space_start[0] = StartIndex;

    // now reshape the memory space data
    dataset_type_trait::Reshape(rData, local_space_dims.data(), local_space_dims.data() + local_dimension);

    if constexpr(std::is_same_v<typename dataset_type_trait::PrimitiveType, int>) {
        KRATOS_ERROR_IF_NOT(HasIntDataType(rPath))
            << "Data type is not int: " << rPath << std::endl;
    } else if constexpr(std::is_same_v<typename dataset_type_trait::PrimitiveType, double>) {
        KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
            << "Data type is not float: " << rPath << std::endl;
    } else {
        static_assert(!std::is_same_v<T, T>, "Unsupported data type.");
    }

    // Set the data type.
    hid_t dtype_id = Internals::GetH5DataType<typename dataset_type_trait::PrimitiveType>();

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
    hid_t mem_space_id = H5Screate_simple(global_dimension, local_space_dims.data(), nullptr);
    KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, local_space_start.data(),
                                        nullptr, local_space_dims.data(), nullptr) < 0)
        << "H5Sselect_hyperslab failed." << std::endl;
    if (dataset_type_trait::Size(rData) > 0) {
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, dxpl_id, dataset_type_trait::GetContiguousData(rData)) < 0)
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

template <class T>
void FileParallel::ReadDataSetMatrixImpl(const std::string& rPath,
                                         Matrix<T>& rData,
                                         unsigned StartIndex,
                                         unsigned BlockSize,
                                         DataTransferMode Mode)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    // Check that full path exists.
    KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
        << "Path is not a data set: " << rPath << std::endl;

    const unsigned ndims = 2;

    // Check consistency of file's data set dimensions.
    std::vector<unsigned> file_space_dims = GetDataDimensions(rPath);
    KRATOS_ERROR_IF(file_space_dims.size() != ndims)
        << "Invalid data set dimension." << std::endl;
    KRATOS_ERROR_IF(StartIndex + BlockSize > file_space_dims[0])
        << "StartIndex (" << StartIndex << ") + BlockSize (" << BlockSize
        << ") > size of data set (" << file_space_dims[0] << ")." << std::endl;

    if (rData.size1() != BlockSize || rData.size2() != file_space_dims[1])
        rData.resize(BlockSize, file_space_dims[1], false);

    hsize_t local_mem_dims[ndims];
    local_mem_dims[0] = rData.size1();
    local_mem_dims[1] = rData.size2();

    // Set global position where local data set starts.
    hsize_t local_start[ndims] = {0};
    local_start[0] = StartIndex;

    // Set the data type.
    hid_t dtype_id = Internals::GetScalarDataType(rData);
    if (dtype_id == H5T_NATIVE_INT)
    {
        KRATOS_ERROR_IF_NOT(HasIntDataType(rPath))
            << "Data type is not int: " << rPath << std::endl;
    }
    else if (dtype_id == H5T_NATIVE_DOUBLE)
    {
        KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
            << "Data type is not float: " << rPath << std::endl;
    }

    hid_t file_id = GetFileId();
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    if (Mode == DataTransferMode::collective)
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    else
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
    hid_t dset_id = H5Dopen(file_id, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
    hid_t file_space_id = H5Dget_space(dset_id);
    hid_t mem_space_id = H5Screate_simple(ndims, local_mem_dims, nullptr);
    KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, local_start,
                                        nullptr, local_mem_dims, nullptr) < 0)
        << "H5Sselect_hyperslab failed." << std::endl;
    if (local_mem_dims[0] > 0)
    {
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id,
                                dxpl_id, &rData(0, 0)) < 0)
            << "H5Dread failed." << std::endl;
    }
    else
    {
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
#define KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(...)                                                                                      \
template void FileParallel::WriteDataSetImpl(const std::string&, const __VA_ARGS__&, DataTransferMode, WriteInfo& rInfo);
#endif

KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<int>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<double>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<array_1d<double, 3>>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<array_1d<double, 4>>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<array_1d<double, 6>>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Vector<array_1d<double, 9>>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Matrix<int>);
KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS(Matrix<double>);

template void FileParallel::ReadDataSetVectorImpl(const std::string& rPath,
                                                  Vector<int>& rData,
                                                  unsigned StartIndex,
                                                  unsigned BlockSize,
                                                  DataTransferMode Mode);
template void FileParallel::ReadDataSetVectorImpl(const std::string& rPath,
                                                  Vector<double>& rData,
                                                  unsigned StartIndex,
                                                  unsigned BlockSize,
                                                  DataTransferMode Mode);
template void FileParallel::ReadDataSetMatrixImpl(const std::string& rPath,
                                                  Matrix<int>& rData,
                                                  unsigned StartIndex,
                                                  unsigned BlockSize,
                                                  DataTransferMode Mode);
template void FileParallel::ReadDataSetMatrixImpl(const std::string& rPath,
                                                  Matrix<double>& rData,
                                                  unsigned StartIndex,
                                                  unsigned BlockSize,
                                                  DataTransferMode Mode);

#undef KRATOS_HDF5_FILE_PARRALLEL_INSTANTIATIONS

} // namespace HDF5.
} // namespace Kratos.