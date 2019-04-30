#include "hdf5_file_parallel.h"

#include "includes/kratos_parameters.h"
#include "utilities/builtin_timer.h"
#include "input_output/logger.h"

namespace Kratos
{
namespace HDF5
{
FileParallel::FileParallel(Parameters& rSettings) : File(rSettings)
{
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSet(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, DataTransferMode::collective, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

void FileParallel::WriteDataSetIndependent(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, DataTransferMode::independent, rInfo);
    KRATOS_CATCH("");
}

unsigned FileParallel::GetPID() const
{
    int rank, ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Comm_rank failed." << std::endl;
    return static_cast<unsigned>(rank);
}

unsigned FileParallel::GetTotalProcesses() const
{
    int num_proc, ierr;
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Comm_size failed." << std::endl;
    return static_cast<unsigned>(num_proc);
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

template <class T>
void FileParallel::WriteDataSetVectorImpl(const std::string& rPath,
                                          const Vector<T>& rData,
                                          DataTransferMode Mode,
                                          WriteInfo& rInfo)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    // Expects a valid free path.
    KRATOS_ERROR_IF(HasPath(rPath)) << "Path already exists: " << rPath << std::endl;

    // Create any missing subpaths.
    auto pos = rPath.find_last_of('/');
    if (pos != 0) // Skip if last '/' is root.
    {
        std::string sub_path = rPath.substr(0, pos);
        AddPath(sub_path);
    }

    // Initialize data space dimensions.
    constexpr bool is_int_type = std::is_same<int, T>::value;
    constexpr bool is_double_type = std::is_same<double, T>::value;
    constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
    constexpr unsigned ndims = (!is_array_1d_type) ? 1 : 2;
    hsize_t local_dims[ndims], global_dims[ndims];
    // Set first data space dimension.
    local_dims[0] = rData.size();
    unsigned send_buf, recv_buf;
    send_buf = local_dims[0];
    MPI_Allreduce(&send_buf, &recv_buf, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    global_dims[0] = recv_buf;
    if (Mode == DataTransferMode::independent)
        if (local_dims[0] > 0)
            KRATOS_ERROR_IF(local_dims[0] != global_dims[0])
                << "Can't perform independent write with MPI. Invalid data."
                << std::endl;
    if (is_array_1d_type)
        local_dims[1] = global_dims[1] = 3; // Set second data space dimension.

    hsize_t local_start[ndims];
    if (Mode == DataTransferMode::collective)
    {
        send_buf = local_dims[0];
        // Get global position where local data set ends.
        MPI_Scan(&send_buf, &recv_buf, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        // Set global position where local data set starts.
        local_start[0] = recv_buf - local_dims[0];
    }
    else
        local_start[0] = 0;

    if (is_array_1d_type)
        local_start[1] = 0;

    // Set the data type.
    hid_t dtype_id;
    if (is_int_type)
        dtype_id = H5T_NATIVE_INT;
    else if (is_double_type)
        dtype_id = H5T_NATIVE_DOUBLE;
    else if (is_array_1d_type)
        dtype_id = H5T_NATIVE_DOUBLE;
    else
        static_assert(is_int_type || is_double_type || is_array_1d_type,
                      "Unsupported data type.");

    // Create and write the data set.
    hid_t file_id = GetFileId();
    hid_t fspace_id = H5Screate_simple(ndims, global_dims, nullptr);
    // H5Dcreate() must be called collectively for both collective and
    // independent write.
    hid_t dset_id = H5Dcreate(file_id, rPath.c_str(), dtype_id, fspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;

    if (Mode == DataTransferMode::collective || local_dims[0] > 0)
    {
        hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
        if (Mode == DataTransferMode::collective)
            H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
        else
            H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
        H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, local_start, nullptr,
                            local_dims, nullptr);
        hid_t mspace_id = H5Screate_simple(ndims, local_dims, nullptr);
        if (local_dims[0] > 0)
        {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, &rData[0]) < 0)
                << "H5Dwrite failed." << std::endl;
        }
        else
        {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, nullptr) < 0) << "H5Dwrite failed. Please ensure global data set is non-empty."
                                                                                                     << std::endl;
        }
        KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(mspace_id) < 0) << "H5Sclose failed." << std::endl;
    }
    KRATOS_ERROR_IF(H5Sclose(fspace_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    // Set the write info.
    rInfo.StartIndex = local_start[0];
    rInfo.BlockSize = local_dims[0];
    rInfo.TotalSize = global_dims[0];

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Write time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

template <class T>
void FileParallel::WriteDataSetMatrixImpl(const std::string& rPath,
                                          const Matrix<T>& rData,
                                          DataTransferMode Mode,
                                          WriteInfo& rInfo)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    // Check that full path does not exist before trying to write data.
    KRATOS_ERROR_IF(HasPath(rPath)) << "Path already exists: " << rPath << std::endl;

    // Create any missing subpaths.
    auto pos = rPath.find_last_of('/');
    if (pos != 0) // Skip if last '/' is root.
    {
        std::string sub_path = rPath.substr(0, pos);
        AddPath(sub_path);
    }

    // Initialize data space dimensions.
    const unsigned ndims = 2;
    hsize_t local_dims[ndims], global_dims[ndims];
    // Set first data space dimension.
    local_dims[0] = rData.size1();
    unsigned send_buf, recv_buf;
    send_buf = local_dims[0];
    MPI_Allreduce(&send_buf, &recv_buf, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    global_dims[0] = recv_buf;
    if (Mode == DataTransferMode::independent)
        if (local_dims[0] > 0)
            KRATOS_ERROR_IF(local_dims[0] != global_dims[0])
                << "Can't perform independent write with MPI. Invalid data."
                << std::endl;
    local_dims[1] = global_dims[1] = rData.size2(); // Set second data space dimension.

    hsize_t local_start[ndims] = {0};
    if (Mode == DataTransferMode::collective)
    {
        send_buf = local_dims[0];
        // Get global position where local data set ends.
        MPI_Scan(&send_buf, &recv_buf, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        // Set global position where local data set starts.
        local_start[0] = recv_buf - local_dims[0];
    }
    else
        local_start[0] = 0;

    // Set the data type.
    hid_t dtype_id = Internals::GetScalarDataType(rData);

    // Create and write the data set.
    hid_t file_id = GetFileId();
    hid_t fspace_id = H5Screate_simple(ndims, global_dims, nullptr);
    // H5Dcreate() must be called collectively for both collective and
    // independent write.
    hid_t dset_id = H5Dcreate(file_id, rPath.c_str(), dtype_id, fspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;

    if (Mode == DataTransferMode::collective || local_dims[0] > 0)
    {
        hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
        if (Mode == DataTransferMode::collective)
            H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
        else
            H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
        H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, local_start, nullptr,
                            local_dims, nullptr);
        hid_t mspace_id = H5Screate_simple(ndims, local_dims, nullptr);
        if (local_dims[0] > 0)
        {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id,
                                     dxpl_id, &rData(0, 0)) < 0)
                << "H5Dwrite failed." << std::endl;
        }
        else
        {
            KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, mspace_id, fspace_id, dxpl_id, nullptr) < 0) << "H5Dwrite failed. Please ensure global data set is non-empty"
                                                                                                     << " and the second matrix dimension is consistent across all processes."
                                                                                                     << std::endl;
        }
        KRATOS_ERROR_IF(H5Pclose(dxpl_id) < 0) << "H5Pclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(mspace_id) < 0) << "H5Sclose failed." << std::endl;
    }
    KRATOS_ERROR_IF(H5Sclose(fspace_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    // Set the write info.
    rInfo.StartIndex = local_start[0];
    rInfo.BlockSize = local_dims[0];
    rInfo.TotalSize = global_dims[0];

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Write time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

template <class T>
void FileParallel::ReadDataSetVectorImpl(const std::string& rPath,
                                         Vector<T>& rData,
                                         unsigned StartIndex,
                                         unsigned BlockSize,
                                         DataTransferMode Mode)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    // Check that full path exists.
    KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
        << "Path is not a data set: " << rPath << std::endl;

    constexpr bool is_int_type = std::is_same<int, T>::value;
    constexpr bool is_double_type = std::is_same<double, T>::value;
    constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
    constexpr unsigned ndims = (!is_array_1d_type) ? 1 : 2;

    // Check consistency of file's data set dimensions.
    std::vector<unsigned> file_space_dims = GetDataDimensions(rPath);
    KRATOS_ERROR_IF(file_space_dims.size() != ndims)
        << "Invalid data set dimension." << std::endl;
    KRATOS_ERROR_IF(StartIndex + BlockSize > file_space_dims[0])
        << "StartIndex (" << StartIndex << ") + BlockSize (" << BlockSize
        << ") > size of data set (" << file_space_dims[0] << ")." << std::endl;
    hsize_t local_mem_dims[ndims];
    // Set first memory space dimension.
    local_mem_dims[0] = BlockSize;
    if (is_array_1d_type)
        local_mem_dims[1] = 3; // Set second dimension.
    if (is_array_1d_type)
        KRATOS_ERROR_IF(file_space_dims[1] != 3)
            << "Invalid data set dimension." << std::endl;

    if (rData.size() != BlockSize)
        rData.resize(BlockSize, false);

    // Set global position where local data set starts.
    hsize_t local_start[ndims];
    local_start[0] = StartIndex;
    if (is_array_1d_type)
        local_start[1] = 0;

    // Set the data type.
    hid_t dtype_id;
    if (is_int_type)
    {
        KRATOS_ERROR_IF_NOT(HasIntDataType(rPath))
            << "Data type is not int: " << rPath << std::endl;
        dtype_id = H5T_NATIVE_INT;
    }
    else if (is_double_type)
    {
        KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
            << "Data type is not float: " << rPath << std::endl;
        dtype_id = H5T_NATIVE_DOUBLE;
    }
    else if (is_array_1d_type)
    {
        KRATOS_ERROR_IF_NOT(HasFloatDataType(rPath))
            << "Data type is not float: " << rPath << std::endl;
        dtype_id = H5T_NATIVE_DOUBLE;
    }
    else
        static_assert(is_int_type || is_double_type || is_array_1d_type,
                      "Unsupported data type.");

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
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, dxpl_id, &rData[0]) < 0)
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

template void FileParallel::WriteDataSetVectorImpl(const std::string& rPath,
                                                   const Vector<int>& rData,
                                                   DataTransferMode Mode,
                                                   WriteInfo& rInfo);
template void FileParallel::WriteDataSetVectorImpl(const std::string& rPath,
                                                   const Vector<double>& rData,
                                                   DataTransferMode Mode,
                                                   WriteInfo& rInfo);
template void FileParallel::WriteDataSetMatrixImpl(const std::string& rPath,
                                                   const Matrix<int>& rData,
                                                   DataTransferMode Mode,
                                                   WriteInfo& rInfo);
template void FileParallel::WriteDataSetMatrixImpl(const std::string& rPath,
                                                   const Matrix<double>& rData,
                                                   DataTransferMode Mode,
                                                   WriteInfo& rInfo);
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

} // namespace HDF5.
} // namespace Kratos.