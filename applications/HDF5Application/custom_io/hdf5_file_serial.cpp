#include "hdf5_file_serial.h"

#include <utility>
#include <vector>
#include "includes/kratos_parameters.h"
#include "utilities/builtin_timer.h"
#include "input_output/logger.h"

namespace
{
namespace hdf5_file_serial_cpp
{
using namespace Kratos;
using namespace HDF5;

hid_t CreateNewDataSet(const hid_t DataTypeId,
                       const std::vector<hsize_t>& rDims,
                       const hid_t FileId,
                       const std::string& rPath)
{
    const hid_t dspace_id = H5Screate_simple(rDims.size(), rDims.data(), nullptr);
    const hid_t dset_id = H5Dcreate(FileId, rPath.c_str(), DataTypeId, dspace_id,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(dspace_id) < 0) << "H5Sclose failed." << std::endl;
    return dset_id;
}

hid_t OpenExistingDataSet(const hid_t FileId,
                          const std::string& rPath)
{
    const hid_t dset_id = H5Dopen(FileId, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
    return dset_id;
}

bool HasMatchingScalarDataType(const HDF5::Vector<int>&, const File& rFile, const std::string& rPath)
{
    return rFile.HasIntDataType(rPath);
}

bool HasMatchingScalarDataType(const HDF5::Vector<double>&,
                               const File& rFile,
                               const std::string& rPath)
{
    return rFile.HasFloatDataType(rPath);
}

bool HasMatchingScalarDataType(const HDF5::Vector<array_1d<double, 3>>&,
                               const File& rFile,
                               const std::string& rPath)
{
    return rFile.HasFloatDataType(rPath);
}

bool HasMatchingScalarDataType(const HDF5::Matrix<int>&, const File& rFile, const std::string& rPath)
{
    return rFile.HasIntDataType(rPath);
}

bool HasMatchingScalarDataType(const HDF5::Matrix<double>&,
                               const File& rFile,
                               const std::string& rPath)
{
    return rFile.HasFloatDataType(rPath);
}

} // namespace hdf5_file_serial_cpp
} // namespace

namespace Kratos
{
namespace HDF5
{
FileSerial::FileSerial(Parameters& rSettings) : File(rSettings)
{
}

FileSerial::FileSerial(FileSerial&& rOther) : File(std::move(rOther))
{
}

FileSerial& FileSerial::operator=(FileSerial&& rOther)
{
    File::operator=(std::move(rOther));
    return *this;
}

void FileSerial::WriteDataSet(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSet(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetVectorImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

void FileSerial::WriteDataSetIndependent(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo)
{
    KRATOS_TRY;
    WriteDataSetMatrixImpl(rPath, rData, rInfo);
    KRATOS_CATCH("");
}

unsigned FileSerial::GetPID() const
{
    return 0;
}

unsigned FileSerial::GetTotalProcesses() const
{
    return 1;
}

void FileSerial::ReadDataSet(const std::string& rPath, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(const std::string& rPath, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSet(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Vector<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Vector<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Vector<array_1d<double, 3>>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetVectorImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Matrix<int>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void FileSerial::ReadDataSetIndependent(const std::string& rPath, Matrix<double>& rData, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    ReadDataSetMatrixImpl(rPath, rData, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

template <class T>
void FileSerial::WriteDataSetVectorImpl(const std::string& rPath,
                                        const Vector<T>& rData,
                                        WriteInfo& rInfo)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    using namespace hdf5_file_serial_cpp;

    // Create any missing subpaths.
    const auto pos = rPath.find_last_of('/');
    if (pos != 0) // Skip if last '/' is root.
    {
        const std::string sub_path = rPath.substr(0, pos);
        AddPath(sub_path);
    }

    // Create and write the data set.
    const hid_t dtype_id = Internals::GetScalarDataType(rData);
    const std::vector<hsize_t> dims = Internals::GetDataDimensions(rData);
    const hid_t file_id = GetFileId();
    hid_t dset_id{};
    if (!HasPath(rPath))
    {
        dset_id = CreateNewDataSet(dtype_id, dims, file_id, rPath);
    }
    else
    {
        KRATOS_ERROR_IF_NOT(HasMatchingScalarDataType(rData, *this, rPath))
            << "Wrong scalar data type: " << rPath << std::endl;
        KRATOS_ERROR_IF(Internals::GetDataDimensions(*this, rPath) != dims)
            << "Wrong dimensions: " << rPath << std::endl;
        dset_id = OpenExistingDataSet(file_id, rPath);
    }

    KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rData[0]) < 0)
        << "H5Dwrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    // Set the write info.
    rInfo.StartIndex = 0;
    rInfo.TotalSize = rInfo.BlockSize = rData.size();

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Write time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

template <class T>
void FileSerial::WriteDataSetMatrixImpl(const std::string& rPath,
                                        const Matrix<T>& rData,
                                        WriteInfo& rInfo)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    using namespace hdf5_file_serial_cpp;

    // Create any missing subpaths.
    const auto pos = rPath.find_last_of('/');
    if (pos != 0) // Skip if last '/' is root.
    {
        const std::string sub_path = rPath.substr(0, pos);
        AddPath(sub_path);
    }

    // Create and write the data set.
    const hid_t dtype_id = Internals::GetScalarDataType(rData);
    const std::vector<hsize_t> dims = Internals::GetDataDimensions(rData);
    const hid_t file_id = GetFileId();
    hid_t dset_id{};
    if (!HasPath(rPath))
    {
        dset_id = CreateNewDataSet(dtype_id, dims, file_id, rPath);
    }
    else
    {
        KRATOS_ERROR_IF_NOT(HasMatchingScalarDataType(rData, *this, rPath))
            << "Wrong scalar data type: " << rPath << std::endl;
        KRATOS_ERROR_IF(Internals::GetDataDimensions(*this, rPath) != dims)
            << "Wrong dimensions: " << rPath << std::endl;
        dset_id = OpenExistingDataSet(file_id, rPath);
    }
    KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rData(0, 0)) < 0)
        << "H5Dwrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;

    // Set the write info.
    rInfo.StartIndex = 0;
    rInfo.TotalSize = rInfo.BlockSize = rData.size1();

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Write time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

template <class T>
void FileSerial::ReadDataSetVectorImpl(const std::string& rPath,
                                       Vector<T>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    using namespace hdf5_file_serial_cpp;
    // Check that full path exists.
    KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
        << "Path is not a data set: " << rPath << std::endl;
    KRATOS_ERROR_IF_NOT(HasMatchingScalarDataType(rData, *this, rPath))
        << "Wrong scalar data type: " << rPath << std::endl;

    constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
    constexpr unsigned ndims = (!is_array_1d_type) ? 1 : 2;

    // Check consistency of file's data set dimensions.
    std::vector<unsigned> file_space_dims = GetDataDimensions(rPath);
    KRATOS_ERROR_IF(file_space_dims.size() != ndims)
        << "Invalid data set dimension." << std::endl;
    KRATOS_ERROR_IF(StartIndex + BlockSize > file_space_dims[0])
        << "StartIndex (" << StartIndex << ") + BlockSize (" << BlockSize
        << ") > size of data set (" << file_space_dims[0] << ")." << std::endl;
    if (is_array_1d_type)
        KRATOS_ERROR_IF(file_space_dims[1] != 3)
            << "Invalid data set dimension." << std::endl;

    if (rData.size() != BlockSize)
        rData.resize(BlockSize, false);

    // Initialize memory space dimensions.
    hsize_t start[ndims] = {0}, mem_dims[ndims];
    start[0] = StartIndex;
    mem_dims[0] = BlockSize; // Set first dimension.
    if (is_array_1d_type)
        mem_dims[1] = 3; // Set second dimension.

    hid_t file_id = GetFileId();
    hid_t dset_id = H5Dopen(file_id, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
    hid_t file_space_id = H5Dget_space(dset_id);
    hid_t mem_space_id = H5Screate_simple(ndims, mem_dims, nullptr);
    KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, nullptr, mem_dims, nullptr) < 0)
        << "H5Sselect_hyperslab failed." << std::endl;
    hid_t dtype_id = Internals::GetScalarDataType(rData);
    KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, H5P_DEFAULT, &rData[0]) < 0)
        << "H5Dread failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(file_space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(mem_space_id) < 0) << "H5Sclose failed." << std::endl;

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Read time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

template <class T>
void FileSerial::ReadDataSetMatrixImpl(const std::string& rPath,
                                       Matrix<T>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    using namespace hdf5_file_serial_cpp;
    // Check that full path exists.
    KRATOS_ERROR_IF_NOT(IsDataSet(rPath))
        << "Path is not a data set: " << rPath << std::endl;
    KRATOS_ERROR_IF_NOT(HasMatchingScalarDataType(rData, *this, rPath))
        << "Wrong scalar data type: " << rPath << std::endl;
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

    // Initialize memory space dimensions.
    hsize_t start[ndims] = {0}, mem_dims[ndims];
    start[0] = StartIndex;
    mem_dims[0] = rData.size1(); // Set first dimension.
    mem_dims[1] = rData.size2(); // Set second dimension.

    // Set the data type.
    hid_t file_id = GetFileId();
    hid_t dset_id = H5Dopen(file_id, rPath.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
    hid_t file_space_id = H5Dget_space(dset_id);
    hid_t mem_space_id = H5Screate_simple(ndims, mem_dims, nullptr);
    KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, nullptr, mem_dims, nullptr) < 0)
        << "H5Sselect_hyperslab failed." << std::endl;
    hid_t dtype_id = Internals::GetScalarDataType(rData);
    KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id,
                            H5P_DEFAULT, &rData(0, 0)) < 0)
        << "H5Dread failed." << std::endl;
    KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(file_space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(mem_space_id) < 0) << "H5Sclose failed." << std::endl;

    KRATOS_INFO_IF("HDF5Application", GetEchoLevel() == 2)
        << "Read time \"" << rPath << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rPath + "\".");
}

template void FileSerial::WriteDataSetVectorImpl(const std::string& rPath,
                                                 const Vector<int>& rData,
                                                 WriteInfo& rInfo);
template void FileSerial::WriteDataSetVectorImpl(const std::string& rPath,
                                                 const Vector<double>& rData,
                                                 WriteInfo& rInfo);
template void FileSerial::WriteDataSetMatrixImpl(const std::string& rPath,
                                                 const Matrix<int>& rData,
                                                 WriteInfo& rInfo);
template void FileSerial::WriteDataSetMatrixImpl(const std::string& rPath,
                                                 const Matrix<double>& rData,
                                                 WriteInfo& rInfo);
template void FileSerial::ReadDataSetVectorImpl(const std::string& rPath,
                                                Vector<int>& rData,
                                                unsigned StartIndex,
                                                unsigned BlockSize);
template void FileSerial::ReadDataSetVectorImpl(const std::string& rPath,
                                                Vector<double>& rData,
                                                unsigned StartIndex,
                                                unsigned BlockSize);
template void FileSerial::ReadDataSetMatrixImpl(const std::string& rPath,
                                                Matrix<int>& rData,
                                                unsigned StartIndex,
                                                unsigned BlockSize);
template void FileSerial::ReadDataSetMatrixImpl(const std::string& rPath,
                                                Matrix<double>& rData,
                                                unsigned StartIndex,
                                                unsigned BlockSize);

} // namespace HDF5.
} // namespace Kratos.