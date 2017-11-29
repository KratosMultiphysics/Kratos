//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_FILE_SERIAL_H_INCLUDED)
#define KRATOS_HDF5_FILE_SERIAL_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// A class for accessing an HDF5 file from a single process.
class FileSerial : public File
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FileSerial);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit FileSerial(Parameters& rParams);

    // Copy constructor.
    FileSerial(const FileSerial& rOther) = delete;

    // Assignment operator.
    FileSerial& operator=(const FileSerial& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void WriteDataSet(std::string Path, const Vector<int>& rData) override;

    void WriteDataSet(std::string Path, const Vector<double>& rData) override;

    void WriteDataSet(std::string Path, const Vector<array_1d<double, 3>>& rData) override;

    void WriteDataSet(std::string Path, const Matrix<int>& rData) override;

    void WriteDataSet(std::string Path, const Matrix<double>& rData) override;

    void WriteDataPartition(std::string Path, const Vector<int>& rData) override;

    void WriteDataPartition(std::string Path, const Vector<double>& rData) override;

    void WriteDataPartition(std::string Path, const Vector<array_1d<double, 3>>& rData) override;

    void WriteDataPartition(std::string Path, const Matrix<int>& rData) override;

    void WriteDataPartition(std::string Path, const Matrix<double>& rData) override;

    void WriteDataSetIndependent(std::string Path, const Vector<int>& rData) override;

    void WriteDataSetIndependent(std::string Path, const Vector<double>& rData) override;

    void WriteDataSetIndependent(std::string Path,
                                 const Vector<array_1d<double, 3>>& rData) override;

    void WriteDataSetIndependent(std::string Path, const Matrix<int>& rData) override;

    void WriteDataSetIndependent(std::string Path, const Matrix<double>& rData) override;

    unsigned GetPID() const override;

    unsigned GetTotalProcesses() const override;

    void ReadDataSet(std::string Path,
                     Vector<int>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSet(std::string Path,
                     Vector<double>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSet(std::string Path,
                     Vector<array_1d<double, 3>>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSet(std::string Path,
                     Matrix<int>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSet(std::string Path,
                     Matrix<double>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSetIndependent(std::string Path,
                                Vector<int>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;

    void ReadDataSetIndependent(std::string Path,
                                Vector<double>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;

    void ReadDataSetIndependent(std::string Path,
                                Vector<array_1d<double, 3>>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;

    void ReadDataSetIndependent(std::string Path,
                                Matrix<int>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;

    void ReadDataSetIndependent(std::string Path,
                                Matrix<double>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;
    ///@}

protected:
    ///@name Protected Operations
    ///@{
    ///@}

private:
    ///@name Member Variables
    ///@{
    ///@}

    ///@name Private Operations
    ///@{
    template <class T>
    void WriteDataSetVectorImpl(std::string Path, const Vector<T>& rData)
    {
        KRATOS_TRY;
        boost::timer timer;
        // Expects a valid free path.
        KRATOS_ERROR_IF(HasPath(Path)) << "Path already exists: " << Path << std::endl;

        // Create any missing subpaths.
        auto pos = Path.find_last_of('/');
        if (pos != 0) // Skip if last '/' is root.
        {
            std::string sub_path = Path.substr(0, pos);
            AddPath(sub_path);
        }

        // Initialize data space dimensions.
        constexpr bool is_int_type = std::is_same<int, T>::value;
        constexpr bool is_double_type = std::is_same<double, T>::value;
        constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
        constexpr unsigned ndims = (!is_array_1d_type) ? 1 : 2;
        hsize_t dims[ndims] = {0};
        dims[0] = rData.size(); // Set first data space dimension.
        if (is_array_1d_type)
            dims[1] = 3; // Set second data space dimension (array_1d<double,3>).

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
        hid_t dspace_id = H5Screate_simple(ndims, dims, nullptr);
        hid_t file_id = GetFileId();
        hid_t dset_id = H5Dcreate(file_id, Path.c_str(), dtype_id, dspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;
        KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rData[0]) < 0)
            << "H5Dwrite failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(dspace_id) < 0) << "H5Sclose failed." << std::endl;

        if (GetEchoLevel() == 1)
            std::cout << "Write time \"" << Path << "\": " << timer.elapsed() << std::endl;
        KRATOS_CATCH("Path: \"" + Path + "\".");
    }

    template <class T>
    void WriteDataSetMatrixImpl(std::string Path, const Matrix<T>& rData)
    {
        KRATOS_TRY;
        boost::timer timer;
        // Check that full path does not exist before trying to write data.
        KRATOS_ERROR_IF(HasPath(Path)) << "Path already exists: " << Path << std::endl;

        // Create any missing subpaths.
        auto pos = Path.find_last_of('/');
        if (pos != 0) // Skip if last '/' is root.
        {
            std::string sub_path = Path.substr(0, pos);
            AddPath(sub_path);
        }

        // Initialize data space dimensions.
        const unsigned ndims = 2;
        hsize_t dims[ndims];
        dims[0] = rData.size1();
        dims[1] = rData.size2();

        // Set the data type.
        hid_t dtype_id = Internals::GetScalarDataType<T>();

        // Create and write the data set.
        hid_t dspace_id = H5Screate_simple(ndims, dims, nullptr);
        hid_t file_id = GetFileId();
        hid_t dset_id = H5Dcreate(file_id, Path.c_str(), dtype_id, dspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dcreate failed." << std::endl;
        KRATOS_ERROR_IF(H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rData(0, 0)) < 0)
            << "H5Dwrite failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(dspace_id) < 0) << "H5Sclose failed." << std::endl;

        if (GetEchoLevel() == 1)
            std::cout << "Write time \"" << Path << "\": " << timer.elapsed() << std::endl;
        KRATOS_CATCH("Path: \"" + Path + "\".");
    }

    template <class T>
    void WriteDataPartitionVectorImpl(std::string Path, const Vector<T>& rData)
    {
        int n = rData.size();
        Vector<int> partition(2); // Serial partition.
        partition[0] = 0;
        partition[1] = n;
       WriteDataSet(Path, partition);
    }

    template <class T>
    void WriteDataPartitionMatrixImpl(std::string Path, const Matrix<T>& rData)
    {
        int n = rData.size1();
        Vector<int> partition(2); // Serial partition.
        partition[0] = 0;
        partition[1] = n;
       WriteDataSet(Path, partition);
    }

    template <class T>
    void ReadDataSetVectorImpl(std::string Path, Vector<T>& rData, unsigned StartIndex, unsigned BlockSize)
    {
        KRATOS_TRY;
        boost::timer timer;
        // Check that full path exists.
        KRATOS_ERROR_IF_NOT(IsDataSet(Path))
            << "Path is not a data set: " << Path << std::endl;

        constexpr bool is_int_type = std::is_same<int, T>::value;
        constexpr bool is_double_type = std::is_same<double, T>::value;
        constexpr bool is_array_1d_type = std::is_same<array_1d<double, 3>, T>::value;
        constexpr unsigned ndims = (!is_array_1d_type) ? 1 : 2;

        // Check consistency of file's data set dimensions.
        std::vector<unsigned> file_space_dims = GetDataDimensions(Path);
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

        // Set the data type.
        hid_t dtype_id;
        if (is_int_type)
        {
            KRATOS_ERROR_IF_NOT(HasIntDataType(Path))
                << "Data type is not int: " << Path << std::endl;
            dtype_id = H5T_NATIVE_INT;
        }
        else if (is_double_type)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(Path))
                << "Data type is not float: " << Path << std::endl;
            dtype_id = H5T_NATIVE_DOUBLE;
        }
        else if (is_array_1d_type)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(Path))
                << "Data type is not float: " << Path << std::endl;
            dtype_id = H5T_NATIVE_DOUBLE;
        }
        else
            static_assert(is_int_type || is_double_type || is_array_1d_type,
                          "Unsupported data type.");

        hid_t file_id = GetFileId();
        hid_t dset_id = H5Dopen(file_id, Path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t file_space_id = H5Dget_space(dset_id);
        hid_t mem_space_id =
            H5Screate_simple(ndims, mem_dims, nullptr);
        KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, nullptr, mem_dims, nullptr) < 0)
            << "H5Sselect_hyperslab failed." << std::endl;
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, H5P_DEFAULT, &rData[0]) < 0)
            << "H5Dread failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(file_space_id) < 0) << "H5Sclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(mem_space_id) < 0) << "H5Sclose failed." << std::endl;

        if (GetEchoLevel() == 1)
            std::cout << "Read time \"" << Path << "\": " << timer.elapsed() << std::endl;
        KRATOS_CATCH("Path: \"" + Path + "\".");
    }

    template <class T>
    void ReadDataSetMatrixImpl(std::string Path, Matrix<T>& rData, unsigned StartIndex, unsigned BlockSize)
    {
        KRATOS_TRY;
        boost::timer timer;
        // Check that full path exists.
        KRATOS_ERROR_IF_NOT(IsDataSet(Path))
            << "Path is not a data set: " << Path << std::endl;

        const unsigned ndims = 2;

        // Check consistency of file's data set dimensions.
        std::vector<unsigned> file_space_dims = GetDataDimensions(Path);
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
        hid_t dtype_id = Internals::GetScalarDataType<T>();
        if (dtype_id == H5T_NATIVE_INT)
        {
            KRATOS_ERROR_IF_NOT(HasIntDataType(Path))
                << "Data type is not int: " << Path << std::endl;
        }
        else if (dtype_id == H5T_NATIVE_DOUBLE)
        {
            KRATOS_ERROR_IF_NOT(HasFloatDataType(Path))
                << "Data type is not float: " << Path << std::endl;
        }

        hid_t file_id = GetFileId();
        hid_t dset_id = H5Dopen(file_id, Path.c_str(), H5P_DEFAULT);
        KRATOS_ERROR_IF(dset_id < 0) << "H5Dopen failed." << std::endl;
        hid_t file_space_id = H5Dget_space(dset_id);
        hid_t mem_space_id =
            H5Screate_simple(ndims, mem_dims, nullptr);
        KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, nullptr, mem_dims, nullptr) < 0)
            << "H5Sselect_hyperslab failed." << std::endl;
        KRATOS_ERROR_IF(H5Dread(dset_id, dtype_id, mem_space_id, file_space_id, H5P_DEFAULT, &rData(0,0)) < 0)
            << "H5Dread failed." << std::endl;
        KRATOS_ERROR_IF(H5Dclose(dset_id) < 0) << "H5Dclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(file_space_id) < 0) << "H5Sclose failed." << std::endl;
        KRATOS_ERROR_IF(H5Sclose(mem_space_id) < 0) << "H5Sclose failed." << std::endl;

        if (GetEchoLevel() == 1)
            std::cout << "Read time \"" << Path << "\": " << timer.elapsed() << std::endl;
        KRATOS_CATCH("Path: \"" + Path + "\".");
    }
    ///@}

};

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_SERIAL_H_INCLUDED defined