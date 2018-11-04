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

#if !defined(KRATOS_HDF5_FILE_H_INCLUDED)
#define KRATOS_HDF5_FILE_H_INCLUDED

// System includes
#include <vector>
#include <string>
#include <type_traits>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "utilities/builtin_timer.h"

// Application includes
#include "hdf5_application_define.h"

namespace Kratos
{

class Parameters;

namespace HDF5
{
///@addtogroup HDF5Application
///@{

namespace Internals
{
/// Check if string is a valid path.
/**
 * Valid paths are similar to linux file system with alphanumeric names
 * and possible underscores separated by '/'. All paths are absolute.
 */
bool IsPath(const std::string& rPath);

/// Return vector of non-empty substrings separated by a delimiter.
std::vector<std::string> Split(const std::string& rPath, char Delimiter);

template <class TScalar>
hid_t GetScalarDataType();

} // namespace Internals.

///@name Kratos Classes
///@{

/// Stores information about a data set written to HDF5.
struct WriteInfo
{
    std::size_t StartIndex = -1;
    std::size_t BlockSize = -1;
    std::size_t TotalSize = -1;
};

/// A base class for reading and writing an HDF5 file.
/**
 * This class stores the file id and is responsible for reading and writing
 * meta data. Reading and writing data sets is the responsibility of the derived
 * class.
 */
class File
{
public:
    ///@name Type Definitions
    ///@{

    template <class T>
    using Vector = HDF5::Vector<T>;

    template <class T>
    using Matrix = HDF5::Matrix<T>;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(File);

    ///@}
    ///@name Life Cycle
    ///@{

    explicit File(Parameters Settings);

    File(const File& rOther) = delete;

    File(File&& rOther);

    virtual ~File();

    File& operator=(const File& rOther) = delete;

    File& operator=(File&& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Check if path exists in HDF5 file.
    bool HasPath(const std::string& rPath) const;

    bool IsGroup(const std::string& rPath) const;

    bool IsDataSet(const std::string& rPath) const;

    bool HasAttribute(const std::string& rObjectPath, const std::string& rName) const;

    std::vector<std::string> GetAttributeNames(const std::string& rObjectPath) const;

    void CreateGroup(const std::string& rPath);

    std::vector<std::string> GetLinkNames(const std::string& rGroupPath) const;

    std::vector<std::string> GetGroupNames(const std::string& rGroupPath) const;

    void AddPath(const std::string& rPath);

    template<class TScalar>
    void WriteAttribute(const std::string& rObjectPath, const std::string& rName, TScalar Value);

    template<class TScalar>
    void WriteAttribute(const std::string& rObjectPath, const std::string& rName, const Vector<TScalar>& rValue);

    template<class TScalar>
    void WriteAttribute(const std::string& rObjectPath, const std::string& rName, const Matrix<TScalar>& rValue);

    void WriteAttribute(const std::string& rObjectPath, const std::string& rName, const std::string& rValue);

    void WriteAttribute(const std::string& rObjectPath, const std::string& rName, const array_1d<double, 3>& rValue);

    /// Write a data set to the HDF5 file.
    /**
     *  Performs collective write in MPI. The data is written blockwise according to
     *  processor rank.
     * 
     * @param[out] rInfo Information about the written data set.
     */
    virtual void WriteDataSet(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo);

    virtual void WriteDataSet(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo);

    virtual void WriteDataSet(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo);

    virtual void WriteDataSet(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo);
    
    virtual void WriteDataSet(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo);

    /// Independently write data set to the HDF5 file.
    /**
     * Performs independent write in MPI. Must be called collectively. Throws 
     * if more than one process has non-empty data.
     * 
     * @param[out] rInfo Information about the written data set.
     */
    virtual void WriteDataSetIndependent(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo);

    virtual void WriteDataSetIndependent(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo);

    virtual void WriteDataSetIndependent(const std::string& rPath,
                                         const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo);

    virtual void WriteDataSetIndependent(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo);

    virtual void WriteDataSetIndependent(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo);

    std::vector<unsigned> GetDataDimensions(const std::string& rPath) const;

    bool HasIntDataType(const std::string& rPath) const;

    bool HasFloatDataType(const std::string& rPath) const;

    void Flush();

    unsigned GetFileSize() const;

    std::string GetFileName() const;

    int GetEchoLevel() const;

    void SetEchoLevel(int Level);

    // Return this process Id with file access.
    virtual unsigned GetPID() const;

    // Return the total number of processes with file access.
    virtual unsigned GetTotalProcesses() const;

    template<class TScalar>
    void ReadAttribute(const std::string& rObjectPath, const std::string& rName, TScalar& rValue);

    template<class TScalar>
    void ReadAttribute(const std::string& rObjectPath, const std::string& rName, Vector<TScalar>& rValue);

    template<class TScalar>
    void ReadAttribute(const std::string& rObjectPath, const std::string& rName, Matrix<TScalar>& rValue);

    void ReadAttribute(const std::string& rObjectPath, const std::string& rName, std::string& rValue);

    void ReadAttribute(const std::string& rObjectPath, const std::string& rName, array_1d<double, 3>& rValue);

    /// Read a data set from the HDF5 file.
    /**
     * Performs collective read in MPI. Throws if out of range.
     */
    virtual void ReadDataSet(const std::string& rPath,
                             Vector<int>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    virtual void ReadDataSet(const std::string& rPath,
                             Vector<double>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    virtual void ReadDataSet(const std::string& rPath,
                             Vector<array_1d<double, 3>>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    virtual void ReadDataSet(const std::string& rPath,
                             Matrix<int>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);
   
    virtual void ReadDataSet(const std::string& rPath,
                             Matrix<double>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    // Independently read data set from the HDF5 file.
    /**
     *  Performs independent read in MPI. Throws if out of range.
     */
    virtual void ReadDataSetIndependent(const std::string& rPath,
                                       Vector<int>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    virtual void ReadDataSetIndependent(const std::string& rPath,
                                       Vector<double>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    virtual void ReadDataSetIndependent(const std::string& rPath,
                                       Vector<array_1d<double, 3>>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    virtual void ReadDataSetIndependent(const std::string& rPath,
                                       Matrix<int>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);
 
    virtual void ReadDataSetIndependent(const std::string& rPath,
                                       Matrix<double>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    unsigned GetOpenObjectsCount() const;
    ///@}

protected:
    ///@name Protected Operations
    ///@{
    hid_t GetFileId() const;
    ///@}

private:
    ///@name Member Variables
    ///@{
    std::string m_file_name;
    hid_t m_file_id = -1; // Default invalid file id.
    int m_echo_level = 0;
    ///@}

    ///@name Private Operations
    ///@{
    void SetFileDriver(const std::string& rDriver, hid_t FileAccessPropertyListId) const;
    ///@}

};
///@} // Kratos Classes

template <class TScalar>
void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, TScalar Value)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    hid_t type_id, space_id, attr_id;

    type_id = Internals::GetScalarDataType<TScalar>();
    space_id = H5Screate(H5S_SCALAR);
    KRATOS_ERROR_IF(space_id < 0) << "H5Screate failed." << std::endl;
    attr_id = H5Acreate_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), type_id,
                                space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Acreate_by_name failed." << std::endl;
    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, &Value) < 0) << "H5Awrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    if (GetEchoLevel() == 2 && GetPID() == 0)
        std::cout << "Write time \"" << rObjectPath << '/' << rName << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

template <class TScalar>
void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, const Vector<TScalar>& rValue)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    hid_t type_id, space_id, attr_id;

    type_id = Internals::GetScalarDataType<TScalar>();
    const hsize_t dim = rValue.size();
    space_id = H5Screate_simple(1, &dim, nullptr);
    KRATOS_ERROR_IF(space_id < 0) << "H5Screate failed." << std::endl;
    attr_id = H5Acreate_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), type_id,
                                space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Acreate_by_name failed." << std::endl;
    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, &rValue[0]) < 0) << "H5Awrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    if (GetEchoLevel() == 2 && GetPID() == 0)
        std::cout << "Write time \"" << rObjectPath << '/' << rName << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

template <class TScalar>
void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, const Matrix<TScalar>& rValue)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    hid_t type_id, space_id, attr_id;

    type_id = Internals::GetScalarDataType<TScalar>();
    const unsigned ndims = 2;
    hsize_t dims[ndims];
    dims[0] = rValue.size1();
    dims[1] = rValue.size2();
    space_id = H5Screate_simple(ndims, dims, nullptr);
    KRATOS_ERROR_IF(space_id < 0) << "H5Screate failed." << std::endl;
    attr_id = H5Acreate_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(), type_id,
                                space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Acreate_by_name failed." << std::endl;
    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, &rValue(0,0)) < 0) << "H5Awrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    if (GetEchoLevel() == 2 && GetPID() == 0)
        std::cout << "Write time \"" << rObjectPath << '/' << rName << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

template <class TScalar>
void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, TScalar& rValue)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    hid_t mem_type_id, attr_type_id, space_id, attr_id;
    int ndims;

    mem_type_id = Internals::GetScalarDataType<TScalar>();
    attr_id = H5Aopen_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(),
                                    H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Aopen_by_name failed." << std::endl;

    // Check data type.
    attr_type_id = H5Aget_type(attr_id);
    KRATOS_ERROR_IF(attr_type_id < 0) << "H5Aget_type failed." << std::endl;
    htri_t is_valid_type = H5Tequal(mem_type_id, attr_type_id);
    KRATOS_ERROR_IF(H5Tclose(attr_type_id) < 0) << "H5Tclose failed." << std::endl; 
    KRATOS_ERROR_IF(is_valid_type < 0) << "H5Tequal failed." << std::endl;
    KRATOS_ERROR_IF(is_valid_type == 0) << "Memory and file data types are different." << std::endl;

    // Check dimensions.
    space_id = H5Aget_space(attr_id);
    KRATOS_ERROR_IF(space_id < 0) << "H5Aget_space failed." << std::endl;
    KRATOS_ERROR_IF((ndims = H5Sget_simple_extent_ndims(space_id)) < 0)
        << "H5Sget_simple_extent_ndims failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(ndims != 0) << "Attribute \"" << rName << "\" is not scalar." << std::endl;
    
    // Read attribute.
    KRATOS_ERROR_IF(H5Aread(attr_id, mem_type_id, &rValue) < 0) << "H5Aread failed." << std::endl; 
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    if (GetEchoLevel() == 2 && GetPID() == 0)
        std::cout << "Read time \"" << rObjectPath << '/' << rName << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

template <class TScalar>
void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, Vector<TScalar>& rValue)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    hid_t mem_type_id, attr_type_id, space_id, attr_id;
    int ndims;
    hsize_t dims[1];

    mem_type_id = Internals::GetScalarDataType<TScalar>();
    attr_id = H5Aopen_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(),
                                    H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Aopen_by_name failed." << std::endl;

    // Check data type.
    attr_type_id = H5Aget_type(attr_id);
    KRATOS_ERROR_IF(attr_type_id < 0) << "H5Aget_type failed." << std::endl;
    htri_t is_valid_type = H5Tequal(mem_type_id, attr_type_id);
    KRATOS_ERROR_IF(H5Tclose(attr_type_id) < 0) << "H5Tclose failed." << std::endl; 
    KRATOS_ERROR_IF(is_valid_type < 0) << "H5Tequal failed." << std::endl;
    KRATOS_ERROR_IF(is_valid_type == 0) << "Memory and file data types are different." << std::endl;

    // Check dimensions.
    space_id = H5Aget_space(attr_id);
    KRATOS_ERROR_IF(space_id < 0) << "H5Aget_space failed." << std::endl;
    KRATOS_ERROR_IF((ndims = H5Sget_simple_extent_ndims(space_id)) < 0)
        << "H5Sget_simple_extent_ndims failed." << std::endl;
    KRATOS_ERROR_IF(ndims != 1) << "Attribute \"" << rName << "\" is not vector." << std::endl;
    KRATOS_ERROR_IF(H5Sget_simple_extent_dims(space_id, dims, nullptr) < 0)
        << "H5Sget_simple_extent_dims failed" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    rValue.resize(dims[0], false);
    // Read attribute.
    KRATOS_ERROR_IF(H5Aread(attr_id, mem_type_id, &rValue[0]) < 0) << "H5Aread failed." << std::endl; 
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    if (GetEchoLevel() == 2 && GetPID() == 0)
        std::cout << "Read time \"" << rObjectPath << '/' << rName << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

template <class TScalar>
void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, Matrix<TScalar>& rValue)
{
    KRATOS_TRY;
    BuiltinTimer timer;
    hid_t mem_type_id, attr_type_id, space_id, attr_id;
    int ndims;
    hsize_t dims[2];

    mem_type_id = Internals::GetScalarDataType<TScalar>();
    attr_id = H5Aopen_by_name(m_file_id, rObjectPath.c_str(), rName.c_str(),
                                    H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Aopen_by_name failed." << std::endl;

    // Check data type.
    attr_type_id = H5Aget_type(attr_id);
    KRATOS_ERROR_IF(attr_type_id < 0) << "H5Aget_type failed." << std::endl;
    htri_t is_valid_type = H5Tequal(mem_type_id, attr_type_id);
    KRATOS_ERROR_IF(H5Tclose(attr_type_id) < 0) << "H5Tclose failed." << std::endl; 
    KRATOS_ERROR_IF(is_valid_type < 0) << "H5Tequal failed." << std::endl;
    KRATOS_ERROR_IF(is_valid_type == 0) << "Memory and file data types are different." << std::endl;

    // Check dimensions.
    space_id = H5Aget_space(attr_id);
    KRATOS_ERROR_IF(space_id < 0) << "H5Aget_space failed." << std::endl;
    KRATOS_ERROR_IF((ndims = H5Sget_simple_extent_ndims(space_id)) < 0)
        << "H5Sget_simple_extent_ndims failed." << std::endl;
    KRATOS_ERROR_IF(ndims != 2) << "Attribute \"" << rName << "\" is not matrix." << std::endl;
    KRATOS_ERROR_IF(H5Sget_simple_extent_dims(space_id, dims, nullptr) < 0)
        << "H5Sget_simple_extent_dims failed" << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    rValue.resize(dims[0], dims[1], false);
    // Read attribute.
    KRATOS_ERROR_IF(H5Aread(attr_id, mem_type_id, &rValue(0,0)) < 0) << "H5Aread failed." << std::endl; 
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    if (GetEchoLevel() == 2 && GetPID() == 0)
        std::cout << "Read time \"" << rObjectPath << '/' << rName << "\": " << timer.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("Path: \"" + rObjectPath + '/' + rName + "\".");
}

namespace Internals
{
template <class TScalar>
hid_t GetScalarDataType()
{
    hid_t type_id;
    constexpr bool is_int_type = std::is_same<int, TScalar>::value;
    constexpr bool is_double_type = std::is_same<double, TScalar>::value;
    if (is_int_type)
        type_id = H5T_NATIVE_INT;
    else if (is_double_type)
        type_id = H5T_NATIVE_DOUBLE;
    else
        static_assert(is_int_type || is_double_type,
                      "Unsupported scalar data type.");

    return type_id;
}
} // namespace Internals.

///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_H_INCLUDED defined