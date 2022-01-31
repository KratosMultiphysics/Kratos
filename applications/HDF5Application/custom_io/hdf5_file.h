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

// External includes

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"

// Application includes
#include "hdf5_application_define.h"

namespace Kratos
{

class Parameters;

namespace HDF5
{
///@addtogroup HDF5Application
///@{

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

    void DeleteAttribute(const std::string& rObjectPath, const std::string& rName);

    std::vector<std::string> GetAttributeNames(const std::string& rObjectPath) const;

    void CreateGroup(const std::string& rPath);

    std::vector<std::string> GetLinkNames(const std::string& rGroupPath) const;

    std::vector<std::string> GetGroupNames(const std::string& rGroupPath) const;

    std::vector<std::string> GetDataSetNames(const std::string& rGroupPath) const;

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

extern template void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, int Value);
extern template void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, double Value);
extern template void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, const Vector<int>& rValue);
extern template void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, const Vector<double>& rValue);
extern template void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, const Matrix<int>& rValue);
extern template void File::WriteAttribute(const std::string& rObjectPath, const std::string& rName, const Matrix<double>& rValue);

extern template void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, int& rValue);
extern template void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, double& rValue);
extern template void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, Vector<int>& rValue);
extern template void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, Vector<double>& rValue);
extern template void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, Matrix<int>& rValue);
extern template void File::ReadAttribute(const std::string& rObjectPath, const std::string& rName, Matrix<double>& rValue);

///@} // Kratos Classes

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

hid_t GetScalarDataType(const Vector<int>&);
hid_t GetScalarDataType(const Vector<double>&);
hid_t GetScalarDataType(const Vector<array_1d<double, 3>>&);
hid_t GetScalarDataType(const Matrix<int>&);
hid_t GetScalarDataType(const Matrix<double>&);
hid_t GetScalarDataType(int);
hid_t GetScalarDataType(double);

std::vector<hsize_t> GetDataDimensions(const Vector<int>& rData);
std::vector<hsize_t> GetDataDimensions(const Vector<double>& rData);
std::vector<hsize_t> GetDataDimensions(const Vector<array_1d<double, 3>>& rData);
std::vector<hsize_t> GetDataDimensions(const Matrix<int>& rData);
std::vector<hsize_t> GetDataDimensions(const Matrix<double>& rData);
std::vector<hsize_t> GetDataDimensions(const File& rFile, const std::string& rPath);
} // namespace Internals.

///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_H_INCLUDED defined