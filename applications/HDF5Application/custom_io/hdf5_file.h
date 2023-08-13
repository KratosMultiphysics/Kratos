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

#pragma once

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "hdf5_application_define.h"

namespace Kratos
{

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
class KRATOS_API(HDF5_APPLICATION) File
{
private:
    ///@name Private enums
    ///@{

    enum class DataTransferMode
    {
        Independent,
        Collective
    };

    ///@}

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

    KRATOS_DEPRECATED_MESSAGE("Please use the constructor with DataCommunicator and Parameters")
    explicit File(Parameters Settings);

    explicit File(
        const DataCommunicator& rDataCommunicator,
        Parameters Settings);

    File(const File& rOther) = delete;

    File(File&& rOther);

    virtual ~File();

    File& operator=(const File& rOther) = delete;

    File& operator=(File&& rOther);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Checks if the given @ref rPath exists.
     *
     * @throws If a string with invalid characters are given.
     *
     * @param rPath     Path to be checked for.
     * @return true     If @ref rPath exists.
     * @return false    If @ref rPath does not exist.
     */
    bool HasPath(const std::string& rPath) const;

    /**
     * @brief Checks if @ref rPath contains a group.
     *
     * @param rPath     Path to be checked for.
     * @return true     If @ref rPath exists and it is a group.
     * @return false    Either @ref rPath does not exist or it is not a group.
     */
    bool IsGroup(const std::string& rPath) const;

    /**
     * @brief Checks if @ref rPath is a data set.
     *
     * @param rPath     Path to be checked for.
     * @return true     If @ref rPath exists and it is a data set.
     * @return false    Either @ref rPath does not exists or it is not a data set.
     */
    bool IsDataSet(const std::string& rPath) const;

    /**
     * @brief Checks if the given @ref rName attribute exists in @ref rObjectPath.
     *
     * @param rObjectPath       Path of a dataset or a group.
     * @param rName             Attribute name.
     * @return true             If attribute exists in the dataset or group.
     * @return false            If attribute does not exist in the dataset or group.
     */
    bool HasAttribute(
        const std::string& rObjectPath,
        const std::string& rName) const;

    /**
     * @brief Deletes the specified attribute from the dataset or group.
     *
     * @param rObjectPath       Path of a dataset or group.
     * @param rName             Attribute name to be deleted.
     */
    void DeleteAttribute(
        const std::string& rObjectPath,
        const std::string& rName);

    /**
     * @brief Get the Attribute Names list.
     *
     * @param rObjectPath                   Path of a dataset or group.
     * @return std::vector<std::string>     List of attribute names.
     */
    std::vector<std::string> GetAttributeNames(const std::string& rObjectPath) const;

    /**
     * @brief Create a Group.
     *
     * @param rPath     Path of the group to be created.
     */
    void CreateGroup(const std::string& rPath);

    /**
     * @brief Get the link names under a group.
     *
     * @param rGroupPath                    Group path.
     * @return std::vector<std::string>     List of link names.
     */
    std::vector<std::string> GetLinkNames(const std::string& rGroupPath) const;

    /**
     * @brief Get the sub group names under a group.
     *
     * @param rGroupPath                    Group path.
     * @return std::vector<std::string>     List of sub group names.
     */
    std::vector<std::string> GetGroupNames(const std::string& rGroupPath) const;

    /**
     * @brief Get the data set names under a group.
     *
     * @param rGroupPath                    Group path.
     * @return std::vector<std::string>     List of sub data set names.
     */
    std::vector<std::string> GetDataSetNames(const std::string& rGroupPath) const;

    /**
     * @brief Add group to the path recursively.
     *
     * This method creates all the parent groups to create the final @ref rPath.
     *
     * @throws If any of the parents in the @ref rPath is not a group.
     *
     * @param rPath     Path to be created.
     */
    void AddPath(const std::string& rPath);

    /**
     * @brief Write attributes to specified group or data set.
     *
     * @tparam TDataType            Data type of the attribute.
     * @param rObjectPath           Group or dataset path.
     * @param rName                 Attribute name.
     * @param rValue                Atribute data.
     */
    template<class TDataType>
    void WriteAttribute(
        const std::string& rObjectPath,
        const std::string& rName,
        const TDataType& rValue);

    /**
     * @brief Write data set to the hDF5 file.
     *
     * Performs collective write in MPI. The data is written blockwise according to
     * processor rank.
     *
     * @tparam TDataType    Data type of the Data.
     * @param rPath         Path to which the data is written.
     * @param rData         Data to be writtent.
     * @param rInfo         Information about the written data (output).
     */
    template<class TDataType>
    void WriteDataSet(
        const std::string& rPath,
        const TDataType& rData,
        WriteInfo& rInfo);

    /**
     * @brief Independently write data set to the HDF5 file.
     *
     * Performs independent write in MPI. Must be called collectively. Throws
     * if more than one process has non-empty data.
     *
     * @tparam TDataType    Data type of the Data.
     * @param rPath         Path to which the data is written.
     * @param rData         Data to be writtent.
     * @param rInfo         Information about the written data (output).
     */
    template<class TDataType>
    void WriteDataSetIndependent(
        const std::string& rPath,
        const TDataType& rData,
        WriteInfo& rInfo);

    /**
     * @brief Get the Data Dimensions at the specified path.
     *
     * @param rPath                     Dataset path.
     * @return std::vector<unsigned>    Dimensions of the data set.
     */
    std::vector<unsigned> GetDataDimensions(const std::string& rPath) const;

    /**
     * @brief Checks if the dataset at path is of int type.
     *
     * @param rPath                     Dataset path.
     * @return true                     If dataset of int type.
     * @return false                    If dataset is not of int type.
     */
    bool HasIntDataType(const std::string& rPath) const;

    /**
     * @brief Checks if the dataset at path is of double type.
     *
     * @param rPath                     Dataset path.
     * @return true                     If dataset of double type.
     * @return false                    If dataset is not of double type.
     */
    bool HasFloatDataType(const std::string& rPath) const;

    /**
     * @brief Checks the datast at path is of @ref TDataType.
     *
     * @tparam TDataType                Datatype to check against.
     * @param rPath                     Dataset path.
     * @return true                     If dataset is of TDataType.
     * @return false                    If dataset is not of TDataType.
     */
    template<class TDataType>
    bool HasDataType(const std::string& rPath) const;

    /**
     * @brief Flush the content to HDF5 file.
     *
     */
    void Flush();

    /**
     * @brief Terminate access to the HDF5 file.
     *
     * @throws If the underlying HDF5 call fails.
     *
     * @return * void
     */
    void Close();

    /**
     * @brief Get the File Size.
     *
     * @return unsigned
     */
    unsigned GetFileSize() const;

    /**
     * @brief Get the File Name.
     *
     * @return std::string
     */
    std::string GetFileName() const;

    /**
     * @brief Get the Echo Level.
     *
     * @return int
     */
    int GetEchoLevel() const;

    /**
     * @brief Set the Echo Level.
     *
     * @param Level
     */
    void SetEchoLevel(int Level);

    /**
     * @brief Get the Data Communicator.
     *
     * @return const DataCommunicator&
     */
    const DataCommunicator& GetDataCommunicator() const;

    /**
     * @brief Get the process rank.
     *
     * @return unsigned
     */
    unsigned GetPID() const;

    /**
     * @brief Get the Total number of processes with file access.
     *
     * @return unsigned
     */
    unsigned GetTotalProcesses() const;

    /**
     * @brief Read attribute from a group or dataset.
     *
     * @tparam TDataType            Datatype of the attribute.
     * @param rObjectPath           Group or dataset path.
     * @param rName                 Attribute name.
     * @param rValue                Attribute value to be written to.
     */
    template<class TDataType>
    void ReadAttribute(
        const std::string& rObjectPath,
        const std::string& rName,
        TDataType& rValue);

    /**
     * @brief Read a data set from the HDF5 file.
     *
     * Performs collective read in MPI. Throws if out of range.
     *
     * @tparam TDataType        Datatype of the read data.
     * @param rPath             Path of the dataset.
     * @param rData             Data to be written to.
     * @param StartIndex        Starting offset of data for this rank.
     * @param BlockSize         Number of data points for this rank.
     */
    template<class TDataType>
    void ReadDataSet(
        const std::string& rPath,
        TDataType& rData,
        const unsigned StartIndex,
        const unsigned BlockSize);

    /**
     * @brief Independently read a data set from the HDF5 file.
     *
     * Performs independent read in MPI. Throws if out of range.
     *
     * @tparam TDataType        Datatype of the read data.
     * @param rPath             Path of the dataset.
     * @param rData             Data to be written to.
     * @param StartIndex        Starting offset of data for this rank.
     * @param BlockSize         Number of data points for this rank.
     */
    template<class TDataType>
    void ReadDataSetIndependent(
        const std::string& rPath,
        TDataType& rData,
        const unsigned StartIndex,
        const unsigned BlockSize);

    unsigned GetOpenObjectsCount() const;
    ///@}

private:
    ///@name Member Variables
    ///@{

    DataCommunicator const * mpDataCommunicator;

    std::string m_file_name;

    hid_t m_file_id = -1; // Default invalid file id.

    int m_echo_level = 0;

    ///@}

    ///@name Private Operations
    ///@{

    void SetFileDriver(
        const std::string& rDriver,
        hid_t FileAccessPropertyListId) const;

    hid_t GetFileId() const;

    template<class TDataType>
    void GetDataSet(
        hid_t& rDataSetId,
        hid_t& rDataSpaceId,
        const std::vector<hsize_t>& rDims,
        const std::string& rPath);

    void CreateNewDataSet(
        hid_t& rDataSetId,
        hid_t& rDataSpaceId,
        const hid_t DataTypeId,
        const std::vector<hsize_t>& rDims,
        const std::string& rPath);

    hid_t OpenExistingDataSet(const std::string& rPath);

    template<class TDataType, DataTransferMode TDataTransferMode>
    void WriteDataSetImpl(
        const std::string& rPath,
        const TDataType& rData,
        WriteInfo& rInfo);

    template<class TDataType, DataTransferMode TDataTransferMode>
    void ReadDataSetImpl(
        const std::string& rPath,
        TDataType& rData,
        const unsigned StartIndex,
        const unsigned BlockSize);

    ///@}
};

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