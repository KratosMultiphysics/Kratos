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

/// Stores information about a dataset written to HDF5.
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
     * @brief Checks if the given @a rPath exists.
     *
     * @throws If a string with invalid characters is given.
     *
     * @param rPath     Path to be checked.
     * @return true     True if @a rPath exists, otherwise false.
     */
    bool HasPath(const std::string& rPath) const;

    /**
     * @brief Checks if @a rPath contains a group.
     *
     * @param rPath     Path to be checked for.
     * @return true     True if @a rPath exists and it is a group., otherwise false.
     */
    bool IsGroup(const std::string& rPath) const;

    /**
     * @brief Checks if @a rPath is a dataset.
     *
     * @param rPath     Path to be checked for.
     * @return true     True if @a rPath exists and it is a dataset, otherwise false.
     */
    bool IsDataSet(const std::string& rPath) const;

    /**
     * @brief Checks if the given @a rName attribute exists in @a rObjectPath.
     *
     * @param rObjectPath       Path of a dataset or a group.
     * @param rName             Attribute name.
     * @return true             True if attribute exists in the dataset or group, otherwise false.
     */
    bool HasAttribute(
        const std::string& rObjectPath,
        const std::string& rName) const;

    /**
     * @brief Checks if the attribute is of the specified type.
     *
     * @tparam TDataType        Data type to check.
     * @param rObjectPath       Dataset or group path.
     * @param rName             Attribute name.
     * @return true             If the attribute has data of type @a TDataType.
     * @return false            If the attribute does not has the data of type @a TDataType.
     */
    template<class TDataType>
    bool HasAttributeType(
        const std::string& rObjectPath,
        const std::string& rName) const;

    /**
     * @brief Get the Attribute Dimensions.
     *
     * @param rObjectPath               Dataset or group path.
     * @param rName                     Attribute name.
     * @return std::vector<hsize_t>     Size in each dimension.
     */
    std::vector<hsize_t> GetAttributeDimensions(
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
     * @brief Get the Attributes' Names list.
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
     * @brief Get the dataset names under a group.
     *
     * @param rGroupPath                    Group path.
     * @return std::vector<std::string>     List of sub dataset names.
     */
    std::vector<std::string> GetDataSetNames(const std::string& rGroupPath) const;

    /**
     * @brief Add a group to the path recursively.
     *
     * This method creates all the parent groups to create the final @a rPath.
     *
     * @throws If any of the parents in the @a rPath is not a group.
     *
     * @param rPath     Path to be created.
     */
    void AddPath(const std::string& rPath);

    /**
     * @brief Write attributes to specified group or dataset.
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
     * @brief Write attributes in a @ref Parameters object to dataset or group.
     *
     * @param rObjectPath           Dataset or group path.
     * @param Attributes            Attributes to be written to.
     */
    void WriteAttribute(
        const std::string& rObjectPath,
        const Parameters Attributes);

    /**
     * @brief Write a dataset to the hDF5 file.
     *
     * Performs collective write in MPI. The data is written blockwise according to
     * processor rank.
     *
     * @tparam TDataType    Data type of the provided data.
     * @param rPath         Path to which the data is written.
     * @param rData         Data to be written.
     * @param rInfo         Information about the written data (output).
     */
    template<class TDataType>
    void WriteDataSet(
        const std::string& rPath,
        const TDataType& rData,
        WriteInfo& rInfo);

    /**
     * @brief Independently write dataset to the HDF5 file.
     *
     * Performs independent write in MPI. Must be called collectively. Throws
     * if more than one process has non-empty data.
     *
     * @tparam TDataType    Data type of the Data.
     * @param rPath         Path to which the data is written.
     * @param rData         Data to be written.
     * @param rInfo         Information about the written data (output).
     */
    template<class TDataType>
    void WriteDataSetIndependent(
        const std::string& rPath,
        const TDataType& rData,
        WriteInfo& rInfo);

    /**
     * @brief Get the shape of the dataset stored at the specified path.
     *
     * @param rPath                     Dataset path.
     * @return std::vector<unsigned>    Dimensions of the dataset.
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
     * @brief Checks the datast at path is of @a TDataType.
     *
     * @tparam TDataType                Datatype to check against.
     * @param rPath                     Dataset path.
     * @return true                     If dataset is of TDataType.
     * @return false                    If dataset is not of TDataType.
     */
    template<class TDataType>
    bool HasDataType(const std::string& rPath) const;

    /// @brief Flush the content to HDF5 file.
    void Flush();

    /**
     * @brief Terminate access to the HDF5 file.
     *
     * @throws If the underlying HDF5 call fails.
     *
     * @return * void
     */
    void Close();

    unsigned GetFileSize() const;

    std::string GetFileName() const;

    int GetEchoLevel() const;

    void SetEchoLevel(int Level);

    const DataCommunicator& GetDataCommunicator() const;

    /// @brief Get the process rank.
    unsigned GetPID() const;

    /// @brief Get the Total number of processes with file access.
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
        TDataType& rValue) const;

    /**
     * @brief Read attributes from a dataset or group to a @ref Parameters object.
     *
     * @param rObjectPath           Dataset or group path.
     * @return Parameters           Parameters object containing attribute name, value pairs.
     */
    Parameters ReadAttribute(const std::string& rObjectPath) const;

    /**
     * @brief Read a dataset from the HDF5 file.
     *
     * Performs collective read in MPI. Throws if out of range.
     *
     * @throws if out of range.
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
        const unsigned BlockSize) const;

    /**
     * @brief Independently read a dataset from the HDF5 file.
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
        const unsigned BlockSize) const;

    unsigned GetOpenObjectsCount() const;
    ///@}

private:
    ///@name Member Variables
    ///@{

    DataCommunicator const * mpDataCommunicator;

    std::string mFileName;

    hid_t mFileId = -1; // Default invalid file id.

    int mEchoLevel = 0;

    ///@}

    ///@name Private Operations
    ///@{

    void SetFileDriver(
        const std::string& rDriver,
        hid_t FileAccessPropertyListId) const;

    hid_t GetFileId() const;

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
        const unsigned BlockSize) const;

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
bool KRATOS_API(HDF5_APPLICATION) IsPath(const std::string& rPath);

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