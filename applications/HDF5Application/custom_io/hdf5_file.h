//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_FILE_H_INCLUDED)
#define KRATOS_HDF5_FILE_H_INCLUDED

// System includes
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <regex>

// External includes
extern "C" {
#include "hdf5.h"
}
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@addtogroup HDF5Application
///@{

/// Provides helper functions that don't depend on the representation of a
/// particular class.
struct HDF5Utils
{
    /// Check if string is a valid path.
    /**
     * Valid paths are similar to linux file system with alphanumeric names
     * and possible underscores separated by '/'. All paths are absolute.
     */
    static bool IsPath(std::string Path)
    {
        return regex_match(Path, std::regex("(/\\w+)+"));
    }

    // Return vector of non-empty substrings separated by a delimiter.
    static std::vector<std::string> Split(std::string Path, char Delimiter)
    {
        std::vector<std::string> result;
        result.reserve(10);
        std::stringstream ss(Path);
        std::string sub_string;
        while (std::getline(ss, sub_string, Delimiter))
          if (sub_string.size() > 0)
            result.push_back(sub_string);

        return result;
    }
};

/// A base class for accessing an HDF5 file.
/**
 * Defines the interface to HDF5 files in Kratos. All functions reading or 
 * writing HDF5 meta data should be defined in this class. Functions reading
 * or writing data sets are defined by the derived class.
 */
class HDF5File
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HDF5File);

    template <class T>
    using Vector = boost::numeric::ublas::vector<T>;

    template <class T>
    using Matrix = boost::numeric::ublas::matrix<T>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit HDF5File(Parameters& rParams);

    // Copy constructor.
    HDF5File(const HDF5File& rOther) = delete;

    /// Destructor.
    virtual ~HDF5File();

    // Assignment operator.
    HDF5File& operator=(const HDF5File& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Check if path exists in HDF5 file.
    bool HasPath(std::string Path) const;

    bool IsGroup(std::string Path) const;

    bool IsDataSet(std::string Path) const;

    bool HasAttribute(std::string ObjectPath, std::string Name) const;

    void CreateGroup(std::string Path);

    void AddPath(std::string Path);

    template<class TScalar>
    void WriteAttribute(std::string ObjectPath, std::string Name, TScalar Value);

    template<class TScalar>
    void WriteAttribute(std::string ObjectPath, std::string Name, const Vector<TScalar>& rValue);

    template<class TScalar>
    void WriteAttribute(std::string ObjectPath, std::string Name, const Matrix<TScalar>& rValue);

    /// Write a data set to the HDF5 file.
    /**
     *  Performs collective write in MPI. The data is written blockwise according to
     *  processor rank.
     */
    virtual void WriteDataSet(std::string Path, const Vector<int>& rData);

    virtual void WriteDataSet(std::string Path, const Vector<double>& rData);

    virtual void WriteDataSet(std::string Path, const Vector<array_1d<double, 3>>& rData);

    virtual void WriteDataSet(std::string Path, const Matrix<int>& rData);
    
    virtual void WriteDataSet(std::string Path, const Matrix<double>& rData);

    /// Write the start and end indices of data blocks (by process rank).
    /**
     * Writes the partition array required to reconstruct a partitioned data set
     * from a file.
     */
    virtual void WriteDataPartition(std::string Path, const Vector<int>& rData);

    virtual void WriteDataPartition(std::string Path, const Vector<double>& rData);

    virtual void WriteDataPartition(std::string Path, const Vector<array_1d<double,3>>& rData);
    
    virtual void WriteDataPartition(std::string Path, const Matrix<int>& rData);
    
    virtual void WriteDataPartition(std::string Path, const Matrix<double>& rData);

    /// Independently write data set to the HDF5 file.
    /**
     * Performs independent write in MPI. Must be called collectively. Throws 
     * if more than one process has non-empty data.
     */
    virtual void WriteDataSetIndependent(std::string Path, const Vector<int>& rData);

    virtual void WriteDataSetIndependent(std::string Path, const Vector<double>& rData);

    virtual void WriteDataSetIndependent(std::string Path,
                                         const Vector<array_1d<double, 3>>& rData);

    virtual void WriteDataSetIndependent(std::string Path, const Matrix<int>& rData);

    virtual void WriteDataSetIndependent(std::string Path, const Matrix<double>& rData);

    std::vector<unsigned> GetDataDimensions(std::string Path) const;

    bool HasIntDataType(std::string Path) const;

    bool HasFloatDataType(std::string Path) const;

    void Flush();

    unsigned GetFileSize() const;

    std::string GetFileName() const;

    int GetEchoLevel() const;

    void SetEchoLevel(int Level);

    // Return this process Id with file access.
    virtual unsigned GetPID() const;

    // Return the total number of processes with file access.
    virtual unsigned GetTotalProcesses() const;

    /// Read a data set from the HDF5 file.
    /**
     * Performs collective read in MPI. Throws if out of range.
     */
    virtual void ReadDataSet(std::string Path,
                             Vector<int>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    virtual void ReadDataSet(std::string Path,
                             Vector<double>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    virtual void ReadDataSet(std::string Path,
                             Vector<array_1d<double, 3>>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    virtual void ReadDataSet(std::string Path,
                             Matrix<int>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);
   
    virtual void ReadDataSet(std::string Path,
                             Matrix<double>& rData,
                             unsigned StartIndex,
                             unsigned BlockSize);

    // Independently read data set from the HDF5 file.
    /**
     *  Performs independent read in MPI. Throws if out of range.
     */
    virtual void ReadDataSetIndependent(std::string Path,
                                       Vector<int>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    virtual void ReadDataSetIndependent(std::string Path,
                                       Vector<double>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    virtual void ReadDataSetIndependent(std::string Path,
                                       Vector<array_1d<double, 3>>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);

    virtual void ReadDataSetIndependent(std::string Path,
                                       Matrix<int>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);
 
    virtual void ReadDataSetIndependent(std::string Path,
                                       Matrix<double>& rData,
                                       unsigned StartIndex,
                                       unsigned BlockSize);
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

template <class TScalar>
void HDF5File::WriteAttribute(std::string ObjectPath, std::string Name, TScalar Value)
{
    KRATOS_TRY;
    hid_t type_id, space_id, attr_id;

    constexpr bool is_int_type = std::is_same<int, TScalar>::value;
    constexpr bool is_double_type = std::is_same<double, TScalar>::value;
    if (is_int_type)
        type_id = H5T_NATIVE_INT;
    else if (is_double_type)
        type_id = H5T_NATIVE_DOUBLE;
    else
        static_assert(is_int_type || is_double_type, "Unsupported data type.");

    space_id = H5Screate(H5S_SCALAR);
    KRATOS_ERROR_IF(space_id < 0) << "H5Screate failed." << std::endl;
    attr_id = H5Acreate_by_name(m_file_id, ObjectPath.c_str(), Name.c_str(), type_id,
                                space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Acreate_by_name failed." << std::endl;
    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, &Value) < 0) << "H5Awrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    KRATOS_CATCH("");
}

template <class TScalar>
void HDF5File::WriteAttribute(std::string ObjectPath, std::string Name, const Vector<TScalar>& rValue)
{
    KRATOS_TRY;
    hid_t type_id, space_id, attr_id;

    constexpr bool is_int_type = std::is_same<int, TScalar>::value;
    constexpr bool is_double_type = std::is_same<double, TScalar>::value;
    if (is_int_type)
        type_id = H5T_NATIVE_INT;
    else if (is_double_type)
        type_id = H5T_NATIVE_DOUBLE;
    else
        static_assert(is_int_type || is_double_type, "Unsupported data type.");

    const hsize_t dim = rValue.size();
    space_id = H5Screate_simple(1, &dim, nullptr);
    KRATOS_ERROR_IF(space_id < 0) << "H5Screate failed." << std::endl;
    attr_id = H5Acreate_by_name(m_file_id, ObjectPath.c_str(), Name.c_str(), type_id,
                                space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Acreate_by_name failed." << std::endl;
    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, &rValue[0]) < 0) << "H5Awrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    KRATOS_CATCH("");
}

template <class TScalar>
void HDF5File::WriteAttribute(std::string ObjectPath, std::string Name, const Matrix<TScalar>& rValue)
{
    KRATOS_TRY;
    hid_t type_id, space_id, attr_id;

    constexpr bool is_int_type = std::is_same<int, TScalar>::value;
    constexpr bool is_double_type = std::is_same<double, TScalar>::value;
    if (is_int_type)
        type_id = H5T_NATIVE_INT;
    else if (is_double_type)
        type_id = H5T_NATIVE_DOUBLE;
    else
        static_assert(is_int_type || is_double_type, "Unsupported data type.");

    const unsigned ndims = 2;
    hsize_t dims[ndims];
    dims[0] = rValue.size1();
    dims[1] = rValue.size2();
    space_id = H5Screate_simple(ndims, dims, nullptr);
    KRATOS_ERROR_IF(space_id < 0) << "H5Screate failed." << std::endl;
    attr_id = H5Acreate_by_name(m_file_id, ObjectPath.c_str(), Name.c_str(), type_id,
                                space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "H5Acreate_by_name failed." << std::endl;
    KRATOS_ERROR_IF(H5Awrite(attr_id, type_id, &rValue(0,0)) < 0) << "H5Awrite failed." << std::endl;
    KRATOS_ERROR_IF(H5Sclose(space_id) < 0) << "H5Sclose failed." << std::endl;
    KRATOS_ERROR_IF(H5Aclose(attr_id) < 0) << "H5Aclose failed." << std::endl;
    KRATOS_CATCH("");
}

///@} addtogroup
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_H_INCLUDED defined