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

    void CreateGroup(std::string Path);

    void AddPath(std::string Path);

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

///@} addtogroup
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_H_INCLUDED defined