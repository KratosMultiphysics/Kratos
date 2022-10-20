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

#if !defined(KRATOS_HDF5_FILE_PARALLEL_H_INCLUDED)
#define KRATOS_HDF5_FILE_PARALLEL_H_INCLUDED

// System includes

// External includes
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif

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

/// A class for accessing a single shared HDF5 file across MPI processes.
/**
 * This class is responsible for reading and writing data sets from all MPI
 * processes.
 */
class KRATOS_API(HDF5_APPLICATION) FileParallel : public File
{
    enum class DataTransferMode { independent, collective };
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FileParallel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit FileParallel(Parameters& rSettings);

    // Copy constructor.
    FileParallel(const FileParallel& rOther) = delete;

    /// Assignment operator.
    FileParallel& operator=(const FileParallel& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void WriteDataSet(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo) override;

    void WriteDataSet(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo) override;

    void WriteDataSet(const std::string& rPath, const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo) override;

    void WriteDataSet(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo) override;

    void WriteDataSet(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo) override;
    
    void WriteDataSetIndependent(const std::string& rPath, const Vector<int>& rData, WriteInfo& rInfo) override;
    
    void WriteDataSetIndependent(const std::string& rPath, const Vector<double>& rData, WriteInfo& rInfo) override;

    void WriteDataSetIndependent(const std::string& rPath,
                                const Vector<array_1d<double, 3>>& rData, WriteInfo& rInfo) override;

    void WriteDataSetIndependent(const std::string& rPath, const Matrix<int>& rData, WriteInfo& rInfo) override;
    
    void WriteDataSetIndependent(const std::string& rPath, const Matrix<double>& rData, WriteInfo& rInfo) override;
    
    unsigned GetPID() const override;

    unsigned GetTotalProcesses() const override;

    void ReadDataSet(const std::string& rPath,
                     Vector<int>& rData,
                     unsigned StartIndex, 
                     unsigned BlockSize) override;

    void ReadDataSet(const std::string& rPath,
                     Vector<double>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSet(const std::string& rPath,
                     Vector<array_1d<double, 3>>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSet(const std::string& rPath,
                     Matrix<int>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSet(const std::string& rPath,
                     Matrix<double>& rData,
                     unsigned StartIndex,
                     unsigned BlockSize) override;

    void ReadDataSetIndependent(const std::string& rPath,
                                Vector<int>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;

    void ReadDataSetIndependent(const std::string& rPath,
                                Vector<double>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;

    void ReadDataSetIndependent(const std::string& rPath,
                                Vector<array_1d<double, 3>>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;

    void ReadDataSetIndependent(const std::string& rPath,
                                Matrix<int>& rData,
                                unsigned StartIndex,
                                unsigned BlockSize) override;

    void ReadDataSetIndependent(const std::string& rPath,
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
    void WriteDataSetVectorImpl(const std::string& rPath,
                                const Vector<T>& rData,
                                DataTransferMode Mode,
                                WriteInfo& rInfo);

    template <class T>
    void WriteDataSetMatrixImpl(const std::string& rPath,
                                const Matrix<T>& rData,
                                DataTransferMode Mode,
                                WriteInfo& rInfo);

    template <class T>
    void ReadDataSetVectorImpl(const std::string& rPath,
                               Vector<T>& rData,
                               unsigned StartIndex,
                               unsigned BlockSize,
                               DataTransferMode Mode);

    template <class T>
    void ReadDataSetMatrixImpl(const std::string& rPath,
                               Matrix<T>& rData,
                               unsigned StartIndex,
                               unsigned BlockSize,
                               DataTransferMode Mode);
    ///@}
};

extern template void FileParallel::WriteDataSetVectorImpl(const std::string& rPath,
                                                        const Vector<int>& rData,
                                                        DataTransferMode Mode,
                                                        WriteInfo& rInfo);
extern template void FileParallel::WriteDataSetVectorImpl(const std::string& rPath,
                                                        const Vector<double>& rData,
                                                        DataTransferMode Mode,
                                                        WriteInfo& rInfo);
extern template void FileParallel::WriteDataSetMatrixImpl(const std::string& rPath,
                                                        const Matrix<int>& rData,
                                                        DataTransferMode Mode,
                                                        WriteInfo& rInfo);
extern template void FileParallel::WriteDataSetMatrixImpl(const std::string& rPath,
                                                        const Matrix<double>& rData,
                                                        DataTransferMode Mode,
                                                        WriteInfo& rInfo);
extern template void FileParallel::ReadDataSetVectorImpl(const std::string& rPath,
                                                       Vector<int>& rData,
                                                       unsigned StartIndex,
                                                       unsigned BlockSize,
                                                       DataTransferMode Mode);
extern template void FileParallel::ReadDataSetVectorImpl(const std::string& rPath,
                                                       Vector<double>& rData,
                                                       unsigned StartIndex,
                                                       unsigned BlockSize,
                                                       DataTransferMode Mode);
extern template void FileParallel::ReadDataSetMatrixImpl(const std::string& rPath,
                                                       Matrix<int>& rData,
                                                       unsigned StartIndex,
                                                       unsigned BlockSize,
                                                       DataTransferMode Mode);
extern template void FileParallel::ReadDataSetMatrixImpl(const std::string& rPath,
                                                       Matrix<double>& rData,
                                                       unsigned StartIndex,
                                                       unsigned BlockSize,
                                                       DataTransferMode Mode);

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_PARALLEL_H_INCLUDED defined