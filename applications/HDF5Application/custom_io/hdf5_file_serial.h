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
/**
 * This class is responsible for reading and writing data sets from a single
 * process.
 */
class KRATOS_API(HDF5_APPLICATION) FileSerial : public File
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FileSerial);

    ///@}
    ///@name Life Cycle
    ///@{

    explicit FileSerial(Parameters& rSettings);

    FileSerial(const FileSerial& rOther) = delete;

    FileSerial(FileSerial&& rOther);

    FileSerial& operator=(const FileSerial& rOther) = delete;

    FileSerial& operator=(FileSerial&& rOther);

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
    void WriteDataSetVectorImpl(const std::string& rPath, const Vector<T>& rData, WriteInfo& rInfo);

    template <class T>
    void WriteDataSetMatrixImpl(const std::string& rPath, const Matrix<T>& rData, WriteInfo& rInfo);

    template <class T>
    void ReadDataSetVectorImpl(const std::string& rPath,
                               Vector<T>& rData,
                               unsigned StartIndex,
                               unsigned BlockSize);

    template <class T>
    void ReadDataSetMatrixImpl(const std::string& rPath,
                               Matrix<T>& rData,
                               unsigned StartIndex,
                               unsigned BlockSize);
    ///@}
};

extern template void FileSerial::WriteDataSetVectorImpl(const std::string& rPath,
                                                        const Vector<int>& rData,
                                                        WriteInfo& rInfo);
extern template void FileSerial::WriteDataSetVectorImpl(const std::string& rPath,
                                                        const Vector<double>& rData,
                                                        WriteInfo& rInfo);
extern template void FileSerial::WriteDataSetMatrixImpl(const std::string& rPath,
                                                        const Matrix<int>& rData,
                                                        WriteInfo& rInfo);
extern template void FileSerial::WriteDataSetMatrixImpl(const std::string& rPath,
                                                        const Matrix<double>& rData,
                                                        WriteInfo& rInfo);
extern template void FileSerial::ReadDataSetVectorImpl(const std::string& rPath,
                                                       Vector<int>& rData,
                                                       unsigned StartIndex,
                                                       unsigned BlockSize);
extern template void FileSerial::ReadDataSetVectorImpl(const std::string& rPath,
                                                       Vector<double>& rData,
                                                       unsigned StartIndex,
                                                       unsigned BlockSize);
extern template void FileSerial::ReadDataSetMatrixImpl(const std::string& rPath,
                                                       Matrix<int>& rData,
                                                       unsigned StartIndex,
                                                       unsigned BlockSize);
extern template void FileSerial::ReadDataSetMatrixImpl(const std::string& rPath,
                                                       Matrix<double>& rData,
                                                       unsigned StartIndex,
                                                       unsigned BlockSize);

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_SERIAL_H_INCLUDED defined