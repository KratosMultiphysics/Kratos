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

    explicit FileParallel(
        const DataCommunicator& rDataCommunicator,
        Parameters Settings);

    // Copy constructor.
    FileParallel(const FileParallel& rOther) = delete;

    /// Assignment operator.
    FileParallel& operator=(const FileParallel& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

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
    template <class TDataSetType>
    void WriteDataSetImpl(
        const std::string& rPath,
        const TDataSetType& rData,
        DataTransferMode Mode,
        WriteInfo& rInfo);

    template <class TDataSetType>
    void ReadDataSetImpl(
        const std::string& rPath,
        TDataSetType& rData,
        unsigned StartIndex,
        unsigned BlockSize,
        DataTransferMode Mode);

    ///@}
};

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_FILE_PARALLEL_H_INCLUDED defined