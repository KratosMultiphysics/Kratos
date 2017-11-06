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

#if !defined(KRATOS_HDF5_NODAL_SOLUTION_STEP_DATA_IO_H_INCLUDED)
#define KRATOS_HDF5_NODAL_SOLUTION_STEP_DATA_IO_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/communicator.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{

///@name Kratos Classes
///@{

/// A class for IO of nodal solution step data in HDF5.
class NodalSolutionStepDataIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NodalSolutionStepDataIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NodalSolutionStepDataIO(Parameters& rParams, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    std::string GetPrefix() const;

    void SetPrefix(std::string const& rPrefix);

    void WriteNodalResults(NodesContainerType const& rNodes, unsigned Step=0);

    void ReadNodalResults(NodesContainerType& rNodes, Communicator& rComm, unsigned Step=0);
    
    ///@}

private:
    ///@name Member Variables
    ///@{
    File::Pointer mpFile;
    std::string mPrefix;
    std::vector<std::string> mVariableNames;
    ///@}
    ///@name Private Operations
    ///@{
    std::tuple<unsigned, unsigned> GetStartIndexAndBlockSize() const;

    template <class TVariableType, class TFileDataType>
    void SetDataBuffer(TVariableType const& rVariable,
                       std::vector<NodeType*> const& rNodes,
                       Vector<TFileDataType>& rData,
                       unsigned Step);

    template <class TVariableType, class TFileDataType>
    void SetNodalSolutionStepData(TVariableType const& rVariable,
                                  Vector<TFileDataType> const& rData,
                                  std::vector<NodeType*>& rNodes,
                                  unsigned Step);
    ///@}

}; // class NodalSolutionStepDataIO.

template <class TVariableType, class TFileDataType>
void NodalSolutionStepDataIO::SetDataBuffer(TVariableType const& rVariable,
                                            std::vector<NodeType*> const& rNodes,
                                            Vector<TFileDataType>& rData,
                                            unsigned Step)
{
    KRATOS_TRY;

    rData.resize(rNodes.size(), false);
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rNodes.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
            rData[i] = rNodes[i]->FastGetSolutionStepValue(rVariable, Step);
    }

    KRATOS_CATCH("");
}

template <class TVariableType, class TFileDataType>
void NodalSolutionStepDataIO::SetNodalSolutionStepData(TVariableType const& rVariable,
                                                       Vector<TFileDataType> const& rData,
                                                       std::vector<NodeType*>& rNodes,
                                                       unsigned Step)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rData.size() != rNodes.size())
        << "File data block size (" << rData.size()
        << ") is not equal to number of nodes (" << rNodes.size() << ")." << std::endl;

    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rNodes.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
            rNodes[i]->FastGetSolutionStepValue(rVariable, Step) = rData[i];
    }

    KRATOS_CATCH("");
}

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_NODAL_SOLUTION_STEP_DATA_IO_H_INCLUDED defined
