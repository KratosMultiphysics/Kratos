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
    bool mDoPartitionedIO;
    std::string mPrefix;
    std::vector<std::string> mVariableNames;
    ///@}
    ///@name Private Operations
    ///@{
    std::tuple<unsigned, unsigned> GetStartIndexAndBlockSize() const;

    /// Divide nodes into local and ghost.
    void DivideNodes(NodesContainerType const& rNodes,
                     std::vector<NodeType*>& rLocalNodes,
                     std::vector<NodeType*>& rGhostNodes);

    void GetLocalNodes(NodesContainerType const& rNodes,
                       std::vector<NodeType*>& rLocalNodes);
    ///@}

}; // class NodalSolutionStepDataIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_NODAL_SOLUTION_STEP_DATA_IO_H_INCLUDED defined
