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

#if !defined(KRATOS_HDF5_SORTED_COORDINATES_PROCESS_H_INCLUDED)
#define KRATOS_HDF5_SORTED_COORDINATES_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"

// Application includes

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// Converts unsorted coordinates to sorted coordinates in HDF5.
/**
 * In partitioned simulations, the node ids are not sorted globally. This
 * function writes the (possibly) unsorted coordinates to an new sorted array
 * so that the connectivities can directly index the geometry without searching.
 * 
 * This class is used after running a Kratos simulation and before
 * post-processing with Xdmf. It expects global node ids to be in the range
 * [1, nnodes] and one-based indices to be used for the connectivities.
 * (see also ReorderConsecutiveModelPartIO).
 */
class SortedCoordinatesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(SortedCoordinatesProcess);
    ///@}
    ///@name Life Cycle
    ///@{
    SortedCoordinatesProcess(std::string FileName, std::string Prefix);

    ~SortedCoordinatesProcess() override
    {
    }
    ///@}
    ///@name Operations
    ///@{
    void Execute() override;
    ///@}
private:
    ///@name Member Variables
    ///@{
    std::string mFileName;
    std::string mPrefix;
    ///@}
}; // class SortedCoordinatesProcess

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_SORTED_COORDINATES_PROCESS_H_INCLUDED defined
