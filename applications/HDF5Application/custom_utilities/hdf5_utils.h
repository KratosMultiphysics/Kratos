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

#if !defined(KRATOS_HDF5_UTILS_H_INCLUDED)
#define KRATOS_HDF5_UTILS_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "hdf5_application_define.h"

namespace Kratos
{
namespace HDF5
{
namespace Detail
{
/// Divide nodes into local and ghost.
void DivideNodes(NodesContainerType const& rNodes,
                 std::vector<NodeType*>& rLocalNodes,
                 std::vector<NodeType*>& rGhostNodes,
                 bool IsPartitioned);

void GetLocalNodes(NodesContainerType const& rNodes,
                   std::vector<NodeType*>& rLocalNodes,
                   bool IsPartitioned);

} // namespace Detail.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_UTILS_H_INCLUDED defined
