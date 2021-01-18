//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    https://github.com/msandre
//

/** @file local_ghost_splitting_utility.h
 *  @brief Methods for splitting nodes into local and ghost containers.
 */

#if !defined(KRATOS_LOCAL_GHOST_SPLITTING_UTILITY_H_INCLUDED)
#define KRATOS_LOCAL_GHOST_SPLITTING_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{

void SplitNodesIntoLocalAndGhost(ModelPart::NodesContainerType const& rNodes,
                                 std::vector<ModelPart::NodeType*>& rLocalNodes,
                                 std::vector<ModelPart::NodeType*>& rGhostNodes);

void GetLocalNodes(ModelPart::NodesContainerType const& rNodes,
                   std::vector<ModelPart::NodeType*>& rLocalNodes);

} // namespace Kratos.

#endif // KRATOS_LOCAL_GHOST_SPLITTING_UTILITY_H_INCLUDED defined
