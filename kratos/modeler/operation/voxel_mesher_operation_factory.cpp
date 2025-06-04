//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "voxel_mesher_operation_factory.h"

#include "includes/kratos_components.h"

#include "modeler/operation/voxel_mesher_operation.h"
#include "modeler/operation/find_contacts_in_skin_model_part.h"
#include "modeler/operation/compute_surrogate_boundary_data.h"

namespace Kratos {

template class KratosComponents<VoxelMesherOperationFactory::RegisteredType>;

void RegisterVoxelMesherOperation()
{
    VoxelMesherOperationFactory::Register<FindContactsInSkinModelPart>("find_contacts_in_skin_model_part");
    VoxelMesherOperationFactory::Register<ComputeSurrogateBoundaryData>("compute_surrogate_boundary_data");
}

}