#include "voxel_mesher_operation_factory.h"

#include "includes/kratos_components.h"

#include "modeler/operation/voxel_mesher_operation.h"
#include "modeler/operation/find_contacts_in_skin_model_part.h"

namespace Kratos {

template class KratosComponents<VoxelMesherOperationFactory::RegisteredType>;

void RegisterVoxelMesherOperation()
{
    VoxelMesherOperationFactory::Register<FindContactsInSkinModelPart>("find_contacts_in_skin_model_part");
}

}