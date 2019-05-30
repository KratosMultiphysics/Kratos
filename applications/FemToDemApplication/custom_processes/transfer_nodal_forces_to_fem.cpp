//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/transfer_nodal_forces_to_fem.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

TransferNodalForcesToFem::TransferNodalForcesToFem(
    ModelPart& rModelPart,
    ModelPart& rDemModelPart)
    : mrModelPart(rModelPart),
      mrDEMModelPart(rDemModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

void TransferNodalForcesToFem::Execute() 
{
	auto& sub_model_conditions = mrModelPart.GetSubModelPart("ContactForcesDEMConditions");

    const auto it_cond_begin = sub_model_conditions.ConditionsBegin();
    //#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(sub_model_conditions.Conditions().size()); i++) {
        auto it_cond = it_cond_begin + i;
        auto& r_geometry = it_cond->GetGeometry();
        const int id_node = r_geometry[0].Id();

        if (mrModelPart.GetNode(id_node).GetValue(IS_DEM)) {
            const array_1d<double,3>& dem_forces = mrDEMModelPart.GetNode(id_node).GetSolutionStepValue(TOTAL_FORCES);
            it_cond->SetValue(FORCE_LOAD, dem_forces);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos