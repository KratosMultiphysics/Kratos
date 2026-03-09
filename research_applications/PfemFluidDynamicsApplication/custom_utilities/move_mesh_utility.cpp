//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B Sautter
//                   Alejandro Cornejo Velazquez
//                   Carlos Eulogio Flores
//

// System includes

// External includes

// Project includes
#include "move_mesh_utility.h"
#include "pfem_fluid_dynamics_application_variables.h"
#include "utilities/parallel_utilities.h"

namespace  Kratos {

void PFEMMoveMeshUtility::ResetPfemKinematicValues(ModelPart& rFluidModelPart)
{
    block_for_each(rFluidModelPart.Nodes(), [&](NodeType& rNode){
        if (rNode.IsNot(RIGID) && rNode.IsNot(SOLID)) { // We update only the fluid part
            array_1d<double, 3>& r_current_disp = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
            array_1d<double, 3>& r_current_vel  = rNode.FastGetSolutionStepValue(VELOCITY, 0);
            array_1d<double, 3>& r_current_acc  = rNode.FastGetSolutionStepValue(ACCELERATION, 0);

            const array_1d<double, 3> old_disp = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);

            noalias(rNode.Coordinates()) = rNode.GetInitialPosition() + old_disp;

            noalias(r_current_disp) = old_disp;
            noalias(r_current_vel)  = rNode.FastGetSolutionStepValue(VELOCITY, 1);
            noalias(r_current_acc)  = rNode.FastGetSolutionStepValue(ACCELERATION, 1);
        }
    });
}

} //  Kratos
