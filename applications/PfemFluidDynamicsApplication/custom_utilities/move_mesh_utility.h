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
//

#ifndef KRATOS_PFEM_MOVE_MESH_UTILITY_H_INCLUDED
#define KRATOS_PFEM_MOVE_MESH_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace  Kratos {
/**
 * @class PFEMMoveMeshUtility
 *
 * @brief This utility resets the nodes of a PFEM fluid
 *        ModelPart, used especially in the PFEM-FEM coupling
 *
 * @author Klaus B Sautter
 */

class  KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) PFEMMoveMeshUtility
{
public:
    typedef Node NodeType;

    KRATOS_CLASS_POINTER_DEFINITION(PFEMMoveMeshUtility);

    /**
     * @brief This function resets the previous kinematic
     *        information of the specified PFEM ModelPart nodes.
     * @param rFluidModelPart The fluid ModelPart
     */
    static void ResetPfemKinematicValues(ModelPart& rFluidModelPart);
}; // PFEMMoveMeshUtility

} //  Kratos

#endif //KRATOS_PFEM_MOVE_MESH_UTILITY_H_INCLUDED
