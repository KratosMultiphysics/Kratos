//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klauss B Sautter
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


namespace  Kratos
{

    /**
     * @class MoveMeshUtility
     *
     * @brief This utilitiy moves the nodes of a pfem fluid
     *        used especially in the pfem-fem coupling
     *
     * @author Klaus B Sautter
     */
    class  KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) MoveMeshUtility
    {
        public:
            typedef ModelPart::NodesContainerType NodesContainerType;

            KRATOS_CLASS_POINTER_DEFINITION(MoveMeshUtility);

            MoveMeshUtility() {};
            ~MoveMeshUtility() {};

            /**
             * @brief This function is the main operation of this utility.
             *        It moves the nodes of a given model part w.r.t. the current DISPLACEMENT
             * @param rNodes The nodes array of a model part
             */
            void MovePfemMesh(NodesContainerType& rNodes) const;

            /**
             * @brief This function resets the previous kinematic
             *        information of the PFEM nodes.
             * @param rFluidModelPart The fluid model part
             */
            void ResetPfemKinematicValues(ModelPart& rFluidModelPart);
        
        /*private:
            ModelPart& mrModelPart;*/
    }; // MoveMeshUtility
} //  Kratos

#endif //KRATOS_PFEM_MOVE_MESH_UTILITY_H_INCLUDED