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
//
//

#if !defined(KRATOS_MOVE_MESH_UTILITY_H_INCLUDED )
#define  KRATOS_MOVE_MESH_UTILITY_H_INCLUDED


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
     * @brief This utilitiy moves the nodes of a dem wall
     *        used especially in the fem-dem coupling
     *
     * @author Klaus B Sautter
     */

    class  KRATOS_API(DEM_APPLICATION) MoveMeshUtility
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
         * @param SetDeltaDisplacement A boolean to also update DELTA_DISPLACEMENT
         */

            void MoveDemMesh(NodesContainerType& rNodes, const bool SetDeltaDisplacement) const;
    }; // MoveMeshUtility
} //  Kratos

#endif //KRATOS_MOVE_MESH_UTILITY_H_INCLUDED