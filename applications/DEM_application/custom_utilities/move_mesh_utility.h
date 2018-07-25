//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klauss B Sautter (based on Riccardo Rossi's work)
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
    class KRATOS_API(DEM_APPLICATION) MoveMeshUtility
    {
        public:
            typedef ModelPart::NodesContainerType NodesContainerType;

            KRATOS_CLASS_POINTER_DEFINITION(MoveMeshUtility);

            MoveMeshUtility() {};
            ~MoveMeshUtility() {};

            void MoveDemMesh(NodesContainerType& rNodes, const bool& rSetDeltaDisplacement) const;

            const bool CheckContact(NodesContainerType& rNodes) const;

            const bool CheckIsNearToWall(NodesContainerType& rNodes,Vector& rMaxDisplacement) const;

    }; // MoveMeshUtility
} //  Kratos

#endif //KRATOS_MOVE_MESH_UTILITY_H_INCLUDED