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
#include "includes/checks.h"




namespace  Kratos
{
    class KRATOS_API(DEM_APPLICATION) MoveMeshUtility
    {
        public:
            typedef ModelPart::NodesContainerType NodesContainerType;
            typedef ModelPart::ElementsContainerType ElementsContainerType;
            typedef ModelPart::ConditionsContainerType ConditionsContainerType;


            KRATOS_CLASS_POINTER_DEFINITION(MoveMeshUtility);

            MoveMeshUtility() {};
            ~MoveMeshUtility() {};

            void CalculateDeltaDispCustom(NodesContainerType& rNodes) const;
            void CalculateDeltaDispCustomFromIntermediatePos(NodesContainerType& rNodes) const;



            void MoveDemMesh(NodesContainerType& rNodes, const bool& rSetDeltaDisplacement) const;

            const bool CheckContact(NodesContainerType& rNodes) const;

            const bool CheckIsNearToWall(NodesContainerType& rNodes,Vector& rMaxDisplacement) const;

            const void SaveCurrentCoordinates(NodesContainerType& rNodes) const;

            const void ResetCoordinates(NodesContainerType& rNodes) const;

            //const void SetInactiveElements(ConditionsContainerType& rElementsDEM, const ElementsContainerType& rElementsFEM);

            const bool CompareTwoArrays(const array_1d<double, 3>& rA1,const array_1d<double, 3>& rA2) const
            {
                KRATOS_TRY;
                const double numerical_limit = std::numeric_limits<double>::epsilon();
                if ((std::abs(rA1[0]-rA2[0])<numerical_limit) && (std::abs(rA1[1]-rA2[1])<numerical_limit) && (std::abs(rA1[2]-rA2[2])<numerical_limit))
                {return true;}
                else return false;
                KRATOS_CATCH("")
            }

            std::vector<Condition*> mDEMConditions;
    }; // MoveMeshUtility
} //  Kratos

#endif //KRATOS_MOVE_MESH_UTILITY_H_INCLUDED