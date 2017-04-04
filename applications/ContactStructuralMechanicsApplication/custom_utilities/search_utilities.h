// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_SEARCH_UTILITIES)
#define KRATOS_SEARCH_UTILITIES

#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "geometries/point.h"
#include "custom_utilities/contact_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
class SearchUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This function checks if there is potential contact between two geometries (two conditions) 
     * @param Geom1: The geometry of the slave 
     * @param Geom2: The geometry of the master 
     * @param ContactNormal1: The normals of the slave
     * @param ContactNormal2: The normals of the master
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     * @return ConditionIsActive: True if the condition is active, false otherwise
     */

    static inline bool ContactChecker(
            GeometryType& Geom1, // SLAVE
            GeometryType& Geom2, // MASTER
            const array_1d<double, 3> & ContactNormal1, // SLAVE
            const array_1d<double, 3> & ContactNormal2, // MASTER
            const double ActiveCheckLength
            )
    {
        // Define the basic information
        const double Tolerance = std::numeric_limits<double>::epsilon();
        
        bool ConditionIsActive = false;
        
        for (unsigned int iNode = 0; iNode < Geom1.size(); iNode++)
        {
            if (Geom1[iNode].Is(ACTIVE) == false)
            {
                Point<3> ProjectedPoint;
                double AuxDistance = 0.0;
                const array_1d<double, 3> Normal = Geom1[iNode].GetValue(NORMAL);
                if (norm_2(Normal) < Tolerance)
                {
                    AuxDistance = ContactUtilities::FastProjectDirection(Geom2, Geom1[iNode], ProjectedPoint, ContactNormal2, ContactNormal1);
                }
                else
                {
                    AuxDistance = ContactUtilities::FastProjectDirection(Geom2, Geom1[iNode], ProjectedPoint, ContactNormal2, Normal);
                }  
                
//                 // Debug
//                 if (Geom1[iNode].GetValue(DISTANCE) > AuxDistance)
//                 {
//                     Geom1[iNode].GetValue(DISTANCE) = AuxDistance;
//                 }
                
                array_1d<double, 3> Result;
                if (AuxDistance <= ActiveCheckLength && Geom2.IsInside(ProjectedPoint, Result)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                { 
                    Geom1[iNode].Set(ACTIVE, true);
                    ConditionIsActive = true;
                }
             }
             else
             {
                 ConditionIsActive = true;
             }
         }
         
         return ConditionIsActive;
    }
    
    /**
     * This function fills the contact_container for the Mortar condition
     * @param ConditionPointers: The vector storing all the potential conditions
     * @param Geom1: The geometry of the slave 
     * @param Geom2: The geometry of the master 
     * @param ContactNormal1: The normals of the slave
     * @param ContactNormal2: The normals of the master
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     * @return ConditionIsActive: True if the condition is active, false otherwise
     */
    
    static inline void ContactContainerFiller(
        std::vector<contact_container> *& ConditionPointers,
        Condition::Pointer & pCond1,       // SLAVE
        const Condition::Pointer & pCond2, // MASTER
        const array_1d<double, 3> & ContactNormal1, // SLAVE
        const array_1d<double, 3> & ContactNormal2, // MASTER
        const double ActiveCheckLength
        )
    {
        const bool ConditionIsActive = ContactChecker(pCond1->GetGeometry(), pCond2->GetGeometry(), ContactNormal1, ContactNormal2, ActiveCheckLength);
        
        if (ConditionIsActive == true)
        {
            pCond1->Set(ACTIVE, true);
            contact_container AuxContactContainer;
            AuxContactContainer.condition   = pCond2;
            AuxContactContainer.active_pair = true;
            ConditionPointers->push_back(AuxContactContainer);
        }
    }
    
    /**
     * This function fills the ConditionMap for the Mortar condition
     * @param ConditionPointers: The map storing the potential contact conditions
     * @param Geom1: The geometry of the slave 
     * @param Geom2: The geometry of the master 
     * @param ContactNormal1: The normals of the slave
     * @param ContactNormal2: The normals of the master
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     * @return ConditionIsActive: True if the condition is active, false otherwise
     */
    
    static inline void ContactContainerFiller(
        ConditionMap *& ConditionPointers,
        Condition::Pointer & pCond1,       // SLAVE
        const Condition::Pointer & pCond2, // MASTER
        const array_1d<double, 3> & ContactNormal1, // SLAVE
        const array_1d<double, 3> & ContactNormal2, // MASTER
        const double ActiveCheckLength
        )
    {
        const bool ConditionIsActive = ContactChecker(pCond1->GetGeometry(), pCond2->GetGeometry(), ContactNormal1, ContactNormal2, ActiveCheckLength);
        
        if (ConditionIsActive == true)
        {
            ConditionPointers->AddNewCondition(pCond2);
        }
    }
    
private:
};// class SearchUtilities
}
#endif /* KRATOS_SEARCH_UTILITIES defined */
