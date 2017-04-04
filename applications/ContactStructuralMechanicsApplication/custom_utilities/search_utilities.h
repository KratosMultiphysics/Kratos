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
     * This function fills the contact_container for the Mortar condition
     * @param ContactPoint: The destination point
     * @param Geom1 and Geom2: The geometries of the slave and master respectively
     * @param ContactNormal1 and ContactNormal2: The normals of the slave and master respectively
     * @param IntegrationOrder: The integration order   
     * @return contact_container: Once has been filled
     */

    static inline bool ContactContainerFiller(
            GeometryType& Geom1, // SLAVE
            GeometryType& Geom2, // MASTER
            const array_1d<double, 3> & ContactNormal1, // SLAVE
            const array_1d<double, 3> & ContactNormal2, // MASTER
            const double ActiveCheckFactor
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
                if (AuxDistance <= ActiveCheckFactor)
//                 if (AuxDistance <= ActiveCheckFactor && Geom2.IsInside(ProjectedPoint, Result)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
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
    
    static inline void ContactContainerFiller(
        std::vector<contact_container> *& ConditionPointers,
        Condition::Pointer & pCond1,       // SLAVE
        const Condition::Pointer & pCond2, // MASTER
        const array_1d<double, 3> & ContactNormal1, // SLAVE
        const array_1d<double, 3> & ContactNormal2, // MASTER
        const double ActiveCheckFactor
        )
    {
        const bool ConditionIsActive = ContactContainerFiller(pCond1->GetGeometry(), pCond2->GetGeometry(), ContactNormal1, ContactNormal2, ActiveCheckFactor);
        
        if (ConditionIsActive == true)
        {
            pCond1->Set(ACTIVE, true);
            contact_container AuxContactContainer;
            AuxContactContainer.condition   = pCond2;
            AuxContactContainer.active_pair = true;
            ConditionPointers->push_back(AuxContactContainer);
        }
    }
    
    static inline void ContactContainerFiller(
        ConditionMap *& ConditionPointers,
        Condition::Pointer & pCond1,       // SLAVE
        const Condition::Pointer & pCond2, // MASTER
        const array_1d<double, 3> & ContactNormal1, // SLAVE
        const array_1d<double, 3> & ContactNormal2, // MASTER
        const double ActiveCheckFactor
        )
    {
        const bool ConditionIsActive = ContactContainerFiller(pCond1->GetGeometry(), pCond2->GetGeometry(), ContactNormal1, ContactNormal2, ActiveCheckFactor);
        
        if (ConditionIsActive == true)
        {
            ConditionPointers->AddNewCondition(pCond2);
        }
    }
    
private:
};// class SearchUtilities
}
#endif /* KRATOS_SEARCH_UTILITIES defined */
