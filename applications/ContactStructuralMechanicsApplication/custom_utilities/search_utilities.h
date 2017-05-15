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

// System includes

// External includes

// Project includes
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
            const array_1d<double, 3>& ContactNormal1, // SLAVE
            const array_1d<double, 3>& ContactNormal2, // MASTER
            const double ActiveCheckLength,
            const bool DualCheck = false,
            const bool StrictCheck = true
            )
    {
        // Define the basic information
//         const double Tolerance = 1.0e-12;
        const double Tolerance = std::numeric_limits<double>::epsilon();
        
        bool ConditionIsActive = false;
        
//         #pragma omp for
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
                
                array_1d<double, 3> Result;
                if (AuxDistance <= ActiveCheckLength && (StrictCheck == true ? Geom2.IsInside(ProjectedPoint, Result, Tolerance) : true)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                {
                    ConditionIsActive = true;
                    
//                     #pragma omp critical
                    Geom1[iNode].Set(ACTIVE, true);
                }
                else if (DualCheck == true)
                {
                    AuxDistance = ContactUtilities::FastProjectDirection(Geom2, Geom1[iNode], ProjectedPoint, ContactNormal2, -ContactNormal2);
                    if (AuxDistance <= ActiveCheckLength && (StrictCheck == true ? Geom2.IsInside(ProjectedPoint, Result, Tolerance) : true)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                    {
                        ConditionIsActive = true;
                        
    //                     #pragma omp critical
                        Geom1[iNode].Set(ACTIVE, true);
                    }
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
     * This function fills the ConditionSet for the Mortar condition
     * @param ConditionPointers: The map storing the potential contact conditions
     * @param Geom1: The geometry of the slave 
     * @param Geom2: The geometry of the master 
     * @param ContactNormal1: The normals of the slave
     * @param ContactNormal2: The normals of the master
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     * @return ConditionIsActive: True if the condition is active, false otherwise
     */
    
    static inline void ContactContainerFiller(
        boost::shared_ptr<ConditionSet>& ConditionPointers,
        Condition::Pointer & pCond1,       // SLAVE
        const Condition::Pointer & pCond2, // MASTER
        const array_1d<double, 3> & ContactNormal1, // SLAVE
        const array_1d<double, 3> & ContactNormal2, // MASTER
        const double ActiveCheckLength,
        const bool DualCheck = false, 
        const bool StrictCheck = true 
        )
    {
        const bool ConditionIsActive = ContactChecker(pCond1->GetGeometry(), pCond2->GetGeometry(), ContactNormal1, ContactNormal2, ActiveCheckLength, DualCheck, StrictCheck);
        
        if (ConditionIsActive == true)
        {
            ConditionPointers->AddNewCondition(pCond2);
        }
    }
    
protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{       

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};// class SearchUtilities

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_SEARCH_UTILITIES defined */
