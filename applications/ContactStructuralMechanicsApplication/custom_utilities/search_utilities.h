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

/* Utilities */
#include "utilities/math_utils.h"
#include "custom_utilities/exact_mortar_segmentation_utility.h"
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
    
    typedef Point<3>                                  PointType;
    typedef Node<3>                                    NodeType;
    typedef Geometry<NodeType>                     GeometryType;
    typedef Geometry<PointType>               GeometryPointType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod   IntegrationMethod;
    
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
        // Initialize geometries
        GeometryType& Geom1 = pCond1->GetGeometry(); // SLAVE
        GeometryType& Geom2 = pCond2->GetGeometry(); // MASTER
        
//         // Initialize variables
//         bool ConditionIsActive = false;
//         const unsigned int Dimension  = Geom1.WorkingSpaceDimension();
//         const unsigned int NumberOfNodes = Geom1.size();
//         
//         if (Dimension == 2 && NumberOfNodes == 2)
//         {
//             ExactMortarIntegrationUtility<2, 2> IntUtil = ExactMortarIntegrationUtility<2, 2>();
//             std::vector<array_1d<PointType,2>>  ConditionsPointsSlave;
//             ConditionIsActive = IntUtil.GetExactIntegration(Geom1, ContactNormal1, Geom2, ContactNormal2, ConditionsPointsSlave);
//         }
//         else if (Dimension == 3 && NumberOfNodes == 3)
//         {
//             ExactMortarIntegrationUtility<3, 3> IntUtil = ExactMortarIntegrationUtility<3, 3>();
//             std::vector<array_1d<PointType,3>>  ConditionsPointsSlave;
//             ConditionIsActive = IntUtil.GetExactIntegration(Geom1, ContactNormal1, Geom2, ContactNormal2, ConditionsPointsSlave);
//         }
//         else if (Dimension == 3 && NumberOfNodes == 4)
//         {
//             ExactMortarIntegrationUtility<3, 4> IntUtil = ExactMortarIntegrationUtility<3, 4>();
//             std::vector<array_1d<PointType,3>>  ConditionsPointsSlave;
//             ConditionIsActive = IntUtil.GetExactIntegration(Geom1, ContactNormal1, Geom2, ContactNormal2, ConditionsPointsSlave);
//         }
//         else
//         {
//             KRATOS_ERROR << "INTEGRATION NOT IMPLEMENTED" << std::endl;
//         }
//         
//         if (ConditionIsActive == true)
//         {
//             for (unsigned int iNode = 0; iNode < Geom1.size(); iNode++)
//             {
//                 Geom1[iNode].Set(ACTIVE, true);
//             }
//         }
        
        // LEGACY WAY
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
                    
                    // Geom1[iNode].SetLock();
                    Geom1[iNode].Set(ACTIVE, true);
                    // Geom1[iNode].UnSetLock();
                }
                else if (DualCheck == true)
                {
                    AuxDistance = ContactUtilities::FastProjectDirection(Geom2, Geom1[iNode], ProjectedPoint, ContactNormal2, -ContactNormal2);
                    if (AuxDistance <= ActiveCheckLength && (StrictCheck == true ? Geom2.IsInside(ProjectedPoint, Result, Tolerance) : true)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                    {
                        ConditionIsActive = true;
                        
                        // Geom1[iNode].SetLock();
                        Geom1[iNode].Set(ACTIVE, true);
                        // Geom1[iNode].UnSetLock();
                    }
                }
             }
             else
             {
                 ConditionIsActive = true;
             }
         }
        
        // If condition is active we add
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
