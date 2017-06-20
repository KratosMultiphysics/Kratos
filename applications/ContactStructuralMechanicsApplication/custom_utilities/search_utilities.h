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
     * @param geom_1: The geometry of the slave 
     * @param geom_2: The geometry of the master 
     * @param ContactNormal1: The normals of the slave
     * @param ContactNormal2: The normals of the master
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     * @return condition_is_active: True if the condition is active, false otherwise
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
        GeometryType& geom_1 = pCond1->GetGeometry(); // SLAVE
        GeometryType& geom_2 = pCond2->GetGeometry(); // MASTER
        
        // LEGACY WAY
        // Define the basic information
//         const double Tolerance = 1.0e-12;
        const double Tolerance = std::numeric_limits<double>::epsilon();
        
        bool condition_is_active = false;
//         #pragma omp for
        for (unsigned int i_node = 0; i_node < geom_1.size(); i_node++)
        {
            if (geom_1[i_node].Is(ACTIVE) == false)
            {
                Point<3> projected_point;
                double aux_distance = 0.0;
                const array_1d<double, 3> normal = geom_1[i_node].GetValue(NORMAL);
                if (norm_2(normal) < Tolerance)
                {
                    aux_distance = ContactUtilities::FastProjectDirection(geom_2, geom_1[i_node], projected_point, ContactNormal2, ContactNormal1);
                }
                else
                {
                    aux_distance = ContactUtilities::FastProjectDirection(geom_2, geom_1[i_node], projected_point, ContactNormal2, normal);
                }  
                
                array_1d<double, 3> result;
                if (aux_distance <= ActiveCheckLength && (StrictCheck == true ? geom_2.IsInside(projected_point, result, Tolerance) : true)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                {
                    condition_is_active = true;
                    
                    // geom_1[i_node].SetLock();
                    geom_1[i_node].Set(ACTIVE, true);
                    // geom_1[i_node].UnSetLock();
                }
                else if (DualCheck == true)
                {
                    aux_distance = ContactUtilities::FastProjectDirection(geom_2, geom_1[i_node], projected_point, ContactNormal2, -ContactNormal2);
                    if (aux_distance <= ActiveCheckLength && (StrictCheck == true ? geom_2.IsInside(projected_point, result, Tolerance) : true)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                    {
                        condition_is_active = true;
                        
                        // geom_1[i_node].SetLock();
                        geom_1[i_node].Set(ACTIVE, true);
                        // geom_1[i_node].UnSetLock();
                    }
                }
             }
             else
             {
                 condition_is_active = true;
             }
         }
        
        // If condition is active we add
        if (condition_is_active == true)
        {
            ConditionPointers->AddNewCondition(pCond2);
        }
        
//         // Initialize variables
//         bool condition_is_active = false;
//         const unsigned int dimension  = geom_1.WorkingSpaceDimension();
//         const unsigned int number_of_nodes = geom_1.size();
//         
//         if (dimension == 2 && number_of_nodes == 2)
//         {
//             ExactMortarIntegrationUtility<2, 2> IntUtil = ExactMortarIntegrationUtility<2, 2>();
//             std::vector<array_1d<PointType,2>> conditions_points_slave;
//             condition_is_active = IntUtil.GetExactIntegration(geom_1, ContactNormal1, geom_2, ContactNormal2, conditions_points_slave);
//             
//             if (condition_is_active == true)
//             {
//                 condition_is_active = CheckPointsInside<2>(conditions_points_slave, geom_1, geom_2, ContactNormal1, ContactNormal2, ActiveCheckLength, StrictCheck);
//             }
//         }
//         else if (dimension == 3 && number_of_nodes == 3)
//         {
//             ExactMortarIntegrationUtility<3, 3> IntUtil = ExactMortarIntegrationUtility<3, 3>();
//             std::vector<array_1d<PointType,3>> conditions_points_slave;
//             condition_is_active = IntUtil.GetExactIntegration(geom_1, ContactNormal1, geom_2, ContactNormal2, conditions_points_slave);
//             
//             if (condition_is_active == true)
//             {
//                 condition_is_active = CheckPointsInside<3>(conditions_points_slave, geom_1, geom_2, ContactNormal1, ContactNormal2, ActiveCheckLength, StrictCheck);
//             }
//         }
//         else if (dimension == 3 && number_of_nodes == 4)
//         {
//             ExactMortarIntegrationUtility<3, 4> IntUtil = ExactMortarIntegrationUtility<3, 4>();
//             std::vector<array_1d<PointType,3>> conditions_points_slave;
//             condition_is_active = IntUtil.GetExactIntegration(geom_1, ContactNormal1, geom_2, ContactNormal2, conditions_points_slave);
//             
//             if (condition_is_active == true)
//             {
//                 condition_is_active = CheckPointsInside<3>(conditions_points_slave, geom_1, geom_2, ContactNormal1, ContactNormal2, ActiveCheckLength, StrictCheck);
//             }
//         }
//         else
//         {
//             KRATOS_ERROR << "INTEGRATION NOT IMPLEMENTED" << std::endl;
//         }
//         
//         // If condition is active we add
//         if (condition_is_active == true)
//         {
//             for (unsigned int i_node = 0; i_node < geom_1.size(); i_node++)
//             {
//                 geom_1[i_node].Set(ACTIVE, true);
//             }
//             
//             ConditionPointers->AddNewCondition(pCond2);
//         }
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
    
    template<unsigned int TDim>
    static inline bool CheckPointsInside(
        std::vector<array_1d<PointType,TDim>> ConditionsPointsSlave,
        GeometryType& Geom1, // SLAVE
        GeometryType& Geom2,  // MASTER
        const array_1d<double, 3> & ContactNormal1, // SLAVE
        const array_1d<double, 3> & ContactNormal2, // MASTER
        const double ActiveCheckLength,
        const bool StrictCheck = true
        )
    {        
        // Define the basic information
        Point<3> projected_point;
        double aux_distance = 0.0;
        Vector N;
        
        for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); i_geom++)
        {
            for (unsigned int i_node = 0; i_node < TDim; i_node++)
            {
                const PointType LocalPoint = ConditionsPointsSlave[i_geom][i_node];
                PointType GlobalPoint;
                Geom1.GlobalCoordinates(GlobalPoint, LocalPoint);
                    
                Geom1.ShapeFunctionsValues( N, LocalPoint );
                const array_1d<double, 3> normal = ContactUtilities::GaussPointNormal(N, Geom1);
                
                if (norm_2(normal) < std::numeric_limits<double>::epsilon())
                {
                    aux_distance = ContactUtilities::FastProjectDirection(Geom2, GlobalPoint, projected_point, ContactNormal2, ContactNormal1);
                }
                else
                {
                    aux_distance = ContactUtilities::FastProjectDirection(Geom2, GlobalPoint, projected_point, ContactNormal2, normal);
                }  
                
                array_1d<double, 3> result;
                if (aux_distance <= ActiveCheckLength && (StrictCheck == true ? Geom2.IsInside(projected_point, result, std::numeric_limits<double>::epsilon()) : true)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                {
                    return true;
                }
            }
        }
        
        return false;
    }

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
