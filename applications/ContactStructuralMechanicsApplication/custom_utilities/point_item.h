// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(POINT_ITEM_DEFINED )
#define  POINT_ITEM_DEFINED

// System includes
#include <iostream>
#include <vector>
#include "boost/smart_ptr.hpp"

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"

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

///@}
///@name Kratos Classes
///@{

/** @brief Custom Point container to be used by the mapper
 */
class PointItem: public Point<3>
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of PointItem
    KRATOS_CLASS_POINTER_DEFINITION( PointItem );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointItem():
        Point<3>()
    {
    }

    PointItem(array_1d<double, 3> Coords):
        Point<3>(Coords)
    {}
    
    PointItem(
        array_1d<double, 3> Coords,
        Condition::Pointer Cond,
        double Radius
    ):
        Point<3>(Coords),
        mpOriginCond(Cond),
        mRadius(Radius)
    {}
    
    PointItem(
        array_1d<double, 3> Coords,
        Node<3>::Pointer Node
    ):
        Point<3>(Coords),
        mpOriginNode(Node)
    {}

    ///Copy constructor  (not really required)
    PointItem(const PointItem& rhs):
        Point<3>(rhs),
        mpOriginCond(rhs.mpOriginCond),
        mpOriginNode(rhs.mpOriginNode),
        mRadius(rhs.mRadius)
    {
    }

    /// Destructor.
    // ~PointItem();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Returns the point
     * @return The point
     */
    Point<3> GetPoint()
    {
        Point<3> Point(this->Coordinates());
        
        return Point;
    }
    
    /**
     * Set the point
     * @param The point
     */
    void SetPoint(Point<3> Point)
    {
        this->Coordinates() = Point.Coordinates();
    }
    
    /**
     * Returns the radius of the condition
     * @return mRadius: The radius of the condition
     */
    double GetRadius()
    {
        return mRadius;
    }
    
    /**
     * Sets the radius of the condition
     * @param Radius: The radius of the condition
     */
    void SetRadius(const double& Radius)
    {
        mRadius = Radius;
    }

    /**
     * Sets the condition associated to the point
     * @param Cond: The pointer to the condition
     */

    void SetCondition(Condition::Pointer Cond)
    {
        mpOriginCond = Cond;
    }
    
    /**
     * Returns the condition associated to the point
     * @return mpOriginCond: The pointer to the condition associated to the point
     */

    Condition::Pointer GetCondition()
    {
        return mpOriginCond;
    }
    
    /**
     * Sets the node associated to the point
     * @param Node: The pointer to the node associated to the point
     */

    void SetNode(Node<3>::Pointer Node)
    {
        mpOriginNode = Node;
    }
    
    /**
     * Returns the condition associated to the point
     * @return mpOriginNode: The pointer to the node associated to the point
     */

    Node<3>::Pointer GetNode()
    {
        return mpOriginNode;
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

    Condition::Pointer mpOriginCond; // Condition pointer
    Node<3>::Pointer   mpOriginNode; // Node pointer
    double                  mRadius; // Radius         
    array_1d<double, 3>     mNormal; // Normal vector      

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
}; // Class PointItem 

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // POINT_ITEM_DEFINED  defined
