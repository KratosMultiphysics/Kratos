// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 

#if !defined(POINT_ITEM_DEFINED )
#define  POINT_ITEM_DEFINED

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "geometries/point.h"

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
class PointItem: public Point
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
        Point()
    {}

    PointItem(const array_1d<double, 3> Coords):
        Point(Coords)
    {}
    
    PointItem(Condition::Pointer Cond):
        mpOriginCond(Cond)
    {
        UpdatePoint(0.0);
    }
    
    PointItem(
        const array_1d<double, 3> Coords,
        Condition::Pointer Cond
    ):
        Point(Coords),
        mpOriginCond(Cond)
    {}

    ///Copy constructor  (not really required)
    PointItem(const PointItem& rhs):
        Point(rhs),
        mpOriginCond(rhs.mpOriginCond)
    {
    }

    /// Destructor.
    ~PointItem() override= default;

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
    Point GetPoint()
    {
        Point Point(this->Coordinates());
        
        return Point;
    }
    
    /**
     * Set the point
     * @param The point
     */
    void SetPoint(const Point Point)
    {
        this->Coordinates() = Point.Coordinates();
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
     * This function updates the database, using as base for the coordinates the condition center
     * @return Coordinates: The coordinates of the item
     */

    void UpdatePoint(const double& DeltaTime)
    {        
        Point center;
        if (mpOriginCond->GetGeometry()[0].SolutionStepsDataHas(VELOCITY_X) == true)
        {
            center = ContactUtilities::GetHalfJumpCenter(mpOriginCond->GetGeometry(), DeltaTime); // NOTE: Center in half delta time
        }
        else
        {
            center = mpOriginCond->GetGeometry().Center(); // NOTE: Real center
        }
        
        this->Coordinates() = center.Coordinates();
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
