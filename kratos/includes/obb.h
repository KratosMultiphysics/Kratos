//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 

#if !defined(OBB_H_DEFINED )
#define  OBB_H_DEFINED

// System includes

// External includes

// Project includes
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

/**
 * @class OBB
 * @ingroup KratosCore
 * @brief This class defines the Oriented bounding box class
 * @details The geometrical definition of the OBB can be done as the
 *                      ` *
 *             direction / `
 *                  *   / half diagonal
 *                  `  *C  *
 *                  ` /   `
 *                  `/  `
 *                  * `
 * For more details https://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf
 * @tparam
 * @author Vicente Mataix Ferrandiz
 */
template<std::size_t TDim>
class OBB
    : public Point
{
public:

    ///@name Type Definitions
    ///@{
    
    typedef Point BaseType; 
    
    /// Counted pointer of OBB
    KRATOS_CLASS_POINTER_DEFINITION( OBB );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructors with point
     * @param rCenterPoint The center of the OBB
     * @param rOrientationVector The orientation vector of the diagonal
     * @param HalfDiagonal The half length of the diagonal
     */
    OBB(
        const BaseType& rCenterPoint,
        const array_1d<double, 3>& rOrientationVector,
        const double HalfDiagonal
        ): BaseType(rPoint),
           mOrientationVector(OrientationVector),
           mHalfDiagonal(HalfDiagonal)
    {}
    
    /**
     * @brief Default constructors with point
     * @param rCenterPoint The center of the OBB
     * @param rOrientationVector The orientation vector of the diagonal
     * @param HalfDiagonal The half length of the diagonal
     */
    OBB(
        const array_1d<double, 3>& rCenterCoords,
        const array_1d<double, 3>& rOrientationVector,
        const double HalfDiagonal
        ): BaseType(rCenterCoords),
           mOrientationVector(rOrientationVector),
           mHalfDiagonal(HalfDiagonal)
    {}

    ///Copy constructor  (not really required)
    OBB(const OBB& rhs):
        BaseType(rhs),
        mOrientationVector(rhs.mOrientationVector)
        mHalfDiagonal(rhs.mHalfDiagonal)
    {
    }

    /// Destructor.
    ~OBB() override= default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the point that defines the center of the OBB
     * @return The center point of the OBB
     */
    BaseType GetCenter()
    {
        BaseType Point(this->Coordinates());

        return Point;
    }

    /**
     * @brief Set the point that defines the center of the OBB
     * @param rCenterPoint The point that defines the center of the OBB
     */
    void SetCenter(const BaseType& rCenterPoint)
    {
        this->Coordinates() = rCenterPoint.Coordinates();
    }

    /**
     * @brief Set the point that defines the center of the OBB
     * @param rCenterCoords The coordinates that defines the center of the OBB
     */
    void SetCenter(const array_1d<double, 3>& rCenterCoords)
    {
        this->Coordinates() = rCenterCoords;
    }

    /**
     * @brief Returns the vector that defines the orientation of the axis
     * @return The orientation vector
     */
    array_1d<double, 3>& GetOrientationVector()
    {
        return mOrientationVector;
    }
    
    /**
     * @brief Set the vector that defines the orientation of the axis
     * @param rOrientationVector The orientation vector
     */
    void SetOrientationVector(const array_1d<double, 3>& rOrientationVector)
    {
        noalias(OrientationVector) = rOrientationVector;
    }

    /**
     * @brief Returns the length of the half of the diagonal
     * @return The length of the half of the diagonal
     */
    double& GetHalfDiagonal()
    {
        return mHalfDiagonal;
    }

    /**
     * @brief Set the length of the half of the diagonal
     * @param HalfDiagonal The length of the half of the diagonal
     */
    void SetHalfDiagonal(const double HalfDiagonal)
    {
        mHalfDiagonal = HalfDiagonal;
    }

    /**
     * @brief Computes the intersection between two OBB (current and new)
     */
    bool HasIntersection(const OBB& rOtherOBB)
    {
        // TODO Finish
        return false;
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

    array_1d<double, 3> mOrientationVector; /// This defines the orientation vector of the OBB
    double mHalfDiagonal;                   /// This defines the half of the distance between the highest point and the lowest point

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
}; // Class OBB

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // OBB_H_DEFINED  defined
