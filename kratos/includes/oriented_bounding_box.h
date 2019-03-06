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

#if !defined(OBB_CLASS_H_DEFINED )
#define  OBB_CLASS_H_DEFINED

// System includes

// External includes

// Project includes
#include "geometries/point.h"
#include "geometries/geometry.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"

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
 * @class OrientedBoundingBox
 * @ingroup KratosCore
 * @brief This class defines the Oriented bounding box class
 * @details The geometrical definition of the OrientedBoundingBox can be done as the
 *                      ` *
 *            directions / `
 *                  * \ / half diagonals
 *                  `  *C  *
 *                  ` / \ `
 *                  `/  `
 *                  * `
 * For more details https://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf
 * @tparam TDim The dimension of the space
 * @author Vicente Mataix Ferrandiz
 */
template<std::size_t TDim>
class OrientedBoundingBox
{
public:

    ///@name Type Definitions
    ///@{

    /// Zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Definition of the output type
    typedef typename std::conditional<TDim == 2, Quadrilateral2D4<Point>, Hexahedra3D8<Point> >::type  OutputType;

    /// Counted pointer of OrientedBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( OrientedBoundingBox );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructors
     * @param rCenterPoint The center of the OrientedBoundingBox
     * @param rOrientationVectors The orientation vector of the diagonal
     * @param rHalfLength The half sides
     */
    OrientedBoundingBox(
        const array_1d<double, 3>& rCenterCoords,
        const array_1d<array_1d<double, 3>, TDim>& rOrientationVectors,
        const array_1d<double, TDim>& rHalfLength
        );

    ///Copy constructor  (not really required)
    OrientedBoundingBox(const OrientedBoundingBox& rhs):
        mPointCenter(rhs.mPointCenter),
        mOrientationVectors(rhs.mOrientationVectors),
        mHalfLength(rhs.mHalfLength)
    {
    }

    /// Destructor.
    ~OrientedBoundingBox()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the point that defines the center of the OrientedBoundingBox
     * @return The center point of the OrientedBoundingBox
     */
    const array_1d<double, 3>& GetCenter() const;

    /**
     * @brief Set the point that defines the center of the OrientedBoundingBox
     * @param rCenterCoords The coordinates that defines the center of the OrientedBoundingBox
     */
    void SetCenter(const array_1d<double, 3>& rCenterCoords);

    /**
     * @brief Returns the vector that defines the orientation of the axis
     * @return The orientation vector
     */
    const array_1d<array_1d<double, 3>, TDim>& GetOrientationVectors() const;

    /**
     * @brief Set the vector that defines the orientation of the axis
     * @param rOrientationVectors The orientation vector
     */
    void SetOrientationVectors(const array_1d<array_1d<double, 3>, TDim>& rOrientationVectors);

    /**
     * @brief Returns the length of the half of the diagonal
     * @return The length of the half of the diagonal
     */
    const array_1d<double, TDim>& GetHalfLength() const;

    /**
     * @brief Set the length of the half of the diagonal
     * @param rHalfLength The length of the half of the diagonal
     */
    void SetHalfLength(const array_1d<double, TDim>& rHalfLength);

    /**
     * @brief Computes the intersection between two OrientedBoundingBox (current and new)
     */
    bool IsInside(const OrientedBoundingBox& rOtherOrientedBoundingBox) const;

    /**
     * @brief Computes the intersection between two OrientedBoundingBox (current and new)
     */
    bool HasIntersection(const OrientedBoundingBox& rOtherOrientedBoundingBox) const;

    /**
     * @brief This method egnerates an equiavelent geometry (debugging)
     * @return Getting the OrientedBoundingBox geometry
     */
    OutputType GetEquivalentGeometry() const;

    /**
     * @brief This method egnerates an equiavelent geometry (debugging)
     * @param rGeometry Geometry to rotate
     */
    void GetEquivalentRotatedGeometry(OutputType& rGeometry);

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

    array_1d<double, 3> mPointCenter;                        /// This defines the center point of the box
    array_1d<array_1d<double, 3>, TDim> mOrientationVectors; /// This defines the orientation vectors of the OrientedBoundingBox
    array_1d<double, TDim> mHalfLength;                      /// This defines the half of the distance which defines the walls of the box of the OrientedBoundingBox

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method does a 2D rotation of a point
     * @param rCoords The coordinates of the point of interest
     */
    void RotateNode2D(array_1d<double, 3>& rCoords) const;

    /**
     * @brief This method does a 3D rotation of a point
     * @param rCoords The coordinates of the point of interest
     * @param rInvertedRotationMatrix The inverted matrix of rotation
     */
    void RotateNode3D(
        array_1d<double, 3>& rCoords,
        BoundedMatrix<double, 4, 4> rInvertedRotationMatrix
        ) const;

    /**
     * @brief This method does a check in 2D if the point is inside the OrientedBoundingBox
     * @param rCoords The coordinates of the point of interest
     * @return True is is inside, false otherwise
     */
    bool CheckIsInside2D(array_1d<double, 3>& rCoords) const;

    /**
     * @brief This method does a check in 3D if the point is inside the OrientedBoundingBox
     * @param rCoords The coordinates of the point of interest
     * @param rInvertedRotationMatrix The inverted matrix of rotation
     * @return True is is inside, false otherwise
     */
    bool CheckIsInside3D(
        array_1d<double, 3>& rCoords,
        BoundedMatrix<double, 4, 4> rInvertedRotationMatrix
        ) const;

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
}; // Class OrientedBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // OBB_CLASS_H_DEFINED  defined
