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
 * @class OBB
 * @ingroup KratosCore
 * @brief This class defines the Oriented bounding box class
 * @details The geometrical definition of the OBB can be done as the
 *                      ` *
 *            directions / `
 *                  * \ / half diagonals
 *                  `  *C  *
 *                  ` / \ `
 *                  `/  `
 *                  * `
 * For more details https://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf
 * @tparam
 * @author Vicente Mataix Ferrandiz
 */
template<std::size_t TDim>
class OBB
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the output type
    typedef typename std::conditional<TDim == 2, Quadrilateral2D4<Point>, Hexahedra3D8<Point> >::type  OutpuType;

    /// Counted pointer of OBB
    KRATOS_CLASS_POINTER_DEFINITION( OBB );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructors
     * @param rCenterPoint The center of the OBB
     * @param rOrientationVectors The orientation vector of the diagonal
     * @param rHalfLength The half sides
     */
    OBB(
        const array_1d<double, 3>& rCenterCoords,
        const array_1d<array_1d<double, 3>, TDim>& rOrientationVectors,
        const array_1d<double, TDim>& rHalfLength
        );

    ///Copy constructor  (not really required)
    OBB(const OBB& rhs):
        mPointCenter(rhs.mPointCenter),
        mOrientationVectors(rhs.mOrientationVectors),
        mHalfLength(rhs.mHalfLength)
    {
    }

    /// Destructor.
    ~OBB()
    {
    }

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
    const array_1d<double, 3>& GetCenter() const;

    /**
     * @brief Set the point that defines the center of the OBB
     * @param rCenterCoords The coordinates that defines the center of the OBB
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
     * @brief Computes the intersection between two OBB (current and new)
     */
    bool IsInside(const OBB& rOtherOBB) const;

    /**
     * @brief Computes the intersection between two OBB (current and new)
     */
    bool HasIntersection(const OBB& rOtherOBB) const;

    /**
     * @brief This method egnerates an equiavelent geometry (debugging)
     */
    OutpuType GetEquiavelentGeometry();

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
    array_1d<array_1d<double, 3>, TDim> mOrientationVectors; /// This defines the orientation vectors of the OBB
    array_1d<double, TDim> mHalfLength;                      /// This defines the half of the distance which defines the walls of the box of the OBB

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
     */
    void RotateNode3D(array_1d<double, 3>& rCoords) const;

    /**
     * @brief This method does a check in 2D if the point is inside the OBB
     * @param rCoords The coordinates of the point of interest
     * @return True is is inside, false otherwise
     */
    bool CheckIsInside2D(array_1d<double, 3>& rCoords) const;

    /**
     * @brief This method does a check in 3D if the point is inside the OBB
     * @param rCoords The coordinates of the point of interest
     * @return True is is inside, false otherwise
     */
    bool CheckIsInside3D(array_1d<double, 3>& rCoords) const;

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

#endif // OBB_CLASS_H_DEFINED  defined
