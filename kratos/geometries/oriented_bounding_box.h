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
#include <iomanip>

// External includes

// Project includes
#include "geometries/point.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"

namespace Kratos {


///@name  Enums
///@{

    /**
     * @brief This enum defines the different types of checks that can be done for the HasIntersection:
     * - Direct: Checks if the nodes are inside the other OBB and if the faces intersect
     * - SeparatingAxisTheorem: This method is more efficient. The Separating Axis Theorem (SAT for short) essentially states if you are able to draw a line to separate two polygons, then they do not collide. It's that simple.
     * @details https://gamedevelopment.tutsplus.com/tutorials/collision-detection-using-the-separating-axis-theorem--gamedev-169
     */
    enum class OBBHasIntersectionType {Direct = 0, SeparatingAxisTheorem = 1};

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
 * For more details
 *      - https://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf
 *      - https://gamedev.stackexchange.com/questions/112883/simple-3d-obb-collision-directx9-c
 * @tparam TDim The dimension of the space
 * @author Vicente Mataix Ferrandiz
 */
template<std::size_t TDim>
class KRATOS_API(KRATOS_CORE) OrientedBoundingBox
{
public:

    ///@name Type Definitions
    ///@{

    // Node definition
    typedef Node                NodeType;

    /// Definition of geometries
    typedef Geometry<NodeType> GeometryType;

    /// Index type definition
    typedef std::size_t           IndexType;

    /// Size type definition
    typedef std::size_t            SizeType;

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

    /**
     * @brief Default constructors
     * @param rCenterPoint The center of the OrientedBoundingBox
     * @param rAxisCoordinates The coordinates that define the axis (orientation vectors and half lengths)
     */
    OrientedBoundingBox(
        const array_1d<double, 3>& rCenterCoords,
        const array_1d<array_1d<double, 3>, TDim>& rAxisCoordinates
        );

    /**
     * @brief Default constructors (with geometry)
     * @param rGeometry The geometry to be considered to build a OBB
     * @param BoundingBoxFactor The bounding box factor
     */
    OrientedBoundingBox(
        const GeometryType& rGeometry,
        const double BoundingBoxFactor,
        const bool BuildFromBoundingBox = true
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
     * @param rOtherOrientedBoundingBox The other Oriented Bounding Box considered
     * @param OBBType The OBB intersection type considered
     */
    bool HasIntersection(
        const OrientedBoundingBox& rOtherOrientedBoundingBox,
        const OBBHasIntersectionType OBBType = OBBHasIntersectionType::SeparatingAxisTheorem
        ) const;

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

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream string_out_coordinates;
        for (std::size_t i = 0; i < TDim; ++i) {
            string_out_coordinates
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(3)
            << std::uppercase
            << "\t" << mPointCenter[i];
        }

        std::stringstream string_out_orientation_vectors;
        for (std::size_t i = 0; i < TDim; ++i) {
            string_out_orientation_vectors << "\nThe orientation axis " << i << " is: ";
            for (std::size_t j = 0; j < TDim; ++j) {
                string_out_orientation_vectors
                << std::setiosflags(std::ios::scientific)
                << std::setprecision(3)
                << std::uppercase
                << "\t" << mOrientationVectors[i][j];
            }
        }

        std::stringstream string_out_half_length;
        for (std::size_t i = 0; i < TDim; ++i) {
            string_out_half_length
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(3)
            << std::uppercase
            << "\t" << mHalfLength[i];
        }

        return "OrientedBoundingBox in " + std::to_string(TDim) + "D space" + "\nWhich center is:" + string_out_coordinates.str() + "\nThe orientation axis are: " + string_out_orientation_vectors.str() + "\nThe half lengths are: " + string_out_half_length.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    array_1d<double, 3> mPointCenter;                        /// This defines the center point of the box
    array_1d<array_1d<double, 3>, TDim> mOrientationVectors; /// This defines the orientation vectors of the OrientedBoundingBox
    array_1d<double, TDim> mHalfLength;                      /// This defines the half of the distance which defines the walls of the box of the OrientedBoundingBox

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Computes the intersection between two OrientedBoundingBox (current and new)
     * @details See https://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf
     * @param rOtherOrientedBoundingBox The other Oriented Bounding Box considered
     */
    bool DirectHasIntersection(const OrientedBoundingBox& rOtherOrientedBoundingBox) const;

    /**
     * @brief Computes the intersection between two OrientedBoundingBox (current and new)
     * @details See: https://gamedev.stackexchange.com/questions/112883/simple-3d-obb-collision-directx9-c
     * @param rOtherOrientedBoundingBox The other Oriented Bounding Box considered
     */
    bool SeparatingAxisTheoremHasIntersection(const OrientedBoundingBox& rOtherOrientedBoundingBox) const;

    /**
     * @brief Check if there's a separating plane in between the selected axis
     * @param rRelativePosition The difference between the center of the both boxes
     * @param rPlane The vector which defines the plane
     * @param rOtherOrientedBoundingBox The other Oriented Bounding Box considered
     */
    bool GetSeparatingPlane(
        const array_1d<double, 3>& rRelativePosition,
        const array_1d<double, 3>& rPlane,
        const OrientedBoundingBox& rOtherOrientedBoundingBox
        ) const;

    /**
     * @brief Check if there's a separating plane in between the selected axis (2D version)
     * @param rRelativePosition The difference between the center of the both boxes
     * @param rPlane The vector which defines the plane
     * @param rOtherOrientedBoundingBox The other Oriented Bounding Box considered
     */
    bool GetSeparatingPlane2D(
        const array_1d<double, 3>& rRelativePosition,
        const array_1d<double, 3>& rPlane,
        const OrientedBoundingBox& rOtherOrientedBoundingBox
        ) const;

    /**
     * @brief Check if there's a separating plane in between the selected axis (3D version)
     * @param rRelativePosition The difference between the center of the both boxes
     * @param rPlane The vector which defines the plane
     * @param rOtherOrientedBoundingBox The other Oriented Bounding Box considered
     */
    bool GetSeparatingPlane3D(
        const array_1d<double, 3>& rRelativePosition,
        const array_1d<double, 3>& rPlane,
        const OrientedBoundingBox& rOtherOrientedBoundingBox
        ) const;

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
        const BoundedMatrix<double, 4, 4>& rInvertedRotationMatrix
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
}; // Class OrientedBoundingBox

///@}

}  // namespace Kratos.

#endif // OBB_CLASS_H_DEFINED  defined
