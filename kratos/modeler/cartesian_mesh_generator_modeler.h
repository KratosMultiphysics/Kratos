//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/geometry.h"
#include "containers/pointer_vector.h"
#include "modeler/modeler.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class CartesianMeshGeneratorModeler
 * @ingroup KratosCore
 * @brief Generates a structured Cartesian mesh filling the bounding box of a
 *        source model part.
 * @details For 3D problems each axis-aligned hexahedral cell is split into
 *          6 tetrahedra using the Freudenthal (Kuhn) decomposition so that
 *          the resulting tetrahedral mesh is conforming across cell boundaries.
 *          For 2D problems quadrilateral elements are generated using a
 *          ray-casting inside test derived from the source boundary elements.
 *
 *          The class is constructed with a reference source model part
 *          (from which the bounding box and boundary are taken) and a uniform
 *          element size.  The generated nodes and elements are placed in the
 *          destination model part passed to GenerateMesh().
 *
 * @author Jordi Cotela Dalmau
 */
class KRATOS_API(KRATOS_CORE) CartesianMeshGeneratorModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CartesianMeshGeneratorModeler
    KRATOS_CLASS_POINTER_DEFINITION(CartesianMeshGeneratorModeler);

    /// The base class type
    using BaseType = Modeler;

    /// The geometry type definition
    using GeometryType = Geometry<Node>;

    /// The nodes vector type definition
    using NodesVectorType = PointerVector<Node>;

    /// The size type definition
    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor (required for registry).
     * @note  Do not call GenerateMesh() on a default-constructed instance.
     */
    CartesianMeshGeneratorModeler();

    /**
     * @brief Constructor for use via the modeler factory / registry.
     * @details Reads the input model part name, output model part name,
     *          element name and element size from @p ModelerParameters.
     *          The default parameters are:
     *          @code
     *          {
     *              "input_model_part_name"  : "",
     *              "output_model_part_name" : "",
     *              "element_name"           : "Element3D4N",
     *              "element_size"           : 1.0
     *          }
     *          @endcode
     *          SetupModelPart() performs the actual mesh generation.
     * @param rModel          Model that owns both the input and output model parts.
     * @param ModelerParameters Parameters block (see default values above).
     */
    CartesianMeshGeneratorModeler(Model& rModel, Parameters ModelerParameters = Parameters());

    /**
     * @brief Legacy direct constructor.
     * @param rSourceModelPart Model part whose geometry defines the bounding
     *        box (and, for 2D, the boundary) used during mesh generation.
     * @param ElementSize Uniform edge length of each voxel cell. The actual
     *        number of cells per direction is computed as
     *        ceil(bounding_box_extent / ElementSize).
     */
    CartesianMeshGeneratorModeler(ModelPart& rSourceModelPart, double ElementSize)
        : Modeler()
        , mpModel(nullptr)
        , mpSourceModelPart(&rSourceModelPart)
        , mElementSize(ElementSize) {}

    /// Destructor.
    ~CartesianMeshGeneratorModeler() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new instance via the modeler factory.
     * @param rModel          Owning model.
     * @param ModelParameters Configuration parameters.
     * @return Shared pointer to the new modeler.
     */
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<CartesianMeshGeneratorModeler>(rModel, ModelParameters);
    }

    /**
     * @brief Returns the default parameters accepted by the Model constructor.
     * @return A Parameters object with default values.
     */
    const Parameters GetDefaultParameters() const override;

    /**
     * @brief Runs mesh generation when the modeler is used through the factory.
     * @details Reads the input/output model part names and element name from
     *          the stored parameters and delegates to GenerateMesh().
     */
    void SetupModelPart() override;

    /**
     * @brief Generates a Cartesian mesh in the destination model part.
     * @details The element type (and therefore the dimension) is determined
     *          by @p rReferenceElement.  For 3D elements six tetrahedra are
     *          created per hex cell (Freudenthal decomposition).  For 2D
     *          elements one quadrilateral is created per active cell.
     * @param rThisModelPart Destination model part that will receive the
     *        generated nodes and elements.
     * @param rReferenceElement Prototype element whose geometry dimension
     *        selects the 2D or 3D code path and whose properties are
     *        reused for every created element.
     */
    void GenerateMesh(ModelPart& rThisModelPart, Element const& rReferenceElement);

    /**
     * @brief Computes outward-pointing normals at every boundary node of
     *        the source model part (2D only).
     * @details The normals are scaled to half the element size and
     *          accumulated over all boundary elements sharing each node.
     *          The result is stored in the NORMAL solution-step variable
     *          and in the internal @c mNormals cache used by
     *          FindNearestNodeIndex().
     */
    void CalculateNormals();

    /**
     * @brief Returns the flat index of the grid node closest to @p rThisPoint
     *        in the direction given by @p rNormal (2D only).
     * @details Starting from the cell that contains the point the method
     *          steps one cell in the dominant normal direction until it
     *          finds a node that is marked as inside.
     * @param rThisPoint Query point in world coordinates.
     * @param rNormal    Outward boundary normal at the query point; used to
     *                   determine which neighbouring cell to prefer.
     * @return Flat index into the internal node array of the nearest inside
     *         node.
     */
    unsigned int FindNearestNodeIndex(Point& rThisPoint, array_1d<double,3>& rNormal);

    /**
     * @brief Marks each grid node as inside or outside the 2D boundary.
     * @details Uses the precomputed horizontal ray–edge intersection table
     *          (populated by CalculateBoundaryIntersections()) to flood-fill
     *          the grid cells between pairs of intersections on each
     *          horizontal scan line.
     * @param rThisModelPart The source model part providing the boundary.
     */
    void CalculateIsInside(ModelPart& rThisModelPart);

    /**
     * @brief Builds the horizontal ray–edge intersection table (2D only).
     * @details For every boundary element in @p rThisModelPart the method
     *          computes the x-coordinates at which the element edge crosses
     *          each horizontal grid line and stores them in @c mIntersections.
     *          Each row of the table is sorted in ascending order after all
     *          edges have been processed.
     * @param rThisModelPart The source model part providing the boundary edges.
     */
    void CalculateBoundaryIntersections(ModelPart& rThisModelPart);

    /**
     * @brief Computes the axis-aligned bounding box of @p rThisModelPart.
     * @details When the model part contains elements the bounding box is
     *          computed from element vertices.  When it contains only nodes
     *          (e.g. a condition-only mesh read from an STL file) the bounding
     *          box is computed directly from the node coordinates.  An empty
     *          model part yields the origin for both corners.
     * @param rThisModelPart The model part to bound.
     * @param rMinPoint      Output: minimum corner of the bounding box.
     * @param rMaxPoint      Output: maximum corner of the bounding box.
     */
    void CalculateBoundingBox(ModelPart& rThisModelPart, Point& rMinPoint, Point& rMaxPoint);

    /**
     * @brief Computes the number of grid cells and node layers per direction.
     * @details Sets @c mSegmentsNumber[i] = ceil(extent_i / mElementSize)
     *          and @c mDivisionsNumber[i] = mSegmentsNumber[i] + 1.
     *          Directions with zero extent are given one segment so that the
     *          degenerate case does not crash.
     */
    void CalculateDivisionNumbers();

    /**
     * @brief Returns the flat element index for cell (i, j, k).
     * @param i Cell index along X.
     * @param j Cell index along Y.
     * @param k Cell index along Z.
     * @return Flat index = i + nx*j + nx*ny*k, where nx = mDivisionsNumber[0]
     *         and ny = mDivisionsNumber[1].
     */
    unsigned int ElementIndex(unsigned int i, unsigned int j, unsigned int k)
    {
        return i + mDivisionsNumber[0] * j + mDivisionsNumber[0] * mDivisionsNumber[1] * k;
    }

    /**
     * @brief Returns the flat node index for grid point (i, j, k).
     * @param i Node index along X.
     * @param j Node index along Y.
     * @param k Node index along Z.
     * @return Flat index = i + sx*j + sx*sy*k, where sx = mSegmentsNumber[0]
     *         and sy = mSegmentsNumber[1].
     */
    unsigned int NodeIndex(unsigned int i, unsigned int j, unsigned int k)
    {
        return i + mSegmentsNumber[0] * j + mSegmentsNumber[0] * mSegmentsNumber[1] * k;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// @brief Returns the class name as a string.
    std::string Info() const override;

    /// @brief Prints the class name to @p rOStream.
    void PrintInfo(std::ostream& rOStream) const override;

    /// @brief Prints object data to @p rOStream (currently empty).
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// Owning model (only set when using the Model constructor).
    Model* mpModel = nullptr;

    /// Source model part whose geometry defines the domain to mesh (nullptr only for the default-constructed prototype).
    ModelPart* mpSourceModelPart;

    /// Uniform voxel cell edge length.
    double mElementSize = 0.0;

    /// Minimum corner of the bounding box of the source model part.
    Point mMinPoint;

    /// Maximum corner of the bounding box of the source model part.
    Point mMaxPoint;

    /// Number of voxel cells per direction (X=0, Y=1, Z=2).
    unsigned int mSegmentsNumber[3] = {0, 0, 0};

    /// Number of node layers per direction: mSegmentsNumber[i] + 1.
    unsigned int mDivisionsNumber[3] = {0, 0, 0};

    /// Sorted x-coordinates of edge–scanline intersections per row (2D only).
    std::vector<std::vector<double>> mIntersections;

    /// Per-node inside flag populated by CalculateIsInside() (2D only).
    std::vector<int> mIsInside;

    /// Accumulated boundary normals at each source-model-part node (2D only).
    std::vector<array_1d<double,3>> mNormals;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Creates equally-spaced intermediate nodes along a geometry edge.
     * @param rThisModelPart Model part in which to create the nodes.
     * @param rGeometry      Two-node edge geometry.
     * @param NumberOfSegments Number of segments to divide the edge into;
     *        (NumberOfSegments - 1) new nodes are created (endpoints excluded).
     * @param StartNodeId    First node ID to assign; incremented for each new node.
     */
    void GenerateNodes(ModelPart& rThisModelPart, GeometryType& rGeometry,
                       SizeType NumberOfSegments, SizeType StartNodeId);

    ///@}
    ///@name Inaccessible methods
    ///@{

    CartesianMeshGeneratorModeler& operator=(CartesianMeshGeneratorModeler const& rOther) = delete;
    CartesianMeshGeneratorModeler(CartesianMeshGeneratorModeler const& rOther) = delete;

    ///@}

}; // Class CartesianMeshGeneratorModeler

///@}

/// @brief Input stream operator (no-op).
inline std::istream& operator>>(std::istream& rIStream, CartesianMeshGeneratorModeler& rThis) { return rIStream; }

/// @brief Output stream operator — prints Info() followed by PrintData().
inline std::ostream& operator<<(std::ostream& rOStream, const CartesianMeshGeneratorModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

