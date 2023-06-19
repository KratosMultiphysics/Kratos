//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Dagmawi Bekel
//                   Ruben Zorrilla
//

#if USE_TRIANGLE_NONFREE_TPL

#if !defined(KRATOS_CAD_TESSELLATION_MODELER_INCLUDED)
#define KRATOS_CAD_TESSELLATION_MODELER_INCLUDED

// System includes

// External includes

// Project includes
#include "modeler.h"
#include "utilities/tessellation_utilities/curve_tessellation.h"
#include "geometries/brep_surface.h"
#include "geometries/brep_curve_on_surface.h"
#include "integration/triangle_gauss_legendre_integration_points.h"

namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class KRATOS_API(KRATOS_CORE) CadTessellationModeler : public Modeler {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CadTessellationModeler
    KRATOS_CLASS_POINTER_DEFINITION(CadTessellationModeler);

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Node NodeType;

    typedef Point EmbeddedNodeType;

    typedef PointerVector<NodeType> ContainerNodeType;

    typedef PointerVector<EmbeddedNodeType> ContainerEmbeddedNodeType;

    typedef BrepSurface<ContainerNodeType, ContainerEmbeddedNodeType> BrepSurfaceType;

    typedef BrepCurveOnSurface<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CadTessellationModeler() : Modeler()
    {
    }

    /**
     * @brief Constructor using a Model and Parameters
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the discretization
     */
    CadTessellationModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~CadTessellationModeler() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates the Modeler Pointer and returns a pointer to a new
     * CadTessellationModeler, created using the given input
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the discretization
     * @return a Pointer to the new Modeler
     */
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters) const override;

    ///@}
    ///@name Stages
    ///@{

    /**
     * @brief Convert the geometry into an analysis suitable ModePart
     */
    void SetupModelPart() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CadTessellationModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;

    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the tessellation of the boundary curves in the
     * geometric space and maps the points into the parametric space of the surface
     * @param rBoundarySegment Segment of the boundary, which is supposed to be tessellated
     * @return tesselled segment in the parametric space of the surface
     */
    std::vector<array_1d<double, 2>> ComputeBoundaryTessellation(const BrepCurveOnSurfaceType& rBoundarySegment);

    /**
     * @brief This method computes the triangulation of the surface in its parametric space
     * @param rSurfaceGeometry reference of the surface, which is supposed to be triangulated
     * @param rBoundaryLoop tessellated closed polygon of the outer boundaries
     * @return a vector of the triangles
     */
    std::vector<BoundedMatrix<double,3,3>> ComputeSurfaceTriangulation(
        const BrepSurfaceType& rSurfaceGeometry,
        const std::vector<array_1d<double, 2>>& rBoundaryLoop);

    /**
     * @brief This method inserts Gauss points into the surface in the parametric space
     * and projects these points onto the exact input surface. Hence, these points lie within
     * the exact surface.
     * @param rSurfaceGeometry reference of the surface on which the points are supposed to be maped
     * @param rPointsCoordinates reference of the list of points coordinates from the triangularization
     * @param rTriangleConnectivities reference of the list of triangles connectivities from the triangularization
     * @return a vector containing each triangle Gauss points mapped onto the exact surface
     */
    std::vector<BoundedMatrix<double,3,3>> InsertGaussPointsExactSurface(
        const BrepSurfaceType& rSurfaceGeometry,
        const std::vector<double>& rPointsCoordinates,
        const std::vector<IndexType>& rTriangleConnectivities);

    /**
     * @brief This method insert maps the triangulation from the parametric into the physical space.
     * Subsequently, Gauss points are inserted into the approximative surface. Hence, the Gauss points
     * lie within the discretization of the surface and not the exact surface
     * @see InsertGaussPointsExactSurface
     * @param rSurfaceGeometry reference of the surface
     * @param rPointsCoordinates reference of the list of points coordinates from the triangularization
     * @param rTriangleConnectivities reference of the list of triangles connectivities from the triangularization
     * @return a vector of each triangle Gauss points maped onto the triangulation of the surface
     */
    std::vector<BoundedMatrix<double,3,3>> InsertGaussPointsApproxSurface(
        const BrepSurfaceType& rSurfaceGeometry,
        const std::vector<double>& rPointsCoordinates,
        const std::vector<IndexType>& rTriangleConnectivities);

    /**
     * @brief This method computes the discretization error of the surface, measured at the Gauss points.
     * This is done by measuring the distance between distinctive points (Gauss points) between the exact
     * surface and the approaximative surface.
     * @see InsertGaussPointsExactSurface
     * @see InsertGaussPointsApproxSurface
     * @param rGaussPointsExact reference of distincitive points in the exact surface
     * @param rGaussPointsApprox reference of distincitive points in the discretization
     * @return the maximum Gauss point error
     */
    double ComputeDiscretizationError(
        const std::vector<BoundedMatrix<double,3,3>>& rGaussPointsExact,
        const std::vector<BoundedMatrix<double,3,3>>& rGaussPointsApprox);

}; // Class CadTessellationModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream,
    CadTessellationModeler& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const CadTessellationModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_CAD_TESSELLATION_MODELER_INCLUDED  defined
#endif