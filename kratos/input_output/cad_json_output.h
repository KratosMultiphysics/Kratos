//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_CAD_JSON_OUTPUT_INCLUDED )
#define  KRATOS_CAD_JSON_OUTPUT_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Geometries
#include "geometries/coupling_geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"

#include "geometries/brep_surface.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/brep_curve.h"

#include "geometries/thb_surface_geometry.h"
#include "geometries/thb_brep_surface.h"

namespace Kratos
{

///@name Kratos Classes
///@{
/// Output for CAD-files.
/** Gives output capabilities for Nurbs based Brep models in the JSON format defined in
https://amses-journal.springeropen.com/articles/10.1186/s40323-018-0109-4. */
class KRATOS_API(KRATOS_CORE) CadJsonOutput
{
    public:

    ///@}
    ///@name Type Definitions
    ///@}

    /// Pointer definition of CadJsonOutput
    KRATOS_CLASS_POINTER_DEFINITION(CadJsonOutput);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node NodeType;
    typedef Point EmbeddedNodeType;

    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    typedef PointerVector<NodeType> ContainerNodeType;
    typedef PointerVector<EmbeddedNodeType> ContainerEmbeddedNodeType;

    typedef CouplingGeometry<NodeType> CouplingGeometryType;

    typedef NurbsSurfaceGeometry<3, ContainerNodeType> NurbsSurfaceType;
    typedef NurbsCurveGeometry<2, ContainerEmbeddedNodeType> NurbsTrimmingCurveType;

    typedef typename NurbsSurfaceType::Pointer NurbsSurfacePointerType;
    typedef typename NurbsTrimmingCurveType::Pointer NurbsTrimmingCurvePointerType;

    typedef BrepSurface<ContainerNodeType, false, ContainerEmbeddedNodeType> BrepSurfaceType;
    typedef BrepCurveOnSurface<ContainerNodeType, false, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;
    typedef BrepCurve<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveType;

    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceArrayType;
    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceLoopType;
    typedef DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>> BrepCurveOnSurfaceLoopArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with path to input file.
    static void GetCadJsonOutput(
        ModelPart& rModelPart,
        std::string& rGeometryFile,
        SizeType EchoLevel = 0)
    {
        Parameters cad_geometry_parameters;

        GetParameters(rModelPart, cad_geometry_parameters, EchoLevel);

        rGeometryFile = cad_geometry_parameters.PrettyPrintJsonString();
    }

    ///@}
    ///@name Python exposed Functions
    ///@{

    /// Adds all CAD geometries to the herin provided model_part.
    static void GetParameters(
        ModelPart& rModelPart, Parameters& rCadGeometry, IndexType EchoLevel);

    /// Returns the paramaters/json of a specific brep_surface
    static void GetBrepSurfaceParameters(
        const typename ModelPart::GeometryIterator& rGeometryIterator, Parameters& rBrepsParameters, IndexType EchoLevel);

    /// Returns the paramaters/json of a loop - typically from a brep_surface
    static void GetBoundaryLoopParameters(
        const BrepCurveOnSurfaceArrayType& rCurveOnSurfaceArray, Parameters& rCadGeometry, IndexType EchoLevel);

    ///@}
}; // Class CadJsonOutput
}  // namespace Kratos.

#endif // KRATOS_CAD_JSON_OUTPUT_INCLUDED  defined
