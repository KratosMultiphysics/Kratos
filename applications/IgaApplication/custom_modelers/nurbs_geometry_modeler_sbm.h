//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//                   Andrea Gorgi
//

#if !defined(KRATOS_NURBS_GEOMETRY_MODELER_SBM_H_INCLUDED )
#define  KRATOS_NURBS_GEOMETRY_MODELER_SBM_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "nurbs_geometry_modeler.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities.h"
#include "geometries/brep_curve_on_surface.h"
#include "utilities/nurbs_utilities/snake_sbm_utilities.h"

namespace Kratos {

class KRATOS_API(IGA_APPLICATION) NurbsGeometryModelerSbm
    : public NurbsGeometryModeler
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( NurbsGeometryModelerSbm );

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceGeometryType;
    typedef typename NurbsSurfaceGeometryType::Pointer NurbsSurfaceGeometryPointerType;

    typedef NurbsVolumeGeometry<PointerVector<NodeType>> NurbsVolumeGeometryType;
    typedef typename NurbsVolumeGeometryType::Pointer NurbsVolumeGeometryPointerType;

    typedef PointerVector<Node> ContainerNodeType;
    typedef PointerVector<Point> ContainerEmbeddedNodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NurbsGeometryModelerSbm()
        : NurbsGeometryModeler() {}

    /// Constructor.
    NurbsGeometryModelerSbm(
        Model & rModel,
        const Parameters ModelerParameters = Parameters())
        : NurbsGeometryModeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~NurbsGeometryModelerSbm() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<NurbsGeometryModelerSbm>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupGeometryModel() override;

    ///@}

protected:

    /**
     * @brief Creates a regular grid composed out of bivariant B-splines.
     * @param PointA Lower point of bounding box.
     * @param PointB Upper point of bounding box.
     * @param Order  Polynomial degree in each direction u,v.
     * @param NumKnotSpans Number of equidistant elements/knot spans in each direction u,v.
     * @note The CP'S are defined as nodes and added to the rModelPart.
     **/
    void CreateAndAddRegularGrid2D( ModelPart& r_model_part, const Point& A_xyz, const Point& B_xyz, const Point& A_uvw, const Point& B_uvw,
        SizeType OrderU, SizeType OrderV, SizeType NumKnotSpansU, SizeType NumKnotSpansV, bool add_surface_to_model_part ) override;

private:

    ///@name Private Member Variables
    ///@{

    Model* mpModel;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Creates a cartesian grid composed out of trivariant B-spline cubes.
     * @param PointA Lower point of bounding box.
     * @param PointB Upper point of bounding box.
     * @param Order  Polynomial degree in each direction u,v,w.
     * @param NumKnotSpans Number of equidistant elements/knot spans in each direction u,v,w.
     * @note The CP'S are defined as nodes and added to the rModelPart.
     **/
    void CreateAndAddRegularGrid3D( ModelPart& r_model_part, const Point& A_xyz, const Point& B_xyz, const Point& A_uvw, const Point& B_uvw,
       SizeType OrderU, SizeType OrderV, SizeType OrderW, SizeType NumKnotSpansU, SizeType NumKnotSpansV, SizeType NumKnotSpansW );

    Parameters ReadParamatersFile(const std::string& rDataFileName) const;   
};

} // End namesapce Kratos
#endif // KRATOS_NURBS_GEOMETRY_MODELER_H_INCLUDED