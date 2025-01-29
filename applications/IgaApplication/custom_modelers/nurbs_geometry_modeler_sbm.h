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

# pragma once

// System includes

// External includes

// Project includes
#include "nurbs_geometry_modeler.h"

namespace Kratos {

class KRATOS_API(IGA_APPLICATION) NurbsGeometryModelerSbm
    : public NurbsGeometryModeler
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( NurbsGeometryModelerSbm );

    using IndexType = std::size_t;
    using SizeType = std::size_t;
    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;
    using GeometryPointerType = GeometryType::Pointer;

    using NurbsSurfaceGeometryType = NurbsSurfaceGeometry<3, PointerVector<NodeType>>;
    using NurbsSurfaceGeometryPointerType = NurbsSurfaceGeometryType::Pointer;

    using NurbsVolumeGeometryType = NurbsVolumeGeometry<PointerVector<NodeType>>;
    using NurbsVolumeGeometryPointerType = NurbsVolumeGeometryType::Pointer;

    using ContainerNodeType = PointerVector<Node>;
    using ContainerEmbeddedNodeType = PointerVector<Point>;

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
    ~NurbsGeometryModelerSbm() = default;

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