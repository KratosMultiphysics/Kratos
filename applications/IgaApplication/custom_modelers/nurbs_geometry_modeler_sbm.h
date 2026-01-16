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
    Modeler::Pointer Create(
        Model& rModel, 
        const Parameters ModelParameters) const override
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
    void CreateAndAddRegularGrid2D(
        ModelPart& rModelPart, 
        const Point& A_xyz, 
        const Point& B_xyz, 
        const Point& A_uvw, 
        const Point& B_uvw,
        const SizeType OrderU, 
        const SizeType OrderV, 
        const SizeType NumKnotSpansU, 
        const SizeType NumKnotSpansV, 
        const bool AddSurfaceToModelPart) override;
    

    /**
     * @brief Create a And Add Regular Grid 3 D object
     * 
     * @param rModelPart 
     * @param A_xyz 
     * @param B_xyz 
     * @param A_uvw 
     * @param B_uvw 
     * @param OrderU 
     * @param OrderV 
     * @param OrderW 
     * @param NumKnotSpansU 
     * @param NumKnotSpansV 
     * @param NumKnotSpansW 
     * @param AddVolumeToModelPart 
     * @return * void 
     */
    void CreateAndAddRegularGrid3D(
        ModelPart& rModelPart, 
        const Point& A_xyz, 
        const Point& B_xyz, 
        const Point& A_uvw, 
        const Point& B_uvw,
        const SizeType OrderU, 
        const SizeType OrderV, 
        const SizeType OrderW, 
        const SizeType NumKnotSpansU, 
        const SizeType NumKnotSpansV, 
        const SizeType NumKnotSpansW, 
        const bool AddVolumeToModelPart) override;

private:

    ///@name Private Member Variables
    ///@{

    Model* mpModel;

    ///@}
    ///@name Private Operations
    ///@{

};

} // End namesapce Kratos