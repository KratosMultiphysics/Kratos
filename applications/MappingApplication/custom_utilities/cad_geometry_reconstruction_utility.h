//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

// Geometries
#include "geometries/coupling_geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/brep_surface.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/brep_curve.h"
#include "geometries/point_on_geometry.h"

namespace Kratos
{
namespace CadGeometryReconstructionUtility
{
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    using NodeType = Node;
    using EmbeddedNodeType = Point;

    using GeometryType = Geometry<NodeType>;
    using GeometryPointerType = typename GeometryType::Pointer;

    using ContainerNodeType = PointerVector<NodeType>;
    using ContainerEmbeddedNodeType = PointerVector<EmbeddedNodeType>;

    using CouplingGeometryType = CouplingGeometry<NodeType>;

    using NurbsSurfaceType = NurbsSurfaceGeometry<3, ContainerNodeType>;
    using NurbsTrimmingCurveType = NurbsCurveGeometry<2, ContainerEmbeddedNodeType>;

    using NurbsSurfacePointerType = typename NurbsSurfaceType::Pointer;
    using NurbsTrimmingCurvePointerType = typename NurbsTrimmingCurveType::Pointer;

    using BrepSurfaceType = BrepSurface<ContainerNodeType, false, ContainerEmbeddedNodeType>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<ContainerNodeType, false, ContainerEmbeddedNodeType>;
    using BrepCurveType = BrepCurve<ContainerNodeType, ContainerEmbeddedNodeType>;
    using PointOnGeometryOnSurfaceType = PointOnGeometry<ContainerNodeType, 3, 2>;
    using PointOnGeometryOnCurveType = PointOnGeometry<ContainerNodeType, 3, 1>;

    using BrepCurveOnSurfaceArrayType = DenseVector<typename BrepCurveOnSurfaceType::Pointer>;
    using BrepCurveOnSurfaceLoopType = DenseVector<typename BrepCurveOnSurfaceType::Pointer>;
    using BrepCurveOnSurfaceLoopArrayType = DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>>;

    /// Reconstructs BrepSurface geometries from a CAD json file using already existing nodes in rModelPart.
    void KRATOS_API(MAPPING_APPLICATION) ReconstructModelPartBrepGeometryFromCadJson(
        const std::string& rDataFileName,
        ModelPart& rModelPart,
        const int EchoLevel = 0);

    /// Same as above, but taking already parsed Parameters.
    void KRATOS_API(MAPPING_APPLICATION) ReconstructModelPartBrepGeometryFromCadJson(
        const Parameters CadJsonParameters,
        ModelPart& rModelPart,
        const int EchoLevel = 0);

    /// Reconstructs BrepSurface geometries from a CAD json file by creating the control-point nodes from the file.
    void KRATOS_API(MAPPING_APPLICATION) ReconstructModelPartBrepGeometryFromCadJsonCreatingNodes(
        const std::string& rDataFileName,
        ModelPart& rModelPart,
        const int EchoLevel = 0);

    /// Same as above, but taking already parsed Parameters.
    void KRATOS_API(MAPPING_APPLICATION) ReconstructModelPartBrepGeometryFromCadJsonCreatingNodes(
        const Parameters CadJsonParameters,
        ModelPart& rModelPart,
        const int EchoLevel = 0);

} // namespace CadGeometryReconstructionUtility
} // namespace Kratos